/********************************************************************************************
 *
 *  Phase 1 of FastK: First the minimal core prefix trie is found over the first 1 Gbp of
 *    the data set and then the entire data set is scanned and partitioned into super-mers
 *    that are sent to file buckets according to the trie.
 *
 *  Author :  Gene Myers
 *  Date   :  October 2020
 *
 ********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <pthread.h>

#include "gene_core.h"
#include "FastK.h"

#define  USE_MAPPING
#undef   HISTOGRAM_TEST
#undef   MAXIMUM_TEST
#undef   FIND_MTHRESH

#undef   DEBUG_PADDED
#define    PAD_LEVEL 0
#undef   DEBUG_SCHEME
#undef   DEBUG_DISTRIBUTE
#undef   SHOW_PACKETS
#define    PACKET 0
#undef   DEBUG_COMPRESS

#define  SKIP_SPLIT
#define  SHOW_DISTRIBUTION

#define THREAD pthread_t

/*******************************************************************************************
 *
 *  Optional homopolymer compressor
 *
 ********************************************************************************************/

static void compress_block(DATA_BLOCK *block)
{ int    nreads = block->nreads;
  int64 *boff   = block->boff+1;
  char  *bases  = block->bases;

  char *s, *t;
  int   q, p, i;
  char *c, x;
  int64 u;
  int   mlen;

  mlen = 0;

  u = 0;
  c = bases;
  t = bases + boff[-1];
  for (i = 0; i < nreads; i++)
    { s = t;
      t = bases + boff[i];

      for (p = 0; p < BC_PREFIX; p++)
        c[u++] = s[p];
      s += BC_PREFIX;
      q = (t-s) - 1;

      c[u++] = x = s[0];
      for (p = 1; p < q; p++)
        if (s[p] != x)
          c[u++] = x = s[p];

      q = u - boff[i-1];
      if (q > mlen)
        mlen = q;

      c[u++] = '\0';
      boff[i] = u;
    }

  block->totlen = u-nreads;
  block->maxlen = mlen;
}

/*******************************************************************************************
 *
 *  DETERMINE PARTITION FOR SUPER K-MER DISTRIBUTION USING FIRST BLOCK
 *     int Determine_Scheme(DATA_BLOCK *block)
 *        Determine base mapping and then the core prefix trie in a series of
 *        scans that measure the tentative core prefix frequencies on the large (1Gbp)
 *        supplied training set in block.  Once the core trie is computed, then
 *        assign prefixes to buckets as evenly as possible with a heuristic.
 *        Return the # of k-mers in the longest super-mer.
 *
 ********************************************************************************************/

#define MIN_LEN      5   // Seed minimizer length
#define MIN_TOT  0x400   // = 4^MIN_LEN

  //  Translation vectors

static IO_UTYPE Dran[256];   //  Dran[x] = 0,1,2,3 for a,c,g,t
static IO_UTYPE Fran[256];   //  fran[x] = 3,2,1,0 for a,c,g,t

static uint64 Tran[256];   //  Tran[x] in 0-3 for a,c,g,t according to frequency  
static uint64 Cran[256];   //  Cran[x] complement of Tran[x] shifted by MIN_L1 or later PAD_L1

  //  Padded Minimizer scheme

static int  PAD;      //  # of base pairs beyond MIN_LEN
static int  PAD2;     //  2*PAD

static int     PAD_LEN;  // MIN_LEN + PAD
static int     PAD_L1;       // = PAD_LEN-1
static int64   PAD_TOT;      // = 4^PAD_LEN
static uint64  PAD_MSK;      // = PAD_TOT-1

static int  Min_States;
static int *Min_Part;      //  The core prefix trie: MP[x] < 0 => goto a-MP[x] on base a
                           //      MP[x] >= 0 => bucket for prefix

typedef struct
  { DATA_BLOCK *block;
    int         beg;
    int         end;
    int64       freq[256];
    int64       nmin;
    int64      *count;
  } Partition_Arg;


  //  Thread to sample the frequency of bases (needed for base mapping Tran/Cran

static void *frequency_thread(void *arg)
{ Partition_Arg *data   = (Partition_Arg *) arg;
  DATA_BLOCK    *block  = data->block;
  int64         *freq   = data->freq;

  int        p, q;
  char      *s;

  for (p = 0; p < 256; p++)
    freq[p] = 0;

  s = block->bases + block->boff[data->beg];
  q = block->boff[data->end] - block->boff[data->beg];
  for (p = 0; p < q; p++)
    freq[(int) s[p]] += 1;                 //  Don't care if a few \0 counts are recorded

  return (NULL);
}

  // for reads [beg,end) compute # of kmers for each core minimizer prefix in proposed scheme

static void *padded_minimizer_thread(void *arg)
{ Partition_Arg *data    = (Partition_Arg *) arg;
  DATA_BLOCK    *block   = data->block;
  int            beg     = data->beg;
  int            end     = data->end;
  int64         *count   = data->count;

  char          *bases   = block->bases;
  int64         *boff    = block->boff+1;

  int   P2M2 = PAD2-2;

#ifdef DEBUG_PADDED
  int   P = (PAD_LEN-1)/2 + 1;
  int   M = (MIN_LEN-1)/2 + 1;
#endif

  int        force;
  int        i, p, q, x;
  uint64     c, u;
  char      *s, *t;

  uint64     min[MOD_LEN];
  uint64     mp, mc;
  int        m, n;
  int        last;
  int        nmin;

  int    b, v;
  int    o;

  nmin = 0;

  t = bases + boff[beg-1];
  for (i = beg; i < end; i++)
    { s = t;
      t = bases + boff[i];

      s += BC_PREFIX;
      q = (t-s) - 1;

      if (q < KMER)
        continue;

      m  = 0;
      mc = PAD_TOT;
      c  = u = 0;
      for (p = 0; p < KMER; p++)
        { x = s[p];
          c = ((c << 2) | Tran[x]) & PAD_MSK;
          u = (u >> 2) | Cran[x];
#ifdef DEBUG_PADDED
          if (PAD == PAD_LEVEL)
            printf(" %5d: %c %0*llx %0*llx",p,x,P,c,P,u);
          fflush(stdout);
#endif
          if (p >= PAD_L1)
            { if (u < c)
                mp = u;
              else
                mp = c;
              min[p] = mp;
              if (mp < mc)
                { m  = p;
                  mc = mp;
                }
#ifdef DEBUG_PADDED
              if (PAD == PAD_LEVEL)
                printf(" %0*llx <%5d:%0*llx>",P,mp,m,P,mc);
              fflush(stdout);
#endif
            }
#ifdef DEBUG_PADDED
          if (PAD == PAD_LEVEL)
            printf("\n");
          fflush(stdout);
#endif
        }
#ifdef DEBUG_PADDED
      if (PAD == PAD_LEVEL)
        printf(" -----\n");
      fflush(stdout);
#endif

      last = KMER-1;
      for (p = KMER; p < q; p++)
        { x = s[p];
          c = ((c << 2) | Tran[x]) & PAD_MSK;
          u = (u >> 2) | Cran[x];
          if (u < c)
            mp = u;
          else
            mp = c;
          min[p & MOD_MSK] = mp;

          force = (p-m >= MAX_SUPER);
          if (force || mp < mc)
            { o = (mc >> PAD2);
              v = count[o];
              b = P2M2;
              while (v < 0)
                { o = ((mc >> b) & 0x3) - v; 
                  v = count[o];
                  b -= 2;
                }
              count[o] += (p-last);
              nmin += 1;

#ifdef DEBUG_PADDED
              if (PAD == PAD_LEVEL)
                { if (force)
                    printf(" >");
                  else
                    printf(" +");
                  printf(" (%d) -> %0*llx / %0*llx -> %d\n",
                         p-last,M,mc>>PAD2,P,mc>>(b+2),o);
                }
              fflush(stdout);
#endif

              if (force)
                { mc = min[(++m) & MOD_MSK];
                  for (n = m+1; n <= p; n++)
                    { mp = min[n & MOD_MSK];
                      if (mp < mc)
                        { m  = n;
                          mc = mp;
                        }
                    }
                }
              else
                { m  = p;
                  mc = mp;
                }
              last = p;
            }

#ifdef DEBUG_PADDED
          if (PAD == PAD_LEVEL)
            printf(" %5d: %c %0*llx %0*llx %0*llx <%5d:%0*llx>\n",p,x,P,c,P,u,P,mp,m,P,mc);
          fflush(stdout);
#endif
        }

#ifdef DEBUG_PADDED
      if (PAD == PAD_LEVEL)
        printf("\n");
      fflush(stdout);
#endif
    }

  data->nmin = nmin;

  return (NULL);
}

  //  assign each core minimzer to a distribution bucket so that the k-mers
  //    will be distributed as evenly as possible.

static int64 *_CNT;

static int PSORT(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);

  if (_CNT[x] > _CNT[y])
    return (-1);
  else if (_CNT[x] == _CNT[y])
    return (0);
  else
    return (1);
}

static void assign_pieces(int nstates, int64 *count, int nparts, int64 pmer)
{ int64 *buck;
  int   *perm;
  int    i, j, n, x;
  int64  p, q;
  int64  v, t;

  perm = Malloc(sizeof(int)*nstates,"Allocating permutation array");
  buck = Malloc(sizeof(int64)*nparts,"Allocating permutation array");

  for (i = 0; i < nstates; i++)
    perm[i] = i;

  _CNT = count;
  qsort(perm,nstates,sizeof(int),PSORT);

#ifdef DEBUG_SCHEME
  printf("\nPieces:\n");
  for (i = 0; i < nstates; i++)
    if (count[perm[i]] > 0)
      printf("  %5d: %7lld -> %5.2f%%\n",perm[i],count[perm[i]],(100.*count[perm[i]])/pmer);
#endif

  for (i = 0; i < nparts; i++)
    buck[i] = 0;
  for (i = 0; i < nstates; i++)
    { x = perm[i];
      p = count[x];
      if (count[x] < 0)
        continue;
      if (p == 0)
        { count[x] = nparts-1;
          continue;
        }
      v = 0;
      for (j = 0; j < nparts; j++)  // if fits in >1 bucket, pick at random weighted
        if (buck[j] + p <= pmer)    //   by how much room each available bucket has
          v += pmer-buck[j];
      if (v == 0)
        { n = 0;
          for (j = 1; j < nparts; j++)
            if (buck[j] < buck[n])
              n = j;
          buck[n] += p;
          count[x] = n;
	}
      else
        { t = v*drand48();
          v = 0;
          for (j = 0; j < nparts; j++)
            if (buck[j] + p <= pmer)
              { v += pmer-buck[j];
                if (v >= t)
                  { buck[j] += p;
                    count[x] = j;
                    break;
                  }
              }
         }
    }

#ifdef DEBUG_SCHEME
  { int64 bmax, bmin;

    printf("\nPacking:\n");
    bmax = bmin = buck[0];
    for (i = 0; i < nparts; i++)
      { printf("  %7llu -> %5.2f%%\n",buck[i],(100.*buck[i])/pmer);
        if (bmax < buck[i])
          bmax = buck[i];
        else if (bmin > buck[i])
          bmin = buck[i];
      }
    printf("  Range %lld - %lld => %g%% diff\n",bmin,bmax,(100.*(bmax-bmin))/pmer);
  }
#endif

  for (i = 0; i < nstates; i++)
    if (perm[i] >= 0)
      { j = i;
        p = count[i];
        do
          { j = perm[j];
            q = count[j];
            count[j] = p;
            p = q;
          }
        while (j != i);
      }

  free(buck);
  free(perm);
}

#ifdef DEBUG_SCHEME

static void _print_tree(int lev, int i, int64 ktot, int64 *count)
{ int j, a;

  if (count[i] >= 0)
    printf(" %10lld %.3f%% (%d)\n",count[i],(100.*count[i])/ktot,lev);
  else
    { printf("\n");
      j = -count[i];
      for (a = 0; a < 4; a++)
        { printf("%*s -> %c:",2*lev,"",(char) Tran[a]);
          _print_tree(lev+1,j+a,ktot,count);
        }
    }
}

static void print_tree(int64 ktot, int64 *count)
{ int i;

  for (i = 0; i < MIN_TOT; i++)
    { printf(" %5d:",i);
      _print_tree(0,i,ktot,count);
    }
}

static void _print_ass(int lev, int i, int *part)
{ int j, a;

  if (part[i] >= 0)
    printf(" %d\n",part[i]);
  else
    { printf("\n");
      j = -part[i];
      for (a = 0; a < 4; a++)
        { printf("%*s -> %c:",2*lev,"",(char) Tran[a]);
          _print_ass(lev+1,j+a,part);
        }
    }
}

static void print_ass(int *part)
{ int i;

  for (i = 0; i < MIN_TOT; i++)
    { printf(" %5d:",i);
      _print_ass(0,i,part);
    }
}

#endif

static void _refine_tree(int lev, int i, int64 kthresh, int64 *count)
{ int j, a;

  lev += 1;
  if (count[i] >= 0)
    { if (count[i] > kthresh)
        {
#ifdef DEBUG_SCHEME
          printf("   U %d->%d: %lld > %lld (%d)\n",i,Min_States,count[i],kthresh,lev);
#endif
          count[i] = -Min_States;
          for (a = 0; a < 4; a++)
            count[Min_States++] = 0;
          if (lev > PAD)
#ifdef FIND_MTHRESH
            PAD += 1;
#else
            PAD += 2;
#endif
        }
      else
        count[i] = 0;
    }
  else
    { j = -count[i];
      for (a = 0; a < 4; a++)
        _refine_tree(lev,j+a,kthresh,count);
    }
}

static void refine_tree(int64 kthresh, int64 *count)
{ int i;

  for (i = 0; i < MIN_TOT; i++)
    _refine_tree(0,i,kthresh,count);
}

#ifdef HISTOGRAM_TEST

int64 *COUNT;

int HSORT(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);
  return (COUNT[y] - COUNT[x]);
}

#endif

static char DNA[4] = { 'a', 'c', 'g', 't' };

  //  Determine the base mapping, the core prefix trie, and the assignment of core prefixes
  //     to buckets (set up in globals PAD and Min_Map, et.c above)

int Determine_Scheme(DATA_BLOCK *block)
{ int64 *count;
  int    nreads, npieces;
  int64  ktot, mtot, kthresh;
  int64  max_count, last_max;
  int    i, j;

#if !(defined(DEBUG_MINIMIZER0) && defined(DEBUG_MINIMIZER1) && defined(DEBUG_PADDED))
  THREAD      threads[NTHREADS];
#endif
  Partition_Arg parmt[NTHREADS];

  (void) DNA;

  if (COMPRESS)
    compress_block(block);

  nreads   = block->nreads;
  npieces  = 2*NPARTS;

  parmt[0].beg = 0;
  for (i = 1; i < NTHREADS; i++)
    parmt[i].beg = parmt[i-1].end = (((int64) nreads) * i) / NTHREADS;
  parmt[NTHREADS-1].end = nreads;
  for (i = 0; i < NTHREADS; i++)
    parmt[i].block = block;

#ifdef DEBUG_SCHEME
  printf("  K = %d, Need %d pieces, Training on %lld bases\n",
         KMER,npieces,block->boff[nreads]-nreads);
#endif

  //  Determine best mapping of bases based on their frequency

  for (i = 1; i < NTHREADS; i++)
    pthread_create(threads+i,NULL,frequency_thread,parmt+i);
  frequency_thread(parmt);

  for (i = 1; i < NTHREADS; i++)
    pthread_join(threads[i],NULL);

  { int64 *freq, f;
    int    c;
#if defined(DEBUG_SCHEME) || defined(HISTOGRAM_TEST)
    int64  m, ftot;
#endif

    freq = parmt[0].freq;
    for (i = 0; i < 256; i++)
      for (j = 0; j < NTHREADS; j++)
	freq[i] += parmt[j].freq[i];

    freq[0] = freq['a'] + freq['A'];
    freq[1] = freq['c'] + freq['C'];
    freq[2] = freq['g'] + freq['G'];
    freq[3] = freq['t'] + freq['T'];

#if defined(DEBUG_SCHEME) || defined(HISTOGRAM_TEST)
    ftot = m = 0;
#endif
    for (i = 0; i < 4; i++)
      { c = 0;
        f = freq[i];
        for (j = 0; j < 4; j++)
          if (freq[j] < f)
            c += 1;
          else if (freq[j] == f && j < i)
            c += 1;
        Tran[i] = c;
#if defined(DEBUG_SCHEME) || defined(HISTOGRAM_TEST)
        if (f > m)
          m = f;
        ftot += f;
#endif
      }

#ifndef USE_MAPPING
    for (i = 0; i < 4; i++)
      Tran[i] = i;
#endif

    Tran['a'] = Tran['A'] = Tran[0];
    Tran['c'] = Tran['C'] = Tran[1];
    Tran['g'] = Tran['G'] = Tran[2];
    Tran['t'] = Tran['T'] = Tran[3];

#if defined(DEBUG_SCHEME) || defined(HISTOGRAM_TEST)
    printf("   Most freq = %lld  gc = %d%%\n",
           m,(int) ((100.*(freq[1]+freq[2]))/ftot));
    for (i = 0; i < 4; i++)
      printf(" Tran[%c] -> %lld (%lld / %7.4f)\n",DNA[i],Tran[i],freq[i],(100.*freq[i])/ftot);
    Tran[Tran['a']] = 'a';
    Tran[Tran['c']] = 'c';
    Tran[Tran['g']] = 'g';
    Tran[Tran['t']] = 't';
#endif
  }

  Dran['a'] = Dran['A'] = 0;
  Dran['c'] = Dran['C'] = 1;
  Dran['g'] = Dran['G'] = 2;
  Dran['t'] = Dran['T'] = 3;

  Fran['a'] = Fran['A'] = 3;
  Fran['c'] = Fran['C'] = 2;
  Fran['g'] = Fran['G'] = 1;
  Fran['t'] = Fran['T'] = 0;

  //  Iteratively determine the padding needed for each MIN_LEN-mer given the target # of pieces

  //  Initial prefix trie is the complete quartenary tree of height MIN_LEN

  count = Malloc(sizeof(int64)*MIN_TOT*NTHREADS,"Allocating count array");
  if (count == NULL)
    exit (1);

  PAD = 0;
  for (i = 0; i < MIN_TOT; i++)
    count[i] = 0;
  Min_States = MIN_TOT;

  last_max = 0;
  while (1)
    { int    o;

      // Compute stats with current prefix trie

      PAD2     = 2*PAD;
      PAD_LEN  = MIN_LEN + PAD;
      PAD_TOT  = (((int64) MIN_TOT) << PAD2);
      PAD_L1   = PAD_LEN - 1;
      PAD_MSK  = PAD_TOT - 1; 

      MAX_SUPER = KMER - PAD_L1;

      Cran['a'] = Cran['A'] = (Tran['t'] << (2*PAD_L1));
      Cran['c'] = Cran['C'] = (Tran['g'] << (2*PAD_L1));
      Cran['g'] = Cran['G'] = (Tran['c'] << (2*PAD_L1));
      Cran['t'] = Cran['T'] = (Tran['a'] << (2*PAD_L1));

      parmt[0].count = count;
      for (i = 1; i < NTHREADS; i++)
        { parmt[i].count = count + i*Min_States;
          memcpy(parmt[i].count,count,Min_States*sizeof(int64));
        }

#ifdef DEBUG_SCHEME
      printf("\n  Scheme: padding = %d(%d), # states = %d, target = %.3f%%\n",
             PAD,MAX_SUPER,Min_States,100./npieces);
#endif

#ifdef DEBUG_PADDED
      for (i = 0; i < NTHREADS; i++)
        padded_minimizer_thread(parmt+i);
#else
      for (i = 1; i < NTHREADS; i++)
        pthread_create(threads+i,NULL,padded_minimizer_thread,parmt+i);
      padded_minimizer_thread(parmt);

      for (i = 1; i < NTHREADS; i++)
        pthread_join(threads[i],NULL);
#endif

      //  For each leaf/core prefix, increase padding if frequency is still too large
      //  First just count how much bigger trie will be

      ktot = 0;
      for (i = 0; i < Min_States; i++)
        if (count[i] >= 0)
          { for (j = 1; j < NTHREADS; j++)
              count[i] += parmt[j].count[i];
            ktot += count[i];
          }
      mtot = 0;
      for (j = 0; j < NTHREADS; j++)
        mtot += parmt[j].nmin;

      if (ktot < (block->totlen - (KMER-1)*block->nreads)/2)
        { fprintf(stderr,"\n%s: Too much of the data is in reads less than the k-mer size\n",
                         Prog_Name);
          exit (1);
        }

      kthresh = ktot/npieces;
     
      max_count = 0;
      o = Min_States;
      for (i = 0; i < Min_States; i++)
        if (count[i] >= 0)
          { if (count[i] > kthresh)
              o += 4;
            if (count[i] > max_count)
              max_count = count[i];
          }

#ifdef MAXIMUM_TEST
      { int64 biggest;

        biggest = 0;
        for (i = 0; i < MIN_TOT; i++)
          if (count[i] > biggest)
            biggest = count[i];
        printf("max = %6.3f => %.2f pieces\n",(100.*biggest)/ktot,ktot/(1.*biggest));
        exit(0);
      }
#endif

#ifdef HISTOGRAM_TEST
      { int perm[MIN_TOT];
        int numc;

        for (i = 0; i < MIN_TOT; i++)
          perm[i] = i;
        COUNT = count;
        qsort(perm,MIN_TOT,sizeof(int),HSORT);
        if (MIN_LEN & 0x1)
          numc = MIN_TOT/2;
        else
          numc = (MIN_TOT + (MIN_TOT>>MIN_LEN))/2;
        for (i = 0; i < numc; i++)
          printf("%5d %5d %6.3f\n",i,perm[i],(100.*count[perm[i]])/ktot);
        printf("max = %6.3f => %.2f pieces\n",(100.*count[perm[0]])/ktot,ktot/(1.*count[perm[0]]));
        exit(0);
      }
#endif

#ifdef DEBUG_SCHEME
      printf("   S_k = %lld, S_m = %lld, super ave = %.1f, kthr = %lld\n",
             ktot,mtot,(1.*ktot)/mtot,kthresh);
      printf("   States needing spliting = %d",(o-Min_States)/4);
      if (PAD > 0)
        printf(", gain = %.3f",(1.*last_max)/max_count);
      printf("\n");
      // print_tree(ktot,count);
#endif
      
      //  If trie is core or improvment is less than 2%, stop

      if (o == Min_States)
        break;
      if (PAD > 0 && last_max < 1.02*max_count)
        { npieces = ktot/max_count+1;
          NPARTS  = npieces/2; 
          fprintf(stderr,"  Even split not possible, dividing into %d parts (%.1f x request)\n",
                         NPARTS,(1.*max_count)/kthresh); 
          break;
        }

      //  Else adjust trie size and actually do the refinement

      count = Realloc(count,sizeof(int64)*o*NTHREADS,"Expanding count array");
      if (count == NULL)
        exit (1);

      refine_tree(kthresh,count);
      Min_States = o;

#ifdef DEBUG_SCHEME
      printf("Max Bucket = %.3f%%  Padding = %d",(100.*max_count)/ktot,PAD);
#endif

      last_max = max_count;
    }

  if (KMER < PAD_LEN)
    { fprintf(stderr,"\n%s: K-mer must be at least %d\n",Prog_Name,PAD_LEN);
      exit (1);
    }

  if (VERBOSE)
    fprintf(stderr,"  Using %dbp of padding\n",PAD);

#ifdef FIND_MTHRESH
  exit (1);
#endif

#ifdef DEBUG_SCHEME
  print_tree(ktot,count);
#endif

  //  Assign buckets to core prefixes and finalize data structure

  assign_pieces(Min_States,count,NPARTS,ktot/NPARTS);

  Min_Part = (int *) count;
  for (i = 0; i < Min_States; i++)
    Min_Part[i] = count[i];

  Min_Part = Realloc(Min_Part,sizeof(int)*Min_States,"Finalizing count array");
  if (Min_Part == NULL)
    exit (1);

#ifdef DEBUG_SCHEME
  printf("\nPadded Assignments: # states = %d\n",Min_States);
  print_ass(Min_Part);
  fflush(stdout);
#endif

  return (MAX_SUPER);
}


/*******************************************************************************************
 *
 *  SUPER K-MER DISTRIBUTOR
 *
 *    void Split_Kmers(Input_Partition *io, char *root)
 *
 *       The input is scanned in parallel by threads working on each partition of it
 *       specified by io, and for each the data is partitioned into super-mers that
 *       are distibuted to blocks according to the core trie encoded in Min_Part.
 *       Split_Kmers calls Scan_All_Input in the io.c module, which spawns the threads,
 *       each feeding:
 *
 *             void Distribute_Block(DATA_BLOCK *block, int tid)
 *
 *       with the thread's input block at a time.  Each thread uses system IO doing its
 *       own buffering and statistics gathering.  The super-mer packets sent to each
 *       stream are bit compacted.
 *
 ********************************************************************************************/

//  Stuff ints and DNA super-mers into bit packed buffer

static inline IO_UTYPE *Stuff_Int(int64 n, int nbits, IO_UTYPE *buf, int *bitp)
{ int      rem;
  IO_UTYPE val;

  rem = *bitp;
  val = (IO_UTYPE) n;
#ifdef DEBUG_COMPRESS
  printf("v = %d/%llx b=%d/%016llx ->",nbits,val,rem,*buf);
#endif
  if (rem > nbits)
    { rem  -= nbits;
      *buf |= (val << rem);
#ifdef DEBUG_COMPRESS
      printf(" b = %d/%016llx\n",rem,*buf);
#endif
    }
  else if (rem == nbits)
    { *buf++ |= val;
      rem = IO_UBITS;
      *buf = 0;
#ifdef DEBUG_COMPRESS
      printf(" b = %d/%016llx\n",rem,buf[-1]);
#endif
    }
  else
    { *buf++ |= (val >> (nbits-rem));
      rem  += IO_UBITS - nbits;
      *buf = (val << rem);
#ifdef DEBUG_COMPRESS
      printf(" b = %d/%016llx|%016llx\n",rem,buf[-1],*buf);
#endif
    }
  *bitp = rem;
  return (buf);
}

static inline IO_UTYPE *Stuff_Seq(char *s, int len, IO_UTYPE *buf, int *bitp, int flip, int *pref)
{ int      rem, i;
  IO_UTYPE val, pre4;

  rem = *bitp;
#ifdef DEBUG_COMPRESS
  printf("seq %d/%d b=%d/%016llx\n",len,flip,rem,*buf);
#endif
  pre4 = 0;
  if (flip)
    { for (i = 1; i <= 4; i++)
        { val = Fran[(int) s[len-i]];
          // val = (IO_UTYPE) (s[len-i] ^ 0x3);
          if (rem > 2)
            { rem -= 2;
              *buf |= (val << rem);
#ifdef DEBUG_COMPRESS
              printf("  %llx: b = %2d/%016llx\n",val,rem,*buf);
#endif
            }
          else if (rem == 2)
            { *buf++ |= val;
              rem = IO_UBITS;
              *buf = 0;
#ifdef DEBUG_COMPRESS
              printf("  %llx: b = %2d/%016llx\n",val,rem,buf[-1]);
#endif
            }
          else
            { *buf++ |= (val >> 1); 
              rem  = IO_UBITS - 1;
              *buf = (val << rem);
#ifdef DEBUG_COMPRESS
              printf("  %llx: b = %2d/%016llx|%016llx\n",val,rem,buf[-1],*buf);
#endif
            }
          pre4 = (pre4 << 2 | val);
        }
      for (i = len-5; i >= 0; i--)
        { val = Fran[(int) s[i]];
          // val = (IO_UTYPE) (s[i] ^ 0x3);
          if (rem > 2)
            { rem -= 2;
              *buf |= (val << rem);
#ifdef DEBUG_COMPRESS
              printf("  %llx: b = %2d/%016llx\n",val,rem,*buf);
#endif
            }
          else if (rem == 2)
            { *buf++ |= val;
              rem = IO_UBITS;
              *buf = 0;
#ifdef DEBUG_COMPRESS
              printf("  %llx: b = %2d/%016llx\n",val,rem,buf[-1]);
#endif
            }
          else
            { *buf++ |= (val >> 1); 
              rem  = IO_UBITS - 1;
              *buf = (val << rem);
#ifdef DEBUG_COMPRESS
              printf("  %llx: b = %2d/%016llx|%016llx\n",val,rem,buf[-1],*buf);
#endif
            }
        }
    }
  else
    { for (i = 0; i < 4; i++)
        { val = Dran[(int) s[i]];
          // val = (IO_UTYPE) s[i];
          if (rem > 2)
            { rem -= 2;
              *buf |= (val << rem);
#ifdef DEBUG_COMPRESS
              printf("  %llx: b = %2d/%016llx\n",val,rem,*buf);
#endif
            }
          else if (rem == 2)
            { *buf++ |= val;
              rem = IO_UBITS;
              *buf = 0;
#ifdef DEBUG_COMPRESS
              printf("  %llx: b = %2d/%016llx\n",val,rem,buf[-1]);
#endif
            }
          else
            { *buf++ |= (val >> 1); 
              rem  = IO_UBITS - 1;
              *buf = (val << rem);
#ifdef DEBUG_COMPRESS
              printf("  %llx: b = %2d/%016llx|%016llx\n",val,rem,buf[-1],*buf);
#endif
            }
          pre4 = (pre4 << 2 | val);
        }
      for (i = 4; i < len; i++)
        { val = Dran[(int) s[i]];
          // val = (IO_UTYPE) s[i];
          if (rem > 2)
            { rem -= 2;
              *buf |= (val << rem);
#ifdef DEBUG_COMPRESS
              printf("  %llx: b = %2d/%016llx\n",val,rem,*buf);
#endif
            }
          else if (rem == 2)
            { *buf++ |= val;
              rem = IO_UBITS;
              *buf = 0;
#ifdef DEBUG_COMPRESS
              printf("  %llx: b = %2d/%016llx\n",val,rem,buf[-1]);
#endif
            }
          else
            { *buf++ |= (val >> 1); 
              rem  = IO_UBITS - 1;
              *buf = (val << rem);
#ifdef DEBUG_COMPRESS
              printf("  %llx: b = %2d/%016llx|%016llx\n",val,rem,buf[-1],*buf);
#endif
            }
        }
    }
#ifdef DEBUG_COMPRESS
  printf("   pref = %llx\n",pre4);
#endif
  *bitp = rem;
  *pref = (int) pre4;
  return (buf);
}

  //  Super-mer output files & buffering

typedef struct
  { int       stream;      //  Open stream
    int64     kmers;       //  Number of k-mers
    int64     nmers;       //  Numer of super-mers
    int64     fours[256];  //  fours[i] = Count of canonical super-mers with first byte i
    int       bbits;       //  Current bit offset in current word
    IO_UTYPE *bptrs;       //  Current word being stuffered in buffer
    IO_UTYPE *data;        //  End of buffer (start is at data-IO_BUF_LEN)
    int       nbits;       //  # of bits for last profile write (if DO_PROFILE)
  } Min_File;

  //  State for each thread between calls with a block of data

static Min_File **ogroup;   //  Vector of bucket files
static int64     *nfirst;   //  # of super-mers generated so far
static int       *nmbits;   //  # of bits currently being used for super-mer indices (if -p)
static int       *totrds;   //  # of reads processed
static int64     *totbps;   //  # of bps processed

void Distribute_Block(DATA_BLOCK *block, int tid)
{ int    nreads  = block->nreads;
  char  *bases   = block->bases;
  int64 *boff    = block->boff+1;

  Min_File *out   = ogroup[tid];
  int64     nidx  = nfirst[tid];
  int       nbits = nmbits[tid];
  int64     nlim  = (0x1ll << (nbits-1));
  int       KM1   = KMER-1;
  int       P2M2  = PAD2-2;
#ifdef DEBUG_DISTRIBUTE
  int       P     = (PAD_LEN-1)/2 + 1;
#endif

  int        force;
  int        i, p, q, x, w;
  uint64     c, u;
  char      *s, *t, *r;

  uint64     min[MOD_LEN];
  uint32     flp[MOD_LEN];
  uint64     mp, mc;
  int        m, n, b, y, o;
  int        last;

  int        pref, *prep = &pref;
  Min_File  *trg;
  IO_UTYPE  *ptr;
  int       *bit;

  if (COMPRESS)
    compress_block(block);

  totrds[tid] += nreads;
  totbps[tid] += block->totlen;

#if defined(SHOW_PACKETS)
  if (DO_PROFILE)
    printf("Index at %lld (%d/%lld)\n",nidx,nbits,nlim);
#endif

  trg = out; 
  ptr = trg->bptrs;
  bit = &(trg->bbits);

  t = bases + boff[-1];
  for (i = 0; i < nreads; i++)
    { s = t;
      t = bases + boff[i];

      s += BC_PREFIX;
      q = (t-s)-1;
      r = s-KM1;

#if defined(DEBUG_DISTRIBUTE) || defined(SHOW_PACKETS)
      printf("READ %d %lld\n",i+1,nidx);
      fflush(stdout);
#endif
      m  = 0;
      mc = PAD_TOT;
      c = u = 0;
      for (p = 0; p < KMER; p++)
        { x = s[p];
          c = ((c << 2) | Tran[x]) & PAD_MSK;
          u = (u >> 2) | Cran[x];
#ifdef DEBUG_DISTRIBUTE
          printf(" %5d: %c %0*llx",p,x,P,c);
          fflush(stdout);
#endif    
          if (p >= PAD_L1)
            { if ((flp[p] = (u < c)))
                mp = u;
              else
                mp = c;
              min[p] = mp;
              if (mp < mc)
                { m  = p;
                  mc = mp;
                }
#ifdef DEBUG_DISTRIBUTE
              printf(" %0*llx <%5d:%0*llx>",P,mp,m,P,mc);
              fflush(stdout);
#endif      
            }
#ifdef DEBUG_DISTRIBUTE
          printf("\n"); fflush(stdout);
          fflush(stdout);
#endif  
        }
#ifdef DEBUG_DISTRIBUTE
      printf(" -----\n");
      fflush(stdout);
#endif      

      last = KMER-1;
      for (p = KMER; p < q; p++)
        { x = s[p];
          c = ((c << 2) | Tran[x]) & PAD_MSK;
          u = (u >> 2) | Cran[x];
          w = p & MOD_MSK;
          if ((flp[w] = (u < c)))
            mp = u;
          else
            mp = c;
          min[w] = mp;

          force = (p-m >= MAX_SUPER);
        one_more:
          if (force || mp < mc)
            { o = (mc >> PAD2);
              b = Min_Part[o];
              y = P2M2;
              while (b < 0)
                { o = ((mc >> y) & 0x3) - b; 
                  b = Min_Part[o];
                  y -= 2;
                }

#ifdef DEBUG_DISTRIBUTE
              if (force)
                printf(" >  (%2d) -> %0*llx [%4d]\n",p-last,P,mc,b);
              else
                printf(" +  (%2d) -> %0*llx [%4d]\n",p-last,P,mc,b);
              fflush(stdout);
#endif

              n = p-last;
              trg = out + b;
              trg->kmers += n--;
              trg->nmers += 1;

              ptr = trg->bptrs;
              bit = &(trg->bbits);

              ptr = Stuff_Int(n,SLEN_BITS,ptr,bit);
              n += KMER;

              if (DO_PROFILE)
                { int   tbits;
                  int64 tlim;

                  ptr = Stuff_Seq(r+last,n,ptr,bit,0,prep);
                  trg->fours[pref] += 1;

                  if (nidx >= nlim)
                    { nlim <<= 1;
                      nbits += 1;
                    }

                  tbits = trg->nbits;
                  while (tbits < nbits)
                    { tlim = (0x1ll << (tbits-1));
                      ptr = Stuff_Int(tlim,tbits,ptr,bit);
                      if (ptr >= trg->data)
                        { IO_UTYPE *start = trg->data - IO_BUF_LEN;
                          write(trg->stream,start,IO_UBYTES*(ptr-start));
                          *start = *ptr;
                          ptr = start;
                        }
                      trg->nbits = ++tbits;
#ifdef SHOW_PACKETS
                      if (PACKET < 0 || b == PACKET)
                        printf("Index bumped to %d (%lld / %d)\n",tbits,tlim<<1,b);
#endif
                    }
                  ptr = Stuff_Int(nidx,nbits,ptr,bit);
                }

              else
                { ptr = Stuff_Seq(r+last,n,ptr,bit,flp[m & MOD_MSK],prep);
                  trg->fours[pref] += 1;
                }

#ifdef SHOW_PACKETS
              if (PACKET < 0 || b == PACKET)
                printf("   %6lld: %5d/%2d: %.*s[%c,%02x]\n",
                       nidx,last,n,n,r+last,flp[m & MOD_MSK]?'-':'+',pref);
              fflush(stdout);
#endif

              nidx += 1;
              if (ptr >= trg->data)
                { IO_UTYPE *start = trg->data - IO_BUF_LEN;
                  write(trg->stream,start,IO_UBYTES*(ptr-start));
                  *start = *ptr;
                  ptr = start;
                }
              trg->bptrs = ptr;

              if (force)
                { mc = min[(++m) & MOD_MSK];
                  for (n = m+1; n <= p; n++)
                    { mp = min[n & MOD_MSK];
                      if (mp < mc)
                        { m  = n;
                          mc = mp;
                        }
                    }
                }
              else
                { m  = p;
                  mc = mp;
                }
              last = p;
            }

#ifdef DEBUG_DISTRIBUTE
          if (p < q)
            printf(" %5d: %c %0*llx %0*llx <%5d:%0*llx>\n",p,x,P,c,P,mp,m,P,mc);
          else
            printf(" <<\n");
          fflush(stdout);
#endif
        }

      if (p == q)
        { mp = mc;
          force = 1;
          goto one_more;
        }

      if (DO_PROFILE)
        { trg->bptrs = Stuff_Int(MAX_SUPER,SLEN_BITS,ptr,bit);
#ifdef SHOW_PACKETS
          if (PACKET < 0 || b == PACKET)
            printf("   EOF\n");
#endif
        }
    }

  nfirst[tid] = nidx;
  nmbits[tid] = nbits;
}


 /********************************************************************************************
 *
 *  SUPER-MER DISTRIBUTION: TOP LEVEL
 *    void Split_Kmers(Input_Partition *io, char *root)
 *       Determine core prefix trie based on frequencies of a first big block, then for all blocks
 *       of the data base, send super-mer packets to partition files.
 *
 *********************************************************************************************/

void Split_Kmers(Input_Partition *io, char *root)
{ int           overflow;
  uint64        nfiles;
  int           nreads;
  int64         totlen;
int64 nids;

  Min_File     *out;
  IO_UTYPE     *buffers;

  //  Allocate output data structures

  nfiles = NPARTS*NTHREADS;

  overflow = IO_BUF_LEN
           + ((64 + 2*SLEN_BITS + 2*(MAX_SUPER+KMER-1)) - 1)/IO_UBITS
           + 1;

  out     = (Min_File *) Malloc(nfiles*sizeof(Min_File),"Allocating buffers");
  buffers = (IO_UTYPE *) Malloc(nfiles*overflow*IO_UBYTES,"Allocating buffers");
  if (out == NULL || buffers == NULL)
    exit (1);

  if (VERBOSE)
    { fprintf(stderr,"\nPhase 1: Partitioning K-mers into %lld Super-mer Files\n",nfiles);
      fflush(stderr);
    }

  //  Setup output data structures

  { int    t, p, n;
    char  *fname;
    int64 _zero = 0, *zero = &_zero;

    fname = (char *) Malloc(strlen(SORT_PATH)+strlen(root)+100,"Allocating file names");
    if (fname == NULL)
      exit (1);

    //  Remove any files that might still exist from a previous run of KMsplit in the
    //    same directory and same DB
  
    sprintf(fname,"rm -f %s/%s.*.T*",SORT_PATH,root);
    system(fname);

    p = 0;
    for (t = 0; t < NTHREADS; t++)
      for (n = 0; n < NPARTS; n++)
        { int f, i;

          sprintf(fname,"%s/%s.%d.T%d",SORT_PATH,root,n,t);
          f = open(fname,O_CREAT|O_TRUNC|O_WRONLY,S_IRWXU);
          if (f == -1)
            { fprintf(stderr,"\n%s: Cannot open external files in %s\n",
                             Prog_Name,SORT_PATH);
              exit (1);
            }

          out[p].stream = f;
          out[p].kmers  = 0;
          out[p].nmers  = 0;
          for (i = 0; i < 256; i++)
            out[p].fours[i] = 0;
          out[p].bptrs  = buffers + p*overflow;
          out[p].data   = buffers + p*overflow + IO_BUF_LEN;
          out[p].bbits  = IO_UBITS;
          out[p].nbits  = 17;
          *(out[p].bptrs) = 0;

#ifdef DEVELOPER
          write(f,zero,sizeof(int64));
          write(f,zero,sizeof(int64));
          write(f,zero,sizeof(int));
          write(f,zero,sizeof(int));
#endif
          write(f,zero,sizeof(int64));
          write(f,zero,sizeof(int64));
          write(f,zero,sizeof(int64));
          write(f,out[p].fours,sizeof(int64)*256);
          p += 1;
        }

    free(fname);
  }

  //  Ready to produce all super-mers and send to distribution/thread specific files

  { int   i;
    int64 val;

    nfirst = Malloc(sizeof(int64)*2*NTHREADS,"Allocating distribution globals");
    nmbits = Malloc(sizeof(int)*2*NTHREADS,"Allocating distribution globals");
    ogroup = Malloc(sizeof(Min_File *)*NTHREADS,"Allocating distribution globals");
    totrds = nmbits + NTHREADS;
    totbps = nfirst + NTHREADS;
    if (nfirst == NULL || nmbits == NULL || ogroup == NULL)
      exit (1);

    for (i = 0; i < NTHREADS; i++)
      { nfirst[i] = 0;
        nmbits[i] = 17;
        ogroup[i] = out + i*NPARTS;
        totrds[i] = 0;
        totbps[i] = 0;
      }

    Scan_All_Input(io);

    nids     = 0;
    nreads   = 0;
    totlen   = 0;
    for (i = 0; i < NTHREADS; i++)
      { nids   += nfirst[i];
        nreads += totrds[i];
        totlen += totbps[i];
      }

    RUN_BITS = 1;
    for (val = nids-1; val > 0; val >>= 1)
      RUN_BITS += 1;
    RUN_BYTES = (RUN_BITS+7) >> 3;

    free(ogroup);
    free(nmbits);
  }

  //  Determine k-mer & super-mer totals for headers, flush buffers,
  //    and set file headers

  { int   t, n, p;
    int   val;
    int64 s, m;
    int64 kmin, ktot, ntot;
    int   kwide, nwide, awide;

    ktot = 0;
    ntot = 0;
    for (n = 0; n < NPARTS; n++)
      { p = n;
        s = m = 0;
        for (t = 0; t < NTHREADS; t++)
          { s += out[p].kmers;
            m += out[p].nmers;
            p += NPARTS;
          }
        ktot += s;
        ntot += m;
      }

    awide = kwide = nwide = 0;
    if (VERBOSE)
      { fprintf(stderr,"  There are ");
        Print_Number((int64) nreads,0,stderr);
        fprintf(stderr," reads totalling ");
        Print_Number(totlen,0,stderr);
        fprintf(stderr," bps\n");
        kwide = Number_Digits(ktot);
        nwide = Number_Digits(ntot);
        awide = Number_Digits(ktot/ntot);
        kwide += (kwide-1)/3;
        nwide += (nwide-1)/3;
        awide += 2;
        if (Number_Digits(KMER)+4 > kwide)
          kwide = Number_Digits(KMER)+4;
        if (10 > nwide)
          nwide = 10;
        if (11 > awide)
          awide = 11;
      }

    NMAX = 0;
    KMAX = 0;
    kmin = SORT_MEMORY;
    if (VERBOSE)
      fprintf(stderr,"\n     Part:%*s%d-mer%*ssuper-mers%*save. length\n",
                     (kwide-Number_Digits(KMER))-2,"",KMER,nwide-8,"",awide-9,"");
    for (n = 0; n < NPARTS; n++)
      { p = n;
        s = m = 0;
        for (t = 0; t < NTHREADS; t++)
          { s += out[p].kmers;
            m += out[p].nmers;
            p += NPARTS;
          }
        if (VERBOSE)
          { fprintf(stderr,"    %5d:  ",n);
            Print_Number(s,kwide,stderr);
            fprintf(stderr,"  ");
            Print_Number(m,nwide,stderr);
            fprintf(stderr,"  %*.1f\n",awide,(1.*s)/m);
          }
        if (KMAX < s)
          KMAX = s;
        if (kmin > s)
          kmin = s;
        if (NMAX < m)
          NMAX = m;
      }
    if (VERBOSE)
      { fprintf(stderr,"      Sum:  ");
        Print_Number(ktot,kwide,stderr);
        fprintf(stderr,"  ");
        Print_Number(ntot,nwide,stderr);
        fprintf(stderr,"  %*.1f\n",awide,(1.*ktot)/ntot);

        fprintf(stderr,"\n      Range ");
        Print_Number(kmin,0,stderr);
        fprintf(stderr," - ");
        Print_Number(KMAX,0,stderr);
        fprintf(stderr," (%.2f%%)\n",(200.*(KMAX-kmin))/(KMAX+kmin));
      }

    KMAX_BYTES = 0;
    for (val = KMAX-1; val > 0; val >>= 8)
      KMAX_BYTES += 1;
    if (KMAX_BYTES < 2)  //  CMER_BYTES = KMAX_BYTES+1 must be 3 or more
      KMAX_BYTES = 2;

    p = 0;
    for (t = 0; t < NTHREADS; t++)
      for (n = 0; n < NPARTS; n++)
        { int       f     = out[p].stream;
          IO_UTYPE *start = out[p].data - IO_BUF_LEN;

          out[p].bptrs = Stuff_Int(0,SLEN_BITS,out[p].bptrs,&(out[p].bbits));

          if (out[p].bptrs > start || out[p].bbits < IO_UBITS)
            write(f,start,IO_UBYTES*((out[p].bptrs-start)+1));

          lseek(f,0,SEEK_SET);
#ifdef DEVELOPER
          write(f,&KMAX,sizeof(int64));
          write(f,&NMAX,sizeof(int64));
          write(f,&KMAX_BYTES,sizeof(int));
          write(f,&RUN_BITS,sizeof(int));
#endif
          write(f,&(out[p].kmers),sizeof(int64));
          write(f,&(out[p].nmers),sizeof(int64));
          write(f,nfirst+t,sizeof(int64));
          write(f,out[p].fours,sizeof(int64)*256);
          close(f);
          p += 1;
        }
  }

  free(nfirst);
  free(buffers);
  free(out);
  free(Min_Part);
}
