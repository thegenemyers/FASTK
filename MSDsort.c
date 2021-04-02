#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>

#include "libfastk.h"
#include "FastK.h"

#undef  IS_SORTED
#undef  SHOW_STUFF
#undef  DEBUG_CANONICAL

#ifdef DEBUG_CANONICAL

static char DNA[4] = { 'a', 'c', 'g', 't' };

static char *fmer[256], _fmer[1280];

#endif

#define SMAX  20
#define NMAX   3

#define THR0 15
#define THR1 15
#define THR2  8
#define GAP1  9
#define GAP2  4

static int S_thr0, S_thr1, S_thr2;
static int S_gap1, S_gap2;

static int    RSIZE;
static int    KSIZE;
static int64 *PARTS;
static uint8 *ARRAY;
static void  (*COUNT)(uint8 *,int64,Range *);

static inline void mycpy(uint8 *a, uint8 *b, int n)
{ while (n--)
    *a++ = *b++;
}

static inline int mycmp(uint8 *a, uint8 *b, int n)
{ while (n--)
    { if (*a++ != *b++)
        return (a[-1] < b[-1] ? -1 : 1);
    }
  return (0);
}

#ifdef IS_SORTED

static inline void sorted(uint8 *array, int64 asize, int digit)
{ int64 p, i;
  int   first = 1;
  int64 beg;

  for (p = RSIZE; p < asize; p += RSIZE)
    if (mycmp(array + (p-RSIZE) + 1,array + p + 1,KSIZE-1) > 0)
      { if (first)
          { first = 0;
            beg = (array-ARRAY)/RSIZE;
            printf("A[%lld-%lld]: %d\n",beg,beg+asize/RSIZE,digit);
          }
        printf("  Not sorted %12lld: ",p/RSIZE);
        for (i = 0; i < KSIZE; i++)
          printf(" %02x",array[(p-RSIZE)+i]);
        printf(" vs ");
        for (i = 0; i < KSIZE; i++)
          printf(" %02x",array[p+i]);
        printf("\n");
      }
}

#endif

static inline void gap_sort(uint8 *array, int asize, int gap, int cmp, int rem)
{ int    i, j;
  uint8  temp[RSIZE];
  uint8 *garray;

  garray = array + gap;
  for (i = gap; i < asize; i += RSIZE)
    { j = i-gap;
      if (mycmp(array+j,array+i,cmp) <= 0)
        continue;
      mycpy(temp,array+i,rem);
      mycpy(array+i,array+j,rem);
      for(j -= gap; j >= 0; j -= gap)
        { if (mycmp(array+j,temp,cmp) <= 0)
            break;
          mycpy(garray+j,array+j,rem);
        }
      mycpy(garray+j,temp,rem);
    }
}

static inline void shell_sort(uint8 *array, int asize, int digit, Range *rng)
{ int    cmp, rem;
  int    i, j;
  uint8 *garray;
  
  cmp    = KSIZE-digit;
  rem    = RSIZE-digit;
  garray = array+digit;

  // if (asize > S_thr1)
    // gap_sort(garray,asize,S_gap1,cmp,rem);
  if (asize > S_thr2)
    gap_sort(garray,asize,S_gap2,cmp,rem);
  gap_sort(garray,asize,RSIZE,cmp,rem);

  j = 0; 
  for (i = RSIZE; i < asize; i += RSIZE)
    if (mycmp(garray+j,garray+i,cmp) != 0)
      { COUNT(array+j,i-j,rng);
        j = i;
      }
  COUNT(array+j,asize-j,rng);
}

static void radix_sort(uint8 *array, int64 asize, int digit, int64 *alive, Range *rng)
{ int64  n, len[256];
  int    y, ntop;
  int    nzero[256];

  { uint8 *end[256];
    uint8 *u, *arrow = array + digit;
    int64  o;
    int    rems;
    int    e, x, z;

    uint8 *off[256];
    uint8  temp[RSIZE];
    uint8 *stack[SMAX];

    while (1)
      { e = arrow[0];
        for (o = RSIZE; o < asize; o += RSIZE)
          { if (arrow[o] != e)
              break;
          }
        if (o < asize)
          break;
        digit += 1;
        if (digit >= KSIZE)
          { COUNT(array,asize,rng);
            return;
          }
        arrow += 1;
      }

    ntop = 1;
    nzero[0] = e;
    alive[e] = o;
    for (; o < asize; o += RSIZE)
      { x = arrow[o];
        if (alive[x] == 0)
          nzero[ntop++] = x;
        alive[x] += RSIZE;
      }

    u = arrow;
    if (ntop <= NMAX)
      { for (y = 1; y < ntop; y++)
          { x = nzero[y];
            for (z = y-1; z >= 0; z--)
              if (nzero[z] < x)
                break;
              else
                nzero[z+1] = nzero[z];
            nzero[z+1] = x;
          }
        for (y = 0; y < ntop; y++)
          { x = nzero[y];
            len[x] = alive[x];
            alive[x] = 0;
            off[x] = u;
            end[x] = u += len[x];
          }
      }
    else
      { ntop = 0;
        for (x = 0; x < 256; x++)
          if (alive[x] > 0)
            { len[x] = alive[x];
              alive[x] = 0;
              off[x] = u;
              end[x] = u += len[x];
              nzero[ntop++] = x;
            }
      }

    rems = RSIZE-digit;
    for (y = 0; y < ntop; y++)
      { uint8   *p;
        int      t, s;
        int      z;

        x = nzero[y];
        while (off[x] < end[x])
          { t = *off[x];

            if (t == x)
              off[x] += RSIZE;
	    else
              { s = 0;
                stack[s++] = off[x];
                while (s < SMAX)
                  if (t == x)
                    { off[x] += RSIZE;
                      break;
                    }
                  else
                    { u = off[t];
                      while ((z = *u) == t)
                        u += RSIZE;
                      off[t] = u+RSIZE;
                      stack[s++] = u;
                      t = z;
                    }

                u = stack[--s];
                mycpy(temp,u,rems);
	        while (s > 0)
                  { p = stack[--s];
                    mycpy(u,p,rems);
                    u = p;
                  }
                mycpy(u,temp,rems);
              }
          }
      }
  }
  
  digit += 1;
  if (digit < KSIZE)
    for (y = 0; y < ntop; y++)
      { n = len[nzero[y]];
        if (n > S_thr0)
          radix_sort(array, n, digit, alive, rng);
        else if (n > RSIZE)
          shell_sort(array, n, digit, rng);
        else if (n > 0)
          COUNT(array,n,rng);
        array += n;
      }
  else
    for (y = 0; y < ntop; y++)
      { n = len[nzero[y]];
        COUNT(array,n,rng);
        array += n;
      }
}

static int INIT_COUNTS;

static void *sort_thread(void *arg) 
{ Range *param = (Range *) arg;

  int      beg   = param->beg;
  int      end   = param->end;
  int64    off   = param->off;
  int64   *khist = param->khist;
  int64   *count = param->count;

  int      x;
  int64    alive[256];

  if (KSIZE <= 0)
    return (NULL);

  for (x = 0; x < 256; x++)
    alive[x] = khist[x] = 0;

  if (INIT_COUNTS)
    { for (x = 0; x < 0x8000; x++)
        count[x] = 0;
      param->max_inst = 0;
    }

  for (x = beg; x < end; x++)
    { if (PARTS[x] == 0)
        continue;

#ifdef SHOW_STUFF
      printf("Thread %3d: %12lld - %12lld\n",x,off,off+PARTS[x]);
#endif

      param->byte1 = x;
      radix_sort(ARRAY + off, PARTS[x], 1, alive, param);

      off += PARTS[x];
    }

  return (NULL);
}

static void msd_sort(uint8 *array, int64 nelem, int rsize, int ksize,
                     int64 *part, int nthreads, Range *parms)
{
#ifndef SHOW_STUFF
  pthread_t     threads[nthreads];
#endif

  int   x, n, beg;
  int64 sum, thr, off;
  int64 asize;

  asize = nelem*rsize;

  ARRAY = array;
  PARTS = part;
  RSIZE = rsize;
  KSIZE = ksize;

  S_thr0 = THR0*RSIZE;
  S_thr1 = THR1*RSIZE;
  S_thr2 = THR2*RSIZE;
  S_gap1 = GAP1*RSIZE;
  S_gap2 = GAP2*RSIZE;

  n   = 0;
  thr = asize / nthreads;
  off = 0;
  sum = 0;
  beg = 0;
  for (x = 0; x < 256; x++)
    { sum += part[x];
      if (sum >= thr)
        { parms[n].end = x+1;
          parms[n].beg = beg;
          parms[n].off = off;
          n  += 1;
          thr = (asize * (n+1))/nthreads;
          beg = x+1;
          off = sum;
        }
    }
  while (n < nthreads)
    { parms[n].end = 256;
      parms[n].beg = 256;
      parms[n].off = off;
      n += 1;
    }

#ifdef SHOW_STUFF
  for (x = 0; x < nthreads; x++)
    sort_thread(parms+x);
#else
  for (x = 1; x < nthreads; x++)
    pthread_create(threads+x,NULL,sort_thread,parms+x);

  sort_thread(parms);

  for (x = 1; x < nthreads; x++)
    pthread_join(threads[x],NULL);
#endif

#ifdef IS_SORTED
  sum = 0;
  for (x = 0; x < 256; x++)
    { sorted(array+sum,part[x],0);
      sum += part[x];
    }
#endif

  array[asize] = 1;
}

static int  ROFF;
static int  RSHIFT;

static inline void count_smers(uint8 *array, int64 asize, Range *rng)
{ int    sln, k;
  uint8 *asp;

  int    fb, fl, fs;
  int    rb, rl, rs;
  int    kb, hb;
  uint8 *f, *r;
  int64 *khist = rng->khist;

  (void) asize;

  asp = array+SMER_BYTES;
  sln = *asp++;
  for (k = 1; k < SLEN_BYTES; k++)
    sln = (sln << 8) + *asp++;

  *array = rng->byte1;

#ifdef DEBUG_CANONICAL
  printf(" %3d:",sln);
  f = array;
  for (k = 0; k < sln+KMER; k += 4)
    printf(" %s",fmer[*f++]);
  printf("\n");
  printf(" %3d:",sln);
  f = array;
  for (k = 0; k < sln+KMER; k += 4)
    printf(" %s",fmer[Comp[*f++]]);
  printf("\n");
#endif

  f  = array;
  fl = *f;
  fs = 0;

  r  = array + ROFF;
  rs = RSHIFT;
  rl = Comp[*r];
  rb = 0;
  if (RSHIFT != 8)
    rb = (rl << 8) | Comp[r[-1]];

  for (k = 0; k <= sln; k += 1)
    { if (fs == 0)
        { kb  = fl;
          fl  = *++f;
          fb  = (kb << 8) | fl;
          fs  = 6;
        }
      else
        { kb = (fb >> fs) & 0xff;
          fs -= 2;
        }
      if (rs == 8)
        { hb  = rl;
          rl  = Comp[*++r];
          rb  = (rl << 8) | hb;
          rs  = 2;
        }
      else
        { hb = (rb >> rs) & 0xff;
          rs += 2;
        }
#ifdef DEBUG_CANONICAL
      printf("   %d / %s%s[%s] / %s\n",fs, fmer[fb>>8], fmer[fb&0xff], fmer[fl], fmer[kb]);
      printf("     %d / %s%s[%s] / %s\n",rs, fmer[rb&0xff], fmer[rb>>8], fmer[rl], fmer[hb]);
#endif
      if (kb < hb)
        khist[kb] += 1;
      else
        khist[hb] += 1;
    }

  *array = 1;
}

void Supermer_Sort(uint8 *array, int64 nelem, int rsize, int ksize,
                   int64 *part, int nthreads, Range *panel)
{ ROFF   = ((KMER-1)>>2);
  RSHIFT = 2*(KMER & 0x3);
  if (RSHIFT == 0)
    RSHIFT = 8;
  COUNT = count_smers;
  INIT_COUNTS = 0;

#ifdef DEBUG_CANONICAL
  { char *t;
    int   i, l3, l2, l1, l0;

    i = 0;
    t = _fmer;
    for (l3 = 0; l3 < 4; l3++)
     for (l2 = 0; l2 < 4; l2++)
      for (l1 = 0; l1 < 4; l1++)
       for (l0 = 0; l0 < 4; l0++)
         { fmer[i] = t;
           *t++ = DNA[l3];
           *t++ = DNA[l2];
           *t++ = DNA[l1];
           *t++ = DNA[l0];
           *t++ = 0;
           i += 1;
         }
  }
#endif

  return (msd_sort(array,nelem,rsize,ksize,part,nthreads,panel));
}

static inline void hist_kmers(uint8 *array, int64 asize, Range *rng)
{ int k, cnt;

  cnt = 0;
  for (k = KMER_BYTES; k < asize; k += RSIZE)
    cnt += *((uint16 *) (array + k));

  if (cnt >= 0x7fff)
    { rng->count[0x7fff] += 1;
      rng->max_inst += cnt;
      cnt = 0x7fff;
    }
  else
    rng->count[cnt] += 1;

  *((uint16 *) (array + KMER_BYTES)) = cnt;

  *array = 1;
}

static inline void invert_kmers(uint8 *array, int64 asize, Range *rng)
{ int    k, cnt;
  int64 *khist = rng->khist;

  cnt = 0;
  for (k = KMER_BYTES; k < asize; k += RSIZE)
    cnt += *((uint16 *) (array + k));

  if (cnt >= 0x7fff)
    { rng->count[0x7fff] += 1;
      rng->max_inst += cnt;
      cnt = 0x7fff;
    }
  else
    rng->count[cnt] += 1;

  for (k = KMER_BYTES; k < asize; k += RSIZE)
    *((uint16 *) (array + k)) = cnt;

  for (k = RSIZE-1; k < asize; k += RSIZE)
    khist[array[k]] += 1;

  *array = 1;
}

void Weighted_Kmer_Sort(uint8 *array, int64 nelem, int rsize, int ksize,
                        int64 *part, int nthreads, Range *panel)
{ if (DO_PROFILE)
    COUNT = invert_kmers;
  else
    COUNT = hist_kmers;
  INIT_COUNTS = 1;
  return (msd_sort(array,nelem,rsize,ksize,part,nthreads,panel));
}
