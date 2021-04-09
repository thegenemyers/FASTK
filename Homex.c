/*********************************************************************************************\
 *
 *  Code to collect statistics on homopolymer error rates
 *
 *  Author:  Gene Myers
 *  Date  :  October, 2020
 *
 *********************************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>

#undef    DEBUG_PARTITION
#undef    DEBUG_QUEUE
#undef    DEBUG_DATA_POINT

#include "libfastk.h"

static char *Usage = "-e<int> -g<int>:<int> <source_root>[.ktab]";

#define MAX_HOMO_LEN 10

/****************************************************************************************
 *
 *  Print & compare utilities
 *
 *****************************************************************************************/

#define  COUNT_OF(p) (*((uint16 *) (p+kbyte)))

static char dna[4] = { 'a', 'c', 'g', 't' };

static char *fmer[256], _fmer[1280];

static void setup_fmer_table()
{ char *t;
  int   i, l3, l2, l1, l0;

  i = 0;
  t = _fmer;
  for (l3 = 0; l3 < 4; l3++)
   for (l2 = 0; l2 < 4; l2++)
    for (l1 = 0; l1 < 4; l1++)
     for (l0 = 0; l0 < 4; l0++)
       { fmer[i] = t;
         *t++ = dna[l3];
         *t++ = dna[l2];
         *t++ = dna[l1];
         *t++ = dna[l0];
         *t++ = 0;
         i += 1;
       }
}

#if defined(DEBUG_PARTITION)

static void print_seq(uint8 *seq, int len)
{ int i, b, k;

  b = (len >> 2);
  for (i = 0; i < b; i++)
    printf("%s",fmer[seq[i]]);
  k = 6;
  for (i = b << 2; i < len; i++)
    { printf("%c",dna[(seq[b] >> k) & 0x3]);
      k -= 2;
    }
}

#endif

static inline void mycpy(uint8 *a, uint8 *b, int n)
{ while (n--)
    *a++ = *b++;
}


/****************************************************************************************
 *
 *  Find Haplotype Pairs
 *
 *****************************************************************************************/

static int ERROR;
static int GOOD_LOW;
static int GOOD_HGH;

static inline int mypref(uint8 *a, uint8 *b, int n)
{ int   i;
  uint8 x, y;
  
  for (i = 0; i <= n; i += 4)
    { if (*a != *b)
        { x = *a;
          y = *b;
          if ((x & 0xf0) != (y & 0xf0))
            if ((x & 0xc0) != (y & 0xc0))
              return (i);
            else
              return (i + 1);
          else
            if ((x & 0xfc) != (y & 0xfc))
              return (i + 2);
            else
              return (i + 3);
        }
      a += 1;
      b += 1;
    }
  return (n+1);
}

static uint8  base[]   = { 0xc0, 0x30, 0x0c, 0x03 };
static uint8  shift[]  = { 6, 4, 2, 0 };

#define SYMBOL(ptr,pos) ((ptr[pos>>2] & base[pos&0x3]) >> shift[pos&0x3])

static inline int mybpcmp(uint8 *a, uint8 *b, int x, int y, int n)
{ int t, u;

  while (n-- > 0)
    { t = SYMBOL(a,x);
      u = SYMBOL(b,y);
      if (t < u)
        return (1);
      else if (t > u)
        return (-1);
      x += 1;
      y += 1;
    }
  return (0);
}

typedef struct
  { int64   correct;
    int64   lessone;
    int64   plusone;
  } Point;

typedef Point Profile[4][MAX_HOMO_LEN+1]; 

Profile *Count_Homopolymer_Errors(Kmer_Stream *T)
{ static Profile profile;

  int    kmer  = T->kmer;
  int    tbyte = T->tbyte;
  int    kbyte = T->kbyte;

  Point *counter;
  int    khalf, klong, kbase, kchkl, kextn;

  uint8  suffs[]  = { 0x00, 0xc0, 0xf0, 0xfc };
  uint8  abyte[]  = { 0x00, 0x55, 0xaa, 0xff };

  uint8 *cache, *cptr, *ctop;

  int    i;
  int64  fing[4];
  int64  fend[4];
  int64  fbeg[5];

  int    hlen, hsym;
  uint8  suffix[kbyte];

  int   a, b, advn[4];
  int   cn[4];
#if defined(DEBUG_PARTITION) || defined(DEBUG_QUEUE)
  int64 ridx;
#endif

  setup_fmer_table();

  for (hsym = 0; hsym < 4; hsym++)
    for (hlen = 0; hlen < MAX_HOMO_LEN; hlen++)
      { profile[hsym][hlen].correct = 0;
        profile[hsym][hlen].lessone = 0;
        profile[hsym][hlen].plusone = 0;
      }

  khalf = kmer/2;
  klong = khalf - (MAX_HOMO_LEN/2);
  if (klong < 10)
    { fprintf(stderr,"%s: A k-mer length of at least %d is needed\n",Prog_Name,20+MAX_HOMO_LEN);
      exit (1);
    }
  klong -= 1;

  cache = Malloc(4097*tbyte,"Allocating entry buffer");
  cptr  = cache;
  ctop  = cache + 4096*tbyte;
  fbeg[4] = 0;

  First_Kmer_Entry(T);
  while (T->csuf != NULL)
    { Current_Entry(T,cache);
      hlen = khalf-1;
      hsym = SYMBOL(cache,hlen);
      for (hlen--; hlen >= klong-1; hlen--)
        if (SYMBOL(cache,hlen) != hsym)
          break;
      hlen += 1;

      if (hlen <= klong)
        { cptr = cache + tbyte;
          for (Next_Kmer_Entry(T); T->csuf != NULL; Next_Kmer_Entry(T))
            { int x = mypref(Current_Entry(T,cptr),cache,khalf); 
              if (x < khalf)
                break;
            }
          mycpy(cache,cptr,tbyte);
          continue;
        }

      hlen = khalf-hlen;
#if defined(DEBUG_PARTITION) || defined(DEBUG_QUEUE)
      ridx = T->cidx;
#endif
#ifdef DEBUG_PARTITION
      printf(" %lld: ",ridx);
      print_seq(cache,kmer);
      printf("  Len = %d  Sym = %c\n",hlen,dna[hsym]);
      fflush(stdout);
#endif

      { int k;

        mycpy(suffix,cache,kbyte);
        k = (khalf>>2);
        suffix[k] = (suffix[k] & suffs[khalf&0x3]) | (abyte[hsym] & ~suffs[khalf&0x3]);
        for (k++; k < kbyte; k++)
          suffix[k] = abyte[hsym];
#ifdef DEBUG_PARTITION
        printf("      ");
        print_seq(suffix,kmer);
        printf("\n");
        fflush(stdout);
#endif
      }

      kbase = khalf + (hlen-1);
      kchkl = khalf + (hlen+2);
      kextn = kmer - kbase;
      for (i = 0; i <= 3; i++)
        fend[i] = -1;
    
      cptr = cache;
      for (; T->csuf != NULL; Next_Kmer_Entry(T))
        { int x = mypref(Current_Entry(T,cptr),suffix,kchkl); 
          if (x < khalf)
            break;
          if (cptr >= ctop)
            { int64 cidx = cptr-cache;
              int64 cmax = ((cidx*14)/(10*tbyte) + 2048)*tbyte;
              cache = Realloc(cache,cmax+tbyte,"Reallocting entry buffer");
              ctop  = cache + cmax;
              cptr  = cache + cidx;
            }
          x -= kbase;
          if (0 <= x && x <= 3)
            { if (fend[x] < 0)
                fbeg[x] = cptr - cache;
              fend[x] = cptr - cache;
            }
          cptr += tbyte;
        }

#ifdef DEBUG_PARTITION
      for (i = 0; i <= 3; i++)
        if (fend[i] >= 0)
          printf(" %lld-%lld",ridx+fbeg[i]/tbyte,ridx+fend[i]/tbyte);
      printf(" >> %lld\n",ridx+(cptr-cache)/tbyte);
      fflush(stdout);
#endif

      if (fend[1] < 0 && fend[2] < 0)
        continue;

      for (i = 3; i >= 0; i--)
        if (fend[i] < 0)
          fing[i] = fend[i] = fbeg[i] = 0;
        else
          { fing[i] = fbeg[i];
            fend[i] += tbyte;
          }

#ifdef DEBUG_PARTITION
      for (i = 0; i <= 3; i++)
        if (fend[i] == 0)
          printf(" ***");
        else
          printf(" %lld-%lld",ridx+fbeg[i]/tbyte,ridx+fend[i]/tbyte);
      printf(" >> %lld\n",ridx+(cptr-cache)/tbyte);
      fflush(stdout);
#endif

#define ADD(i) advn[a++] = i;

#define SET(i)	\
{ a = 0;	\
  b = i;	\
  ADD(i);	\
}

      counter = profile[hsym];
      hlen  <<= 1;

      while (1)
        { for (i = 0; i <= 3; i++)
            if (fing[i] < fend[i])
              break;
          if (i > 3)
            break;
          SET(i);
          for (i++; i <= 3; i++)
            if (fing[i] < fend[i])
              { int v = mybpcmp(cache+fing[b],cache+fing[i],kbase+b,kbase+i,kextn-i);
                if (v == 0)
                  ADD(i)
                else if (v < 0)
                  SET(i)
              }

#ifdef DEBUG_QUEUE
          for (i = 0; i < a; i++)
            printf(" %d(%d) %lld",advn[i],COUNT_OF(cache+fing[advn[i]]),ridx+fing[advn[i]]/tbyte);
#endif

          cn[0] = cn[1] = cn[2] = cn[3] = 0;
          for (i = 0; i < a; i++)
            { b = advn[i];
              cn[b] = COUNT_OF(cache+fing[b]);
              fing[b] += tbyte;
              if (fing[b] == fbeg[b+1])
                fing[b] = fend[b+1];
            }

          if (GOOD_LOW <= cn[1] && cn[1] <= GOOD_HGH && cn[0] <= ERROR && cn[2] <= ERROR)
            { counter[hlen].correct += cn[1];
              counter[hlen].lessone += cn[0]; 
              counter[hlen].plusone += cn[2]; 
#ifdef DEBUG_DATA_POINT
              printf(" -> %d%c %d:%d:%d\n",hlen,dna[hsym],cn[0],cn[1],cn[2]);
              fflush(stdout);
#endif
            }
          else if (GOOD_LOW <= cn[2] && cn[2] <= GOOD_HGH && cn[1] <= ERROR && cn[3] <= ERROR)
            { if (hlen < MAX_HOMO_LEN)
                { counter[hlen+1].correct += cn[2];
                  counter[hlen+1].lessone += cn[1]; 
                  counter[hlen+1].plusone += cn[3]; 
#ifdef DEBUG_DATA_POINT
                  printf(" -> %d%c %d:%d:%d\n",hlen+1,dna[hsym],cn[1],cn[2],cn[3]);
              fflush(stdout);
#endif
                }
            }
#ifdef DEBUG_QUEUE
          else
            printf("\n");
#endif
        }
    }

  return (&profile);
}


/****************************************************************************************
 *
 *  Main
 *
 *****************************************************************************************/

int main(int argc, char *argv[])
{ Kmer_Stream *T;
  Profile     *P;

  { int    i, j, k;
    int    flags[128];
    char  *eptr, *fptr;

    (void) flags;

    ARG_INIT("Homex");

    ERROR    = -1;
    GOOD_LOW = -1;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
            break;
          case 'e':
            ERROR = strtol(argv[i]+2,&eptr,10);
            if (eptr > argv[i]+2 && *eptr == '\0')
              { if (ERROR < 1 || ERROR > 0x7fff)
                  { fprintf(stderr,"%s: Error threshold %d is out of range\n",
                                   Prog_Name,ERROR);
                    exit (1);
                  }
                break;
              }
            fprintf(stderr,"%s: Syntax of -e option invalid -e<int>\n",Prog_Name);
            exit (1);
          case 'g':
            GOOD_LOW = strtol(argv[i]+2,&eptr,10);
            if (eptr > argv[i]+2)
              { if (GOOD_LOW < 1 || GOOD_LOW > 0x7fff)
                  { fprintf(stderr,"%s: Minimum valid count %d is out of range\n",
                                   Prog_Name,GOOD_LOW);
                    exit (1);
                  }
                if (*eptr == ':')
                  { GOOD_HGH = strtol(eptr+1,&fptr,10);
                    if (fptr > eptr+1 && *fptr == '\0')
                      { if (GOOD_HGH < 1 || GOOD_HGH > 0x7fff)
                          { fprintf(stderr,"%s: Maximum valid count %d is out of range\n",
                                           Prog_Name,GOOD_HGH);
                            exit (1);
                          }
                        if (GOOD_LOW > GOOD_HGH)
                          { fprintf(stderr,"%s: Good count range is invalid\n",Prog_Name);
                            exit (1);
                          }
                        break;
                      }
                  }
              }
            fprintf(stderr,"%s: Syntax of -g option invalid -g<int>:<int>\n",Prog_Name);
            exit (1);
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -e: Counts <= this value are considered errors.\n");
        fprintf(stderr,"      -g: Counts in this range are considered correct.\n");
        exit (1);
      }

    if (ERROR < 0)
      { fprintf(stderr,"%s: Must give error count threshold -e\n",Prog_Name);
        exit (1);
      }
    if (GOOD_LOW < 0)
      { fprintf(stderr,"%s: Must give good count range -g\n",Prog_Name);
        exit (1);
      }
  }

  T = Open_Kmer_Stream(argv[1]);

  P = Count_Homopolymer_Errors(T);

  Free_Kmer_Stream(T);

  { int h;

    printf("\n              -1      Good          +1      Error Rate\n\n");
    for (h = 2; h <= MAX_HOMO_LEN; h++)
      { int64 cc = (*P)[0][h].correct + (*P)[3][h].correct;
        int64 cl = (*P)[0][h].lessone + (*P)[3][h].lessone;
        int64 cp = (*P)[0][h].plusone + (*P)[3][h].plusone;

        printf(" %2d at: %10lld %10lld %10lld -> %.1f%%\n",h,cl,cc,cp,(100.*(cl+cp))/(cc+cl+cp));
      }

    printf("\n");
    for (h = 2; h <= MAX_HOMO_LEN; h++)
      { int64 cc = (*P)[1][h].correct + (*P)[2][h].correct;
        int64 cl = (*P)[1][h].lessone + (*P)[2][h].lessone;
        int64 cp = (*P)[1][h].plusone + (*P)[2][h].plusone;

        printf(" %2d cg: %10lld %10lld %10lld -> %.1f%%\n",h,cl,cc,cp,(100.*(cl+cp))/(cc+cl+cp));
      }
  }

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
