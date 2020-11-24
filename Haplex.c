/*********************************************************************************************\
 *
 *  Code for produce all potential haplotype k-mers with a single SNP in the center of
 *    the k-mers
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

#undef DEBUG_PARTITION

#include "libfastk.h"

static char *Usage = "[-h<int>:<int>] <source>[.ktab]";

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

static void print_seq(uint8 *seq, int len)
{ int i, b, k;

  b = len >> 2;
  for (i = 0; i < b; i++)
    printf("%s",fmer[seq[i]]);
  k = 6;
  for (i = b << 2; i < len; i++)
    { printf("%c",dna[seq[b] >> k]);
      k -= 2;
    }
}

static inline int mycmp(uint8 *a, uint8 *b, int n)
{ while (n-- > 0)
    { if (*a++ != *b++)
        return (a[-1] < b[-1] ? -1 : 1);
    }
  return (0);
}

static inline void mycpy(uint8 *a, uint8 *b, int n)
{ while (n--)
    *a++ = *b++;
}


/****************************************************************************************
 *
 *  Find Haplotype Pairs
 *
 *****************************************************************************************/

static int HAPLO_LOW, HAPLO_HGH;

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

void Find_Haplo_Pairs(Kmer_Stream *T)
{ int    kmer  = T->kmer;
  int    tbyte = T->tbyte;
  int    kbyte = T->kbyte;

  int    khalf;
  uint8  prefs[] = { 0x3f, 0x0f, 0x03, 0x00 };
  int    mask, offs, rem;

  uint8 *iptr;
  uint8 *cache, *cptr, *ctop;

  int    f, i;
  int    index[4];
  uint8 *finger[4];
  uint8 *flimit[4];

  int    a, advn[4];
  int    c, good[4];
  int    mc, hc;
  uint8 *mr, *hr;

  setup_fmer_table();

  khalf = kmer/2;
  mask  = prefs[khalf&0x3]; 
  offs  = (khalf >> 2) + 1;
  rem   = ((kmer+3) >> 2) - offs;

#ifdef DEBUG_PARTITION
  printf("Extension = K[%d]&%02x . K[%d..%d)\n",offs-1,mask,offs,offs+rem);
#endif

  cache = Malloc(4096*tbyte,"Allocating entry buffer");
  cptr  = cache;
  ctop  = cache + 4096*tbyte;

  iptr = First_Kmer_Entry(T);
  while (iptr != NULL)
    { f = 0;
      cptr = cache;
      index[f++] = 0;
      mycpy(cptr,iptr,tbyte);
      for (iptr = Next_Kmer_Entry(T); iptr != NULL; iptr = Next_Kmer_Entry(T))
        { int x = mypref(cptr,iptr,khalf); 
          cptr += tbyte;
          if (x < khalf)
            break;
          if (x == khalf)
            index[f++] = cptr-cache;
          if (cptr >= ctop)
            { int64 cidx = ctop-cache;
              int64 cmax = cidx*1.4 + 2048*tbyte; 
              cache = Realloc(cache,cmax,"Reallocting entry buffer");
              ctop  = cache + cmax;
              cptr  = cache + cidx;
            }
          mycpy(cptr,iptr,tbyte);
        }

#ifdef DEBUG_PARTITION
      printf("part %d",f);
      for (i = 0; i < f; i++)
        printf(" %d",index[i]/tbyte);
      printf(" %ld\n",(cptr-cache)/tbyte);
#endif

      if (f <= 1)
        continue;

      for (i = 0; i < f; i++)
        finger[i] = cache + index[i];
      for (i = 1; i < f; i++)
        flimit[i-1] = finger[i];
      flimit[f-1] = cptr;

#define ADD(i)				\
{ int cn = COUNT_OF(finger[i]);		\
  advn[a++] = i;			\
  if (cn <= HAPLO_HGH)			\
    good[c++] = i;			\
}

#define SET(i)	\
{ mc = hc;	\
  mr = hr;	\
  a = c = 0;	\
  ADD(i);	\
}

      while (1)
        { for (i = 0; i < f; i++)
            if (finger[i] < flimit[i])
              break;
          if (i >= f)
            break;
          hr = finger[i]+offs;
          hc = hr[-1] & mask;
          SET(i);
          for (i++; i < f; i++)
            if (finger[i] < flimit[i])
              { hr = finger[i]+offs;
                hc = hr[-1] & mask;
                if (hc == mc)
                  { int v = mycmp(hr,mr,rem);
                    if (v == 0)
                      ADD(i)
                    else if (v < 0)
                      SET(i)
                  }
                else if (hc < mc)
                  SET(i)
              }

#ifdef DEBUG_PARTITION
          { int j;

            j = 0;
            for (i = 0; i < a; i++)
              { printf(" %d",advn[i]);
                if (j < c && advn[i] == good[j])
                  { printf("+");
                    j += 1;
                  }
              }
            printf("\n");
          }
#endif

          if (c > 1) 
            { for (i = 0; i < c; i++)
                { uint8 *fp = finger[good[i]];
                  print_seq(fp,kmer);
                  printf(" %d\n",COUNT_OF(fp));
                }
              printf("\n");
            }
          for (i = 0; i < a; i++)
            finger[advn[i]] += tbyte;

#ifdef DEBUG_PARTITION
          for (i = 0; i < f; i++)
            printf(" %ld",(finger[i]-cache)/tbyte);
          printf("\n");
#endif
        }
    }
}
          

/****************************************************************************************
 *
 *  Main
 *
 *****************************************************************************************/

int main(int argc, char *argv[])
{ Kmer_Stream *T;

  { int    i, j, k;
    int    flags[128];
    char  *eptr, *fptr;

    (void) flags;

    ARG_INIT("Haplex");

    HAPLO_LOW = 1;
    HAPLO_HGH = 0x7fff;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
            break;
          case 'h':
            HAPLO_LOW = strtol(argv[i]+2,&eptr,10);
            if (eptr > argv[i]+2)
              { if (HAPLO_LOW < 1 || HAPLO_LOW > 0x7fff)
                  { fprintf(stderr,"%s: Haplotype minimum count %d is out of range\n",
                                   Prog_Name,HAPLO_LOW);
                    exit (1);
                  }
                if (*eptr == ':')
                  { HAPLO_HGH = strtol(eptr+1,&fptr,10);
                    if (fptr > eptr+1 && *fptr == '\0')
                      { if (HAPLO_HGH < 1 || HAPLO_HGH > 0x7fff)
                          { fprintf(stderr,"%s: Haplotype maximum count %d is out of range\n",
                                           Prog_Name,HAPLO_HGH);
                            exit (1);
                          }
                        if (HAPLO_LOW > HAPLO_HGH)
                          { fprintf(stderr,"%s: Histogram range is invalid\n",Prog_Name);
                            exit (1);
                          }
                        break;
                      }
                  }
              }
            fprintf(stderr,"%s: Syntax of -h option invalid -h<int>:<int>\n",Prog_Name);
            exit (1);
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -h: Accept only haplotypes with count in given range (inclusive).\n");
        exit (1);
      }
  }

  T = Open_Kmer_Stream(argv[1],HAPLO_LOW);

  Find_Haplo_Pairs(T);

  Free_Kmer_Stream(T);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
