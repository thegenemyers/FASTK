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

#include "gene_core.h"

static char *Usage = "[-h<int>:<int>] <source_root>.K<k>";


/****************************************************************************************
 *
 *  The interface you may want to lift for your own use
 *
 *****************************************************************************************/

typedef struct
  { int    kmer;    //  Kmer length
    int    kbyte;   //  Kmer encoding in bytes
    int    tbyte;   //  Kmer+count entry in bytes
    int64  nels;    //  # of unique, sorted k-mers in the table
    uint8 *table;   //  The (huge) table in memory
  } Kmer_Table;

#define  KMER(i)  (table+(i)*tbyte)
#define  COUNT(i) (*((uint16 *) (table+(i)*tbyte+kbyte)))
#define  COUNT_OF(p) (*((uint16 *) (p+kbyte)))

Kmer_Table *Load_Kmer_Table(char *name, int cut_freq);
void        Free_Kmer_Table(Kmer_Table *T);

void        List_Kmer_Table(Kmer_Table *T);
void        Check_Kmer_Table(Kmer_Table *T);
int         Find_Kmer(Kmer_Table *T, char *kseq);


/****************************************************************************************
 *
 *  Print & compare utilities
 *
 *****************************************************************************************/

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

static void print_pack(uint8 *seq, int len)
{ int i;

  for (i = 0; i < (len+3)/4; i++)
    printf(" %02x",seq[i]);
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
 *  Read in a table and return as Kmer_Table object
 *
 *****************************************************************************************/

Kmer_Table *Load_Kmer_Table(char *name, int cut_freq)
{ Kmer_Table *T;
  int         kmer, tbyte, kbyte;
  int64       nels;
  uint8      *table;

  FILE  *f;
  int    p;
  int64  n;

  //  Find all parts and accumulate total size

  nels = 0;
  for (p = 1; 1; p++)
    { f = fopen(Catenate(name,Numbered_Suffix(".T",p,""),"",""),"r");
      if (f == NULL)
        break;
      fread(&kmer,sizeof(int),1,f);
      fread(&n,sizeof(int64),1,f);
      nels += n;
      fclose(f);
    }
  if (p == 1)
    { fprintf(stderr,"%s: Cannot find table files for %s\n",Prog_Name,name);
      exit (1);
    }

  //  Allocate in-memory table

  kbyte = (kmer+3)>>2;
  tbyte = kbyte+2;
  table = Malloc(nels*tbyte,"Allocating k-mer table\n");
  if (table == NULL)
    exit (1);

  //  Load the parts into memory

  fprintf(stderr,"Loading %d-mer table with ",kmer);
  Print_Number(nels,0,stderr);
  fprintf(stderr," entries in %d parts\n",p-1);
  fflush(stderr);

  nels = 0;
  for (p = 1; 1; p++)
    { f = fopen(Catenate(name,Numbered_Suffix(".T",p,""),"",""),"r");
      if (f == NULL)
        break;
      fread(&kmer,sizeof(int),1,f);
      fread(&n,sizeof(int64),1,f);
      fread(KMER(nels),n*tbyte,1,f);
      nels += n;
      fclose(f);
    }

  if (cut_freq > 1)
    { int64 i, j;
      uint8 *iptr, *jptr;

      jptr = table;
      for (i = 0; i < nels; i++)
        { iptr = KMER(i);
          if (COUNT_OF(iptr) >= cut_freq)
            { mycpy(jptr,iptr,tbyte);
              jptr += tbyte;
            }
        }
      j = (jptr-table)/tbyte;
      if (j < nels)
        { nels = j;
          table = Realloc(table,nels*tbyte,"Reallocating table");
        }
    }

  T = Malloc(sizeof(Kmer_Table),"Allocating table record");
  if (T == NULL)
    exit (1);

  T->kmer  = kmer;
  T->tbyte = tbyte;
  T->kbyte = kbyte;
  T->nels  = nels;
  T->table = table;

  return (T);
}

void Free_Kmer_Table(Kmer_Table *T)
{ free(T->table);
  free(T);
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

void Find_Haplo_Pairs(Kmer_Table *T)
{ int    kmer  = T->kmer;
  int    tbyte = T->tbyte;
  int    kbyte = T->kbyte;
  int64  nels  = T->nels;
  uint8 *table = T->table;

  int    khalf;
  uint8  prefs[] = { 0x3f, 0x0f, 0x03, 0x00 };
  int    mask, offs, rem;

  uint8 *iptr, *nptr;

  int    f, i;
  uint8 *finger[5];
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

  nptr = KMER(nels);
  iptr = table + 368645*tbyte;
  while (iptr < nptr)
    { f = 0;
      finger[f++] = iptr;
      for (iptr += tbyte; iptr < nptr; iptr += tbyte)
        { int x = mypref(iptr-tbyte,iptr,khalf); 
          if (x < khalf)
            break;
          if (x == khalf)
            finger[f++] = iptr;
        }

#ifdef DEBUG_PARTITION
      printf("part %d",f);
      for (i = 0; i < f; i++)
        printf(" %ld",(finger[i]-table)/tbyte);
      printf(" %ld\n",(iptr-table)/tbyte);
#endif

      if (f <= 1)
        continue;
      finger[f] = iptr;
      for (i = 0; i < f; i++)
        flimit[i] = finger[i+1];

#define ADD(i)					\
{ int cn = COUNT_OF(finger[i]);			\
  advn[a++] = i;				\
  if (HAPLO_LOW <= cn && cn <= HAPLO_HGH)	\
    good[c++] = i;				\
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
            printf(" %ld",(finger[i]-table)/tbyte);
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
{ Kmer_Table *T;

  (void) print_pack;

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
            fprintf(stderr,"%s: Syntax of -h option invalid -h[<int(1)>:]<int>\n",Prog_Name);
            exit (1);
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, output statistics/progress.\n");
        fprintf(stderr,"      -h: Accept only haplotypes with count in given range (inclusive).\n");
        exit (1);
      }
  }

  T = Load_Kmer_Table(argv[1],1);

  Find_Haplo_Pairs(T);

  Free_Kmer_Table(T);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
