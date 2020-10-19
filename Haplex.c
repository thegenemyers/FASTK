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

#include "gene_core.h"

static char *Usage = "[-t<int>] <source_root>.K<k> (LIST|CHECK|(k-mer:string>) ...";


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
{ while (n--)
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

  printf("Loading %d-mer table with ",kmer);
  Print_Number(nels,0,stdout);
  printf(" entries in %d parts\n",p-1);
  fflush(stdout);

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


/****************************************************************************************
 *
 *  Free, Check, and List a Table
 *
 *****************************************************************************************/

void Free_Kmer_Table(Kmer_Table *T)
{ free(T->table);
  free(T);
}

void Check_Table(Kmer_Table *T)
{ int    kmer  = T->kmer;
  int    tbyte = T->tbyte;
  int    kbyte = T->kbyte;
  int64  nels  = T->nels;
  uint8 *table = T->table;
  
  int    i;
  uint8 *iptr;

  setup_fmer_table();

  printf("\n");
  iptr = KMER(0);
  for (i = 1; i < nels; i++)
    { if (mycmp(iptr,iptr+tbyte,kbyte) >= 0)
        { printf("Out of Order\n");
          printf(" %9d:",i-1);
          print_pack(iptr,kmer);
          printf("  ");
          print_seq(iptr,kmer);
          printf(" = %4d\n",COUNT_OF(iptr));
          iptr += tbyte;
          printf(" %9d:",i);
          print_pack(iptr,kmer);
          printf("  ");
          print_seq(iptr,kmer);
          printf(" = %4d\n",COUNT_OF(iptr));
          break;
        }
      iptr += tbyte;
    }
  if (i >= nels)
    printf("Table is OK\n");
}

void List_Table(Kmer_Table *T)
{ int    kmer  = T->kmer;
  int    tbyte = T->tbyte;
  int    kbyte = T->kbyte;
  int64  nels  = T->nels;
  uint8 *table = T->table;

  int    i;
  uint8 *iptr;

  setup_fmer_table();

  printf("\nElement Bytes = %d  Kmer Bytes = %d\n",tbyte,kbyte);

  printf(" %9d: ",0);
  iptr = KMER(0);
  print_seq(iptr,kmer);
  printf(" = %4d\n",COUNT_OF(iptr));
  for (i = 1; i < nels; i++)
    { iptr += tbyte;
      if (mycmp(iptr-tbyte,iptr,kbyte) >= 0)
        printf("Out of Order\n");
      printf(" %9d: ",i);
      print_seq(iptr,kmer);
      printf(" = %4d\n",COUNT_OF(iptr));
    }
}

static inline int mypref(uint8 *a, uint8 *b, int n)
{ int   i;
  uint8 x, y;
  
  for (i = 0; i <= n; i += 4)
    { if (*a != *b)
        { x = *a;
          y = *b;
          if ((x & 0xf0) == (y & 0xf0))
            if ((x & 0xc0) == (y & 0xc0))
              return (i + 1);
            else
              return (i);
          else
            if ((x & 0xfc) == (y & 0xfc))
              return (i + 3);
            else
              return (i + 2);
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
  uint8 *iptr, *jptr, *nptr;
  uint8 *finger[5];
  int    n, i, c, f, x, v;
  int    mc, hc;
  uint8 *mr, *hr;
  int    mask, offs, rem;

  khalf = kmer/2;

  rem = offs = mask = 0;  //  FIX

  nptr = KMER(nels);
  for (iptr = table; iptr < nptr; iptr = jptr)
    { f = 0;
      finger[f++] = iptr;
      for (jptr = iptr+tbyte; jptr < nptr; jptr += tbyte)
        { x = mypref(jptr-tbyte,jptr,khalf); 
          if (x < khalf)
            break;
          if (x == khalf)
            finger[f++] = jptr;
        }
      if (f <= 1)
        continue;

      finger[f] = jptr;
      for (n = (jptr-iptr)/tbyte; n > 0; n--)
        { mr = finger[0]+offs;
          mc = mr[-1] & mask;
          c  = 1;
          x  = 0;
          for (i = 1; i < f; i++)
            if (finger[i] < finger[i+1])
              { hr = finger[i]+offs;
                hc = hr[-1] & mask;
                if (hc == mc)
                  { v = mycmp(hr,mr,rem);
                    if (v == 0)
                      c += 1;
                    else if (v < 0)
                      { mc = hc;
                        mr = hr;
                        c = 1;
                        x = i;
                      }
                  }
                else if (hc < mc)
                  { mc = hc;
                    mr = hr;
                    c = 1;
                    x = i;
                  }
              }
          if (c > 1)
            { for (i = 0; i < f; i++)
                if (finger[i] < finger[i+1])
                  { hr = finger[i]+offs;
                    hc = hr[-1] & mask;
                    if (hc == mc && mycmp(hr,mr,rem) == 0)
                      { if (c > 1)
                          { print_seq(finger[i],kmer);
                            printf(" %d\n",COUNT_OF(finger[i]));
                          }
                        finger[i] += tbyte;
                      }
                  }
                printf("\n");
              }
          else
            finger[x] += tbyte;
        }
    }
}
          
      

/****************************************************************************************
 *
 *  Find k-mer in table
 *
 *****************************************************************************************/

static uint8 code[128] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

static uint8 comp[128] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 3, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 3, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

static int is_minimal(char *seq, int len)
{ int j, k;
  int x, y;
  
  for (k = 0, j = len-1; k < j; k++, j--)
    { x = code[(int) seq[k]];
      y = comp[(int) seq[j]];
      if (x < y)
        return (1);
      if (x > y)
        return (0);
    }
  if (k <= j)
    { x = code[(int) seq[k]];
      if (x < 2)
        return (1);
      else
        return (0);
    }
  else
    return (1);
}

static void compress(char *s, int len, uint8 *t)
{ int    i;
  char   c, d, e;
  char  *s0, *s1, *s2, *s3;

  s0 = s;
  s1 = s0+1;
  s2 = s1+1;
  s3 = s2+1;

  c = s0[len];
  d = s1[len];
  e = s2[len];
  s0[len] = s1[len] = s2[len] = 0;

  for (i = 0; i < len; i += 4)
    *t++ = ((code[(int) s0[i]] << 6) | (code[(int) s1[i]] << 4)
         |  (code[(int) s2[i]] << 2) | code[(int) s3[i]] );

  s0[len] = c;
  s1[len] = d;
  s2[len] = e;
}

static void compress_comp(char *s, int len, uint8 *t)
{ int    i;
  char   c, d, e;
  char  *s0, *s1, *s2, *s3;

  s0 = s;
  s1 = s0-1;
  s2 = s1-1;
  s3 = s2-1;

  c = s0[0];
  d = s1[0];
  e = s2[0];
  s0[len] = s1[len] = s2[len] = 3;

  for (i = len-1; i >= 0; i -= 4)
    *t++ = ((comp[(int) s0[i]] << 6) | (comp[(int) s1[i]] << 4)
         |  (comp[(int) s2[i]] << 2) | comp[(int) s3[i]] );

  s0[0] = c;
  s1[0] = d;
  s2[0] = e;
}

int Find_Kmer(Kmer_Table *T, char *kseq)
{ int    kmer  = T->kmer;
  int    tbyte = T->tbyte;
  int    kbyte = T->kbyte;
  int64  nels  = T->nels;
  uint8 *table = T->table;

  uint8  cmp[kbyte];
  int64  l, r, m;

  //  kseq must be at least kmer bp long

  if (is_minimal(kseq,kmer))
    compress(kseq,kmer,cmp);
  else
    compress_comp(kseq,kmer,cmp);

  // smallest l s.t. KMER(l) >= (kmer) cmp  (or nels if does not exist)

  l = 0;
  r = nels;
  while (l < r)
    { m = ((l+r) >> 1);
      if (mycmp(KMER(m),cmp,kbyte) < 0)
        l = m+1;
      else
        r = m;
    }

  if (l >= nels || mycmp(KMER(l),cmp,kbyte) != 0)
    return (0);

  return (COUNT(l));
}


/****************************************************************************************
 *
 *  Test Stub
 *
 *****************************************************************************************/

int main(int argc, char *argv[])
{ Kmer_Table *T;
  int         CUT;

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("Tabex");

    CUT = 1;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
            break;
          case 't':
            ARG_POSITIVE(CUT,"Cutoff for k-mer table")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc != 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  T = Load_Kmer_Table(argv[1],CUT);

  { int c, cnt;

    for (c = 2; c < argc; c++)
      if (strcmp(argv[c],"LIST") == 0)
        List_Table(T);
      else if (strcmp(argv[c],"CHECK") == 0)
        Check_Table(T);
      else
        { if ((int) strlen(argv[c]) != T->kmer)
            printf("%*s: Not a %d-mer\n",T->kmer,argv[c],T->kmer);
          else
            { cnt = Find_Kmer(T,argv[c]);
              if (cnt == 0)
                printf("%*s: Not found\n",T->kmer,argv[c]);
              else
                printf("%*s: %5d\n",T->kmer,argv[c],cnt);
            }
        }
  }

  Free_Kmer_Table(T);
  exit (0);
}
