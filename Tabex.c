/*********************************************************************************************\
 *
 *  Example code for reading, listing, and searching a kmer-count table produced by FastK.
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

#undef UNTESTED

#include "gene_core.h"

static char *Usage = "[-t<int>] <source_root>[.ktab] (LIST|CHECK|(k-mer:string>) ...";


/****************************************************************************************
 *
 *  The interface you may want to lift for your own use
 *
 *****************************************************************************************/

typedef struct
  { int     kmer;    //  Kmer length
    int     kbyte;   //  Kmer encoding in bytes
    int     tbyte;   //  Kmer+count entry in bytes
    int64   nels;    //  # of unique, sorted k-mers in the table
    uint8  *table;   //  The (huge) table in memory
    int64  *index;   //  Search accelerator if needed
  } Kmer_Table;

#define  KMER(i)  (table+(i)*tbyte)
#define  COUNT(i) (*((uint16 *) (table+(i)*tbyte+kbyte)))

Kmer_Table *Load_Kmer_Table(char *name, int cut_freq, int smer, int nthreads);
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

Kmer_Table *Load_Kmer_Table(char *name, int cut_freq, int smer, int nthreads)
{ Kmer_Table *T;
  int         kmer, tbyte, kbyte;
  int64       nels;
  uint8      *table;

  FILE  *f;
  int    p;
  int64  n;

  //  Find all parts and accumulate total size

  nels = 0;
  for (p = 1; p <= nthreads; p++)
    { f = fopen(Catenate(name,Numbered_Suffix(".ktab.",p,""),"",""),"r");
      if (f == NULL)
        { fprintf(stderr,"%s: Table part %s.ktab.%d is misssing ?\n",Prog_Name,name,p);
          exit (1);
        }
      fread(&kmer,sizeof(int),1,f);
      fread(&n,sizeof(int64),1,f);
      nels += n;
      if (kmer != smer)
        { fprintf(stderr,"%s: Table part %s.ktab.%d does not have k-mer length matching stub ?\n",
                         Prog_Name,name,p);
          exit (1);
        }
      fclose(f);
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
  for (p = 1; p <= nthreads; p++)
    { f = fopen(Catenate(name,Numbered_Suffix(".ktab.",p,""),"",""),"r");
      fread(&kmer,sizeof(int),1,f);
      fread(&n,sizeof(int64),1,f);
      fread(KMER(nels),n*tbyte,1,f);
      nels += n;
      fclose(f);
    }

  if (cut_freq > 1)
    { int64 i, j;

      j = 0;
      for (i = 0; i < nels; i++)
        if (COUNT(i) >= cut_freq)
          mycpy(KMER(j++),KMER(i),tbyte);
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
  T->index = NULL;

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
  
  int i;

  setup_fmer_table();

  printf("\n");
  for (i = 1; i < nels; i++)
    { if (mycmp(KMER(i-1),KMER(i),kbyte) >= 0)
        { printf("Out of Order\n");
          printf(" %9d:",i-1);
          print_pack(KMER(i-1),kmer);
          printf("  ");
          print_seq(KMER(i-1),kmer);
          printf(" = %4d\n",COUNT(i-1));
          printf(" %9d:",i);
          print_pack(KMER(i),kmer);
          printf("  ");
          print_seq(KMER(i),kmer);
          printf(" = %4d\n",COUNT(i));
          break;
        }
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

  int i;

  setup_fmer_table();

  printf("\nElement Bytes = %d  Kmer Bytes = %d\n",tbyte,kbyte);

  printf(" %9d: ",0);
  print_seq(KMER(0),kmer);
  printf(" = %4d\n",COUNT(0));
  for (i = 1; i < nels; i++)
    { if (mycmp(KMER(i-1),KMER(i),kbyte) >= 0)
        printf("Out of Order\n");
      printf(" %9d: ",i);
      print_seq(KMER(i),kmer);
      printf(" = %4d\n",COUNT(i));
    }
}


/****************************************************************************************
 *
 *  Find k-mer in table
 *
 *****************************************************************************************/

void set_up_accelerator(Kmer_Table *T)
{ int     tbyte = T->tbyte;
  int     kbyte = T->kbyte;
  int64   nels  = T->nels;
  uint8  *table = T->table;
  int64  *index;

  uint8 *iptr, *nptr;
  int64  i;
  int    idx, val;

  index = Malloc(sizeof(uint8 *)*0x1000001,"Allocating acceleraator");
  if (index == NULL)
    exit (1);

  idx  = 1;
  iptr = table;
  nptr = KMER(nels);
  for (i = 1, iptr += tbyte; iptr < nptr; i++, iptr += tbyte)
    { if (mycmp(iptr,iptr-tbyte,kbyte) == 0)
        continue;
      val = (iptr[0] << 16) | (iptr[1] << 8) | iptr[0];
      while (idx <= val)
        index[idx++] = i;
    }

  index[0] = 0;
  index[0x1000000] = nels;

  T->index = index;
}

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
  s0[0] = s1[0] = s2[0] = 3;

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

#ifdef UNTESTED
  if (kbyte >= 3)
    { int64 *index = T->index;
      if (index == NULL)
        { set_up_accelerator(T);
          index = T->index;
        }
      m = (cmp[0] << 16) | (cmp[1] << 8) | cmp[2];
      l = index[m];
      r = index[m+1];
    }
  else
#endif

    { l = 0;
      r = nels;
    }

  // smallest l s.t. KMER(l) >= (kmer) cmp  (or nels if does not exist)

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
{ char       *name;
  int         smer, nthreads;
  Kmer_Table *T;
  int         CUT;

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    (void) flags;

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

    if (argc < 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  { FILE *f;
    char *dir, *root;

    dir  = PathTo(argv[1]);
    root = Root(argv[1],".ktab");
    name = Strdup(Catenate(dir,"/.",root,""),NULL);
    f = fopen(Catenate(dir,"/",root,".ktab"),"r");
    if (f == NULL)
      { fprintf(stderr,"%s: Cannot open %s for reading\n",Prog_Name,Catenate(dir,"/",root,".ktab"));
        exit (1);
      }
    fread(&smer,sizeof(int),1,f);
    fread(&nthreads,sizeof(int),1,f);
    fclose(f);
    free(root);
    free(dir);
  }

  T = Load_Kmer_Table(name,CUT,smer,nthreads);

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

  free(name);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
