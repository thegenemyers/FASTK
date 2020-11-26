/*******************************************************************************************
 *
 *  C library routines to access and operate upon FastK histogram, k-mer tables, and profiles
 *
 *  Author:  Gene Myers
 *  Date  :  November 2020
 *
 *******************************************************************************************/

#include "libfastk.h"

#include "gene_core.c"

/*********************************************************************************************\
 *
 *  HISTOGRAM CODE
 *
 *********************************************************************************************/

Histogram *Load_Histogram(char *name)
{ Histogram *H;
  int        kmer, low, high;
  int64     *hist;
  char      *dir, *root;
  int        f;

  dir  = PathTo(name);
  root = Root(name,".hist");
  f = open(Catenate(dir,"/",root,".hist"),O_RDONLY);
  if (f < 0)
    return (NULL);
  free(root);
  free(dir);

  read(f,&kmer,sizeof(int));
  read(f,&low,sizeof(int));
  read(f,&high,sizeof(int));

  H    = Malloc(sizeof(Histogram),"Allocating histogram");
  hist = Malloc(sizeof(int64)*((high-low)+1),"Allocating histogram");
  if (H == NULL || hist == NULL)
    exit (1);

  read(f,hist,sizeof(int64)*((high-low)+1));
    
  close(f);

  H->kmer = kmer;
  H->low  = low;
  H->high = high;
  H->hist = hist-low;

  return (H);
}

void Subrange_Histogram(Histogram *H, int low, int high)
{ int64  under, over;
  int    i;
 
  if (H->low > low || H->high < high)
    return;

  under = 0;
  for (i = H->low; i < low; i++)
    under += H->hist[i];
  over = 0;
  for (i = high+1; i < H->high; i++)
    over += H->hist[i];

  if (low != H->low)
    memmove(H->hist+H->low,H->hist+low,((high-low)+1)*sizeof(int64));

  H->hist += H->low;
  H->hist  = Realloc(H->hist,((high-low)+1)*sizeof(int64),"Reallocating histogram");
  H->hist -= low;

  H->hist[low]  += under;
  H->hist[high] += over;
  H->low  = low;
  H->high = high;
}

void Free_Histogram(Histogram *H)
{ free(H->hist+H->low);
  free(H);
}

/****************************************************************************************
 *
 *  K-MER TABLE CODE
 *
 *****************************************************************************************/

typedef struct
  { int     kmer;         //  Kmer length
    int     minval;       //  The minimum count of a k-mer in the table
    int     kbyte;        //  Kmer encoding in bytes
    int     tbyte;        //  Kmer+count entry in bytes
    int64   nels;         //  # of unique, sorted k-mers in the table
    uint8  *table;        //  The (huge) table in memory
    void   *index;        //  Accelerator index for searches
  } _Kmer_Table;

/****************************************************************************************
 *
 *  Print & compare utilities
 *
 *****************************************************************************************/

#define  KMER(i)  (table+(i)*tbyte)
#define  COUNT(i) (*((uint16 *) (table+(i)*tbyte+kbyte)))
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

static void print_seq(FILE *out, uint8 *seq, int len)
{ int i, b, k;

  b = len >> 2;
  for (i = 0; i < b; i++)
    fprintf(out,"%s",fmer[seq[i]]);
  k = 6;
  for (i = b << 2; i < len; i++)
    { fprintf(out,"%c",dna[seq[b] >> k]);
      k -= 2;
    }
}

static void print_pack(FILE *out, uint8 *seq, int len)
{ int i;

  for (i = 0; i < (len+3)/4; i++)
    fprintf(out," %02x",seq[i]);
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

static inline int64 big_read(int f, uint8 *buffer, int64 bytes)
{ int64 v;

  v = 0;
  while (bytes > 0x70000000)
    { v += read(f,buffer,0x70000000);
      bytes  -= 0x70000000;
      buffer += 0x70000000;
    }
  v += read(f,buffer,bytes);
  return (v);
}

Kmer_Table *Load_Kmer_Table(char *name, int cut_off)
{ Kmer_Table  *T;
  Kmer_Stream *S;
  int          kmer, tbyte, kbyte, minval;
  int64        nels;
  uint8       *table;

  int    f;
  char  *dir, *root, *hidden;
  int    smer, nthreads;

  setup_fmer_table();

  //  Open stub file and get # of parts

  dir    = PathTo(name);
  root   = Root(name,".ktab");
  hidden = Strdup(Catenate(dir,"/.",root,""),NULL);
  f = open(Catenate(dir,"/",root,".ktab"),O_RDONLY);
  free(root);
  free(dir);
  if (f < 0)
    return (NULL);
  read(f,&smer,sizeof(int));
  read(f,&nthreads,sizeof(int));
  read(f,&minval,sizeof(int));
  close(f);

  kmer  = smer;
  kbyte = (kmer+3)>>2;
  tbyte = kbyte+2;

  //  Find all parts and accumulate total size

  nels = 0;
  if (cut_off > minval)

    { uint8 *iptr;

      S = Open_Kmer_Stream(name,cut_off);
      for (iptr = First_Kmer_Entry(S); iptr != NULL; iptr = Next_Kmer_Entry(S))
        nels += 1;
    }

  else

    { int    p;
      int64  n;
 
      for (p = 1; p <= nthreads; p++)
        { f = open(Catenate(hidden,Numbered_Suffix(".ktab.",p,""),"",""),O_RDONLY);
          if (f < 0)
            { fprintf(stderr,"Table part %s.ktab.%d is misssing ?\n",hidden,p);
              exit (1);
            }
          read(f,&kmer,sizeof(int));
          read(f,&n,sizeof(int64));
          nels += n;
          if (kmer != smer)
            { fprintf(stderr,"Table part %s.ktab.%d does not have k-mer length matching stub ?\n",
                             hidden,p);
              exit (1);
            }
          close(f);
        }
    }

  //  Allocate in-memory table

  T     = Malloc(sizeof(Kmer_Table),"Allocating table record");
  table = Malloc(nels*tbyte,"Allocating k-mer table\n");
  if ( T == NULL || table == NULL)
    exit (1);

  //  Load the parts into memory

  if (cut_off > minval)

    { uint8 *iptr, *jptr;

      jptr = table;
      for (iptr = First_Kmer_Entry(S); iptr != NULL; iptr = Next_Kmer_Entry(S))
        { mycpy(jptr,iptr,tbyte);
          jptr += tbyte;
        }
      Free_Kmer_Stream(S);

      minval = cut_off;
    }

  else

    { int    p;
      int64  n;
 
      nels = 0;
      for (p = 1; p <= nthreads; p++)
        { f = open(Catenate(hidden,Numbered_Suffix(".ktab.",p,""),"",""),O_RDONLY);
          read(f,&kmer,sizeof(int));
          read(f,&n,sizeof(int64));
          big_read(f,KMER(nels),n*tbyte);
          nels += n;
          close(f);
        }
    }

  free(hidden);

  T->kmer   = kmer;
  T->minval = minval;
  T->tbyte  = tbyte;
  T->kbyte  = kbyte;
  T->nels   = nels;
  T->table  = table;
  ((_Kmer_Table *) T)->index = NULL;

  return (T);
}


/****************************************************************************************
 *
 *  Free, Fetch, Check, and List a Table
 *
 *****************************************************************************************/

char *Fetch_Kmer(Kmer_Table *T, int i)
{ static char *seq = NULL;
  static int   max = 0;

  int    kmer  = T->kmer;
  int    tbyte = T->tbyte;
  uint8 *table = T->table;

  if (T == NULL)
    { free(seq);
      max = 0;
      return (NULL);
    }

  if (kmer > max)
    { max = kmer;
      seq = Realloc(seq,max+1,"Reallocating k-mer buffer");
      if (seq == NULL)
        exit (1);
    }

  { int    j, k;
    uint8 *a;

    a = KMER(i);
    for (j = 0; j < kmer; j += 4)
      sprintf(seq+j,"%s",fmer[*a++]);
    k = 6;
    for (j -= 4; j < kmer; j++)
      { seq[j] = dna[*a >> k];
        k -= 2;
      }
  }

  return (seq);
}

int Fetch_Count(Kmer_Table *T, int i)
{ int    kbyte = T->kbyte;
  int    tbyte = T->tbyte;
  uint8 *table = T->table;

  return (COUNT(i));
}

void Free_Kmer_Table(Kmer_Table *T)
{ free(T->table);
  free(((_Kmer_Table *) T)->index);
  free(T);
}

int Check_Kmer_Table(Kmer_Table *T)
{ int    kmer  = T->kmer;
  int    tbyte = T->tbyte;
  int    kbyte = T->kbyte;
  int64  nels  = T->nels;
  uint8 *table = T->table;
  
  int i;

  for (i = 1; i < nels; i++)
    { if (mycmp(KMER(i-1),KMER(i),kbyte) >= 0)
        { fprintf(stderr,"\nOut of Order\n");
          fprintf(stderr," %9d:",i-1);
          print_pack(stderr,KMER(i-1),kmer);
          fprintf(stderr,"  ");
          print_seq(stderr,KMER(i-1),kmer);
          fprintf(stderr," = %4d\n",COUNT(i-1));
          fprintf(stderr," %9d:",i);
          print_pack(stderr,KMER(i),kmer);
          fprintf(stderr,"  ");
          print_seq(stderr,KMER(i),kmer);
          fprintf(stderr," = %4d\n",COUNT(i));
          return (0);
        }
    }
  return (1);
}

void List_Kmer_Table(Kmer_Table *T, FILE *out)
{ int    kmer  = T->kmer;
  int    tbyte = T->tbyte;
  int    kbyte = T->kbyte;
  int64  nels  = T->nels;
  uint8 *table = T->table;

  int i;

  fprintf(out,"\nElement Bytes = %d  Kmer Bytes = %d\n",tbyte,kbyte);

  fprintf(out," %9d: ",0);
  print_seq(out,KMER(0),kmer);
  fprintf(out," = %4d\n",COUNT(0));
  for (i = 1; i < nels; i++)
    { if (mycmp(KMER(i-1),KMER(i),kbyte) >= 0)
        fprintf(out,"Out of Order\n");
      fprintf(out," %9d: ",i);
      print_seq(out,KMER(i),kmer);
      fprintf(out," = %4d\n",COUNT(i));
    }
}


/****************************************************************************************
 *
 *  Find k-mer in table
 *
 *****************************************************************************************/

static void set_up_accelerator(_Kmer_Table *T)
{ int     tbyte = T->tbyte;
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
  index[0] = 0;
  for (i = 1, iptr += tbyte; iptr < nptr; i++, iptr += tbyte)
    { if (mycmp(iptr,iptr-tbyte,3) == 0)
        continue;
      val = (iptr[0] << 16) | (iptr[1] << 8) | iptr[2];
      while (idx <= val)
        index[idx++] = i;
    }
  while (idx <= 0x1000000)
    index[idx++] = nels;

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

static void compress_norm(char *s, int len, uint8 *t)
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
    compress_norm(kseq,kmer,cmp);
  else
    compress_comp(kseq,kmer,cmp);

  if (kbyte >= 3)
    { int64 *index = ((_Kmer_Table *) T)->index;
      if (index == NULL)
        { set_up_accelerator((_Kmer_Table *) T);
          index = ((_Kmer_Table *) T)->index;
        }
      m = (cmp[0] << 16) | (cmp[1] << 8) | cmp[2];
      l = index[m];
      r = index[m+1];
    }
  else
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
 *  K-MER STREAM CODE
 *
 *****************************************************************************************/

typedef struct
  { int    kmer;    //  Kmer length
    int    minval;  //  The minimum count of a k-mer in the stream
    int    kbyte;   //  Kmer encoding in bytes
    int    tbyte;   //  Kmer+count entry in bytes
    int64  nels;    //  # of elements in entire table

    int    copn;    //  File currently open
    int    part;    //  Thread # of file currently open
    int    nthr;    //  # of thread parts
    int    cut;     //  need to cut by minval
    char  *name;    //  Root name for table
    uint8 *table;   //  The (huge) table in memory
    uint8 *celm;    //  Current element
    uint8 *ctop;    //  Ptr top of current table block
  } _Kmer_Stream;

/****************************************************************************************
 *
 *  Open a table and return as Kmer_Stream object
 *
 *****************************************************************************************/

static void More_Kmer_Stream(_Kmer_Stream *S)
{ int    tbyte = S->tbyte;
  uint8 *table = S->table;
  int    copn  = S->copn;
  uint8 *ctop;
  int64  rels;
  int    kmer;

  while (1)
    { ctop = table + read(copn,table,1024*tbyte);
      if (ctop <= table)
        close(copn);
      else if (S->cut)
        { int    minv  = S->minval;
          int    kbyte = S->kbyte;
          uint8 *iptr, *jptr;

          jptr = table;
          for (iptr = table; iptr < ctop; iptr += tbyte)
            { if (COUNT_OF(iptr) >= minv)
                { mycpy(jptr,iptr,tbyte);
                  jptr += tbyte;
                }
            }
          ctop = jptr;
          if (ctop <= table)
            continue;
          else
            break;
        }
      else
        break;
      S->part += 1;
      if (S->part > S->nthr)
        { S->celm = NULL;
          return;
        }
      copn = open(Catenate(S->name,Numbered_Suffix(".ktab.",S->part,""),"",""),O_RDONLY);
      read(copn,&kmer,sizeof(int));
      read(copn,&rels,sizeof(int64));
    }
  S->celm = table;
  S->ctop = ctop;
  S->copn = copn;
}

Kmer_Stream *Open_Kmer_Stream(char *name, int cut_off)
{ _Kmer_Stream *S;
  int           kmer, tbyte, kbyte, minval;
  int64         nels;
  int           copn;

  int    f;
  char  *dir, *root, *hidden;
  int    smer, nthreads;
  int    p;
  int64  n;

  setup_fmer_table();

  //  Open stub file and get # of parts

  dir    = PathTo(name);
  root   = Root(name,".ktab");
  hidden = Strdup(Catenate(dir,"/.",root,""),NULL);
  f = open(Catenate(dir,"/",root,".ktab"),O_RDONLY);
  free(root);
  free(dir);
  if (f < 0)
    return (NULL);
  read(f,&smer,sizeof(int));
  read(f,&nthreads,sizeof(int));
  read(f,&minval,sizeof(int));
  close(f);

  kbyte   = (smer+3)>>2;
  tbyte    = kbyte+2;

  //  Find all parts and accumulate total size

  nels = 0;
  for (p = 1; p <= nthreads; p++)
    { copn = open(Catenate(hidden,Numbered_Suffix(".ktab.",p,""),"",""),O_RDONLY);
      if (copn < 0)
        { fprintf(stderr,"%s: Table part %s.ktab.%d is misssing ?\n",Prog_Name,hidden,p);
          exit (1);
        }
      read(copn,&kmer,sizeof(int));
      read(copn,&n,sizeof(int64));
      nels += n;
      if (kmer != smer)
        { fprintf(stderr,"%s: Table part %s.ktab.%d does not have k-mer length matching stub ?\n",
                         Prog_Name,hidden,p);
          exit (1);
        }
      close(copn);
    }

  //  Allocate in-memory buffer & establish initial condition

  S        = Malloc(sizeof(Kmer_Stream),"Allocating table record");
  S->name  = Strdup(hidden,"Allocating name");
  S->table = Malloc(1024ll*tbyte,"Allocating k-mer table\n");
  if (S == NULL || S->name == NULL || S->table == NULL)
    exit (1);

  free(hidden);

  S->kmer   = kmer;
  if (minval < cut_off)
    S->minval = cut_off;
  else
    S->minval = minval;
  S->tbyte  = tbyte;
  S->kbyte  = kbyte;
  S->nels   = nels;

  copn = open(Catenate(hidden,".ktab.1","",""),O_RDONLY);
  read(copn,&kmer,sizeof(int));
  read(copn,&n,sizeof(int64));

  S->copn  = copn;
  S->part  = 1;
  S->nthr  = nthreads;
  S->cut   = (minval < cut_off);

  More_Kmer_Stream(S);

  return ((Kmer_Stream *) S);
}

/****************************************************************************************
 *
 *  Free, Iterate, and Get Stream Entries
 *
 *****************************************************************************************/

inline uint8 *First_Kmer_Entry(Kmer_Stream *_S)
{ _Kmer_Stream *S = (_Kmer_Stream *) _S;

  if (S->part != 1)
    { if (S->part <= S->nthr)
        close(S->copn);
      S->copn = open(Catenate(S->name,".ktab.1","",""),O_RDONLY);
      S->part = 1;
    }

  lseek(S->copn,sizeof(int)+sizeof(int64),SEEK_SET);

  More_Kmer_Stream(S);

  return (S->celm);
}

#define REAL(S) ((_Kmer_Stream *) S)

inline uint8 *Next_Kmer_Entry(Kmer_Stream *S)
{ REAL(S)->celm += REAL(S)->tbyte;
  if (REAL(S)->celm >= REAL(S)->ctop)
    More_Kmer_Stream(REAL(S));
  return (REAL(S)->celm);
}

char *Current_Kmer(Kmer_Stream *S)
{ static char *seq = NULL;
  static int   max = 0;

  int kmer = S->kmer;

  if (S == NULL)
    { free(seq);
      max = 0;
      return (NULL);
    }

  if (kmer > max)
    { max = kmer;
      seq = Realloc(seq,max+1,"Reallocating k-mer buffer");
      if (seq == NULL)
        exit (1);
    }

  { int    j, k;
    uint8 *a;

    a = ((_Kmer_Stream *) S)->celm;
    for (j = 0; j < kmer; j += 4)
      sprintf(seq+j,"%s",fmer[*a++]);
    k = 6;
    for (j -= 4; j < kmer; j++)
      { seq[j] = dna[*a >> k];
        k -= 2;
      }
  }

  return (seq);
}

int Current_Count(Kmer_Stream *S)
{ int kbyte = S->kbyte;
  return (COUNT_OF(((_Kmer_Stream *) S)->celm));
}

void Free_Kmer_Stream(Kmer_Stream *_S)
{ _Kmer_Stream *S = (_Kmer_Stream *) _S;

  free(S->name);
  free(S->table);
  if (S->copn >= 0)
    close(S->copn);
  free(_S);
}


/*********************************************************************************************\
 *
 *  PROFILE CODE
 *
 *****************************************************************************************/

/****************************************************************************************
 *
 *  Open a profile as a Profile_Index.  Index to compressed profiles is in memory,
 *    but compressed profiles are left on disk and reaad only when requested.
 *
 *****************************************************************************************/

Profile_Index *Open_Profiles(char *name)
{ Profile_Index *P;
  int            kmer, nparts;
  int64          nreads, *nbase, *index;
  int           *nfile;

  int    f;
  char  *dir, *root, *hidden;
  int    smer, nthreads;
  int64  n;

  //  Open stub file and get # of parts

  dir    = PathTo(name);
  root   = Root(name,".prof");
  hidden = Strdup(Catenate(dir,"/.",root,""),NULL);
  f = open(Catenate(dir,"/",root,".prof"),O_RDONLY);
  free(root);
  free(dir);
  if (f < 0)
    return (NULL);
  read(f,&smer,sizeof(int));
  read(f,&nthreads,sizeof(int));
  close(f);

  //  Find all parts and accumulate total size

  nreads = 0;
  for (nparts = 0; nparts < nthreads; nparts++)
    { f = open(Catenate(hidden,Numbered_Suffix(".pidx.",nparts+1,""),"",""),O_RDONLY);
      if (f < 0)
        { fprintf(stderr,"Profile part %s.pidx.%d is misssing ?\n",hidden,nparts+1);
          exit (1);
        }
      read(f,&kmer,sizeof(int));
      read(f,&n,sizeof(int64));
      read(f,&n,sizeof(int64));
      nreads += n;
      if (kmer != smer)
        { fprintf(stderr,"Profile part %s.pidx.%d does not have k-mer length matching stub ?\n",
                         hidden,nparts+1);
          exit (1);
        }
      close(f);
    }

  //  Allocate in-memory table

  P     = Malloc(sizeof(Profile_Index),"Allocating profile record");
  index = Malloc((nreads+1)*sizeof(int64),"Allocating profile index");
  nbase = Malloc(nparts*sizeof(int64),"Allocating profile index");
  nfile = Malloc(nparts*sizeof(FILE *),"Allocating profile index");
  if (P == NULL || index == NULL || nbase == NULL || nfile == NULL)
    exit (1);

  nreads = 0;
  index[0] = 0;
  for (nparts = 0; nparts < nthreads; nparts++)
    { f = open(Catenate(hidden,Numbered_Suffix(".pidx.",nparts+1,""),"",""),O_RDONLY);
      read(f,&kmer,sizeof(int));
      read(f,&n,sizeof(int64));
      read(f,&n,sizeof(int64));
      read(f,index+(nreads+1),n*sizeof(int64));
      nreads += n;
      nbase[nparts] = nreads;
      close(f);

      f = open(Catenate(hidden,Numbered_Suffix(".prof.",nparts+1,""),"",""),O_RDONLY);
      if (f < 0)
        { fprintf(stderr,"Profile part %s.prof.%d is misssing ?\n",hidden,nparts+1);
          exit (1);
        }
      nfile[nparts] = f;
    }

  free(hidden);

  P->kmer   = kmer;
  P->nparts = nparts;
  P->nreads = nreads;
  P->index  = index;
  P->nbase  = nbase;
  P->nfile  = nfile;

  return (P);
}

/****************************************************************************************
 *
 *  Free a Profile_Index and fetch a profile
 *
 *****************************************************************************************/

#undef SHOW_RUN

void Free_Profiles(Profile_Index *P)
{ int i;

  free(P->index);
  free(P->nbase);
  for (i = 0; i < P->nparts; i++)
    close(P->nfile[i]);
  free(P->nfile);
  free(P);
}

  //  Places uncompressed profile for read id (0-based) in profile of length plen.
  //    Returns the length of the uncompressed profile.  If the plen is less than
  //    this then only the first plen counts are uncompressed into profile

#define PROF_BUF0 4096
#define PROF_BUF1 4095

int Fetch_Profile(Profile_Index *P, int64 id, int plen, uint16 *profile)
{ uint8 count[PROF_BUF0], *cend = count+PROF_BUF1;
  int    f;
  int    w, len;
  uint8 *p, *q;
  uint16 x, d, i;
  int    n;

  for (w = 0; w < P->nparts; w++)
    if (id < P->nbase[w])
      break;
  if (w >= P->nparts)
    { fprintf(stderr,"Id %lld is out of range [1,%lld]\n",id,P->nbase[P->nparts-1]);
      exit (1);
    }
  f = P->nfile[w];

  if (id == 0 || (w > 0 && id == P->nbase[w-1]))
    { lseek(f,0,SEEK_SET);
      len = P->index[id+1];
    }
  else
    { int64 off = P->index[id];
      lseek(f,off,SEEK_SET);
      len = P->index[id+1] - off;
    }

  if (len == 0)
    return (len);

  read(f,count,PROF_BUF0);

  p = count;
  q = count + len;

  x = *p++;
  if ((x & 0x80) != 0)
    d = ((x & 0x7f) << 8) | *p++;
  else
    d = x;
  n = 1;

  if (plen > 0)
    { profile[0] = d;
#ifdef SHOW_RUN
  printf(" %d\n",d);
#endif

      while (p < q)
        { if (p >= cend)
            { if (p == cend)
                { *count = *p; 
                  read(f,count+1,PROF_BUF1);
                  q -= PROF_BUF1;
                }
              else
                { read(f,count,PROF_BUF0);
                  q -= PROF_BUF0;
                }
              p = count;
            }
          x = *p++;
          if ((x & 0xc0) == 0)
            { if (n+x > plen)
                { n += x;
                  break;
                }
              for (i = 0; i < x; i++)
                profile[n++] = d;
#ifdef SHOW_RUN
              printf(" [%hu]\n",x);
#endif
            }
          else
            { if ((x & 0x80) != 0)
                { if ((x & 0x40) != 0)
                    x <<= 8;
                  else
                    x = (x << 8) & 0x7fff;
                  x |= *p++;
                  d = (d+x) & 0x7fff;
#ifdef SHOW_RUN
                  printf(" %hd+(%d)\n",x,d);
#endif
                }
              else
                { if ((x & 0x20) != 0)
                    d += (x & 0x1fu) | 0xffe0u;
                  else
                    d += (x & 0x1fu);
#ifdef SHOW_RUN
                  if ((x & 0x20) != 0)
                    printf(" -%d(%d)\n",32-(x&0x1fu),d);
                  else
                    printf(" +%d(%d)\n",x&0x1fu,d);
#endif
                }
              if (n >= plen)
                { n += 1;
                  break;
                }
              profile[n++] = d;
            }
        }
    }

  while (p < q)
    { if (p >= cend)
        { if (p == cend)
            { *count = *p; 
              read(f,count+1,PROF_BUF1);
              q -= PROF_BUF1;
            }
          else
            { read(f,count,PROF_BUF0);
              q -= PROF_BUF0;
            }
          p = count;
        } 
      x = *p++;
      if ((x & 0xc0) == 0)
        n += x;
      else
        { if ((x & 0x80) != 0)
            p += 1;
          n += 1;
        }
    }

  return (n);
}
