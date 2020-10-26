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
#include <math.h>
#include <fcntl.h>

#undef DEBUG_PARTITION

#include "gene_core.h"

static char *Usage = "[-h[<int(1)>:]<int(100)>] <source_root1>.K<k> <source_root2>.K<k> ...";


/****************************************************************************************
 *
 *  The interface you may want to lift for your own use
 *
 *****************************************************************************************/

typedef struct
  { FILE  *copn;    //  File currently open
    int    part;    //  Thread # of file currently open
    char  *name;    //  Root name for table
    int    kmer;    //  Kmer length
    int    kbyte;   //  Kmer encoding in bytes
    int    tbyte;   //  Kmer+count entry in bytes
    int    nels;    //  # of unique, sorted k-mers in the table
    uint8 *table;   //  The (huge) table in memory
  } Kmer_Table;

#define  KMER(i)  (table+(i)*tbyte)
#define  COUNT(i) (*((uint16 *) (table+(i)*tbyte+kbyte)))
#define  COUNT_OF(p) (*((uint16 *) (p+kbyte)))

Kmer_Table *Open_Kmer_Table(char *name, int cut_freq);
int         More_Kmer_Table(Kmer_Table *T);
void        Free_Kmer_Table(Kmer_Table *T);


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

Kmer_Table *Open_Kmer_Table(char *name, int cut_freq)
{ Kmer_Table *T;
  int         kmer, tbyte, kbyte;
  int64       nels;
  uint8      *table;
  FILE       *copn;
  int         part;
  int64  n;

  //  Find all parts and accumulate total size

  nels = 0;
  for (part = 1; 1; part++)
    { copn = fopen(Catenate(name,Numbered_Suffix(".T",part,""),"",""),"r");
      if (copn == NULL)
        break;
      fread(&kmer,sizeof(int),1,copn);
      fread(&n,sizeof(int64),1,copn);
      nels += n;
      fclose(copn);
    }
  if (part == 1)
    { fprintf(stderr,"%s: Cannot find table files for %s\n",Prog_Name,name);
      exit (1);
    }

  //  Allocate in-memory table

  kbyte = (kmer+3)>>2;
  tbyte = kbyte+2;
  table = Malloc(1000l*tbyte,"Allocating k-mer table\n");
  if (table == NULL)
    exit (1);

  //  Load the parts into memory

  fprintf(stderr,"Loading %d-mer table with ",kmer);
  Print_Number(nels,0,stderr);
  fprintf(stderr," entries in %d parts\n",part-1);
  fflush(stderr);

  copn = fopen(Catenate(name,".T1","",""),"r");
  fread(&kmer,sizeof(int),1,copn);
  fread(&n,sizeof(int64),1,copn);
  nels = fread(KMER(0),1000*tbyte,1,copn) / tbyte;

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
  T->name  = Strdup(name,"Allocating name");
  if (T == NULL || T->name == NULL)
    exit (1);

  T->copn  = copn;
  T->part  = 1;
  T->kmer  = kmer;
  T->tbyte = tbyte;
  T->kbyte = kbyte;
  T->nels  = nels;
  T->table = table;

  return (T);
}

int More_Kmer_Table(Kmer_Table *T)
{ int    tbyte = T->tbyte;
  FILE  *copn  = T->copn;
  int    nels, n;
 
  nels = fread(T->table,1000*tbyte,1,copn) / tbyte;
  while (nels == 0)
    { T->part += 1;
      copn = fopen(Catenate(T->name,Numbered_Suffix(".T",T->part,""),"",""),"r");
      if (copn == NULL)
        break;
      fread(&n,sizeof(int),1,copn);
      fread(&n,sizeof(int64),1,copn);
      nels = fread(T->table,1000*tbyte,1,copn) / tbyte;
    }
  if (copn == NULL)
    return (0);
  T->nels = nels;
  T->copn = copn;
  return (nels);
}

void Free_Kmer_Table(Kmer_Table *T)
{ free(T->name);
  free(T->table);
  free(T);
}


/****************************************************************************************
 *
 *  Find Venn Histograms
 *
 *****************************************************************************************/

static int HIST_LOW, HIST_HGH;

void Venn2(Kmer_Table **Tv, int64 **comb)
{ int    tbyte = Tv[0]->tbyte;
  int    kbyte = Tv[0]->kbyte;

  Kmer_Table *T, *U;
  int64      *Inter, *AminB, *BminA;

  int    i, j;
  uint8 *iptr, *jptr;
  int64 *h;
  int    c, d, v;

  T = Tv[0];
  U = Tv[1];

  AminB = comb[0];
  BminA = comb[1];
  Inter = comb[2];

  setup_fmer_table();

  i = 0;
  j = 0;
  iptr = T->table;
  jptr = U->table;
  while (1)
    { if (i >= T->nels)
        { if (More_Kmer_Table(T) == 0)
            { i = j;
              T = U;
              iptr = jptr;
              AminB = BminA;
              break;
            }
          i = 0;
          iptr = T->table;
        }
      if (j >= U->nels)
        { if (More_Kmer_Table(U) == 0)
            break;
          j = 0;
          jptr = U->table;
        }
      v = mycmp(iptr,jptr,kbyte);
      if (v == 0)
        { h = Inter;
          c = COUNT_OF(iptr);
          d = COUNT_OF(jptr);
          if (c > d)
            c = d;
          iptr += tbyte;
          jptr += tbyte;
        }
      else if (v < 0)
        { h = AminB;
          c = COUNT_OF(iptr);
          iptr += tbyte;
        }
      else
        { h = BminA;
          c = COUNT_OF(jptr);
          jptr += tbyte;
        }
      if (c <= HIST_LOW)
        h[HIST_LOW] += 1;
      else if (c >= HIST_HGH)
        h[HIST_HGH] += 1;
      else
        h[c] += 1;
    }

  while (1)
    { if (i >= T->nels)
        { if (More_Kmer_Table(T) == 0)
            break;
          i = 0;
          iptr = T->table;
        }
      c = COUNT_OF(iptr);
      iptr += tbyte;
      if (c <= HIST_LOW)
        AminB[HIST_LOW] += 1;
      else if (c >= HIST_HGH)
        AminB[HIST_HGH] += 1;
      else
        AminB[c] += 1;
    }
}

/*

void Venn(Kmer_Table **T, int nway)
{ int    kmer  = T[0]->kmer;
  int    tbyte = T[0]->tbyte;
  int    kbyte = T[0]->kbyte;

  int    khalf;
  uint8  prefs[] = { 0x3f, 0x0f, 0x03, 0x00 };
  int    mask, offs, rem;

  uint8 *iptr, *nptr;

  int    f;
  uint8 *finger[5];
  uint8 *flimit[4];

  int    a, advn[4];
  int    c, good[4];
  int    mc, hc;
  uint8 *mr, *hr;

  setup_fmer_table();

  for (c = 0; c < nway; c++)
    { idx[c] = 0;
      ptr[c] = T[c]->table;
    }
  while (1)
    { for (c = 0; c < nway; c++)
        if (idx[c] >= T[c]->nels)
          { if (More_Kmer_Table(T[c]) == 0)
              idx[c] = impossible
            else
              { idx[c] = 0;
                ptr[c] = T[c]->table;
              }
          }
      v = mycmp(iptr,jptr,kbyte);
      if (v == 0)
        { if (op == U)
            out iptr cnt1+cnt2
          else if (op == X)
            out iptr min(cnt1,cnt2)
          Inter[COUNT_OF(iptr)] += 1;
          iptr += tbyte;
          jptr += tbyte;
        }
      else if (v < 0)
        { if (op == U || op == A-B)
            out iptr cnt
          AminB[COUNT_OF(iptr)] += 1;
          iptr += tbyte;
        }
      else
        { if (op == U || op == B-A)
            out jptr cnt
          BminA[COUNT_OF(jptr)] += 1;
          jptr += tbyte;
        }
    }
}

*/
          

/****************************************************************************************
 *
 *  Main
 *
 *****************************************************************************************/

int main(int argc, char *argv[])
{ int nway;

  (void) print_pack;

  { int    i, j, k;
    int    flags[128];
    char  *eptr, *fptr;

    (void) flags;
    (void) print_seq;

    ARG_INIT("Vennex");

    HIST_LOW    = 1;
    HIST_HGH    = 100;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
            break;
          case 'h':
            HIST_LOW = strtol(argv[i]+2,&eptr,10);
            if (eptr > argv[i]+2)
              { if (HIST_LOW < 1 || HIST_LOW > 0x7fff)
                  { fprintf(stderr,"%s: Histogram count %d is out of range\n",
                                   Prog_Name,HIST_LOW);
                    exit (1);
                  }
                if (*eptr == ':')
                  { HIST_HGH = strtol(eptr+1,&fptr,10);
                    if (fptr > eptr+1 && *fptr == '\0')
                      { if (HIST_LOW > HIST_HGH)
                          { fprintf(stderr,"%s: Histogram range is invalid\n",Prog_Name);
                            exit (1);
                          }
                        break;
                      }
                  }
                else if (*eptr == '\0')
                  { HIST_HGH = HIST_LOW;
                    HIST_LOW = 1;
                    break;
                  }
              }
            fprintf(stderr,"%s: Syntax of -h option invalid -h[<int(1)>:]<int>\n",Prog_Name);
            exit (1);
        }
      else
        argv[j++] = argv[i];
    argc = j;

    nway = argc-1;
    if (nway < 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  { Kmer_Table *T[nway];
    char       *upp[nway];
    char       *low[nway];
    int64     **comb;
    char       *name;
    int         kmer, ncomb;

    { int64 *hist;
      int    i;

      ncomb = (1 << nway) - 1;
      hist  = Malloc(sizeof(int64)*((HIST_HGH-HIST_LOW)+1),"Allocating histograms");
      comb  = Malloc(sizeof(int64 *)*ncomb,"Allocating histograms");

      comb[0] = hist-HIST_LOW;
      for (i = 1; i < ncomb; i++)
        comb[i] = comb[i-1] + (HIST_HGH-HIST_LOW) + 1;

      bzero(hist,sizeof(int64)*((HIST_HGH-HIST_LOW)+1));
    }

    { int c;

      kmer = 0;
      for (c = 0; c < nway; c++)
        { T[c] = Open_Kmer_Table(argv[c+1],1);
          if (kmer == 0)
            kmer = T[c]->kmer;
          else if (T[c]->kmer != kmer)
            { fprintf(stderr,"%s: K-mer tables do not involve the same K\n",Prog_Name);
              exit (1);
            }
        }
    }

    { int   nlen;
      char *p, *n;
      int   c;

      nlen = nway + 10;
      for (c = 0; c < nway; c++)
        { n = argv[c+1];
          p = index(n,'.');
          if (p != NULL)
            *p = '\0';
          upp[c] = Strdup(n,"Allocating upper case name");
          low[c] = Strdup(n,"Allocating lower case name");
          nlen  += strlen(n);
          if (p != NULL)
            *p = '.';
        }
      name = Malloc(nlen,"Allocating name string");
    }

    if (nway == 2)
      Venn2(T,comb);
    // else
      // Venn(T,nway,comb,ncomb);

    { int   i, b, f, c;
      char *a;

      for (i = 1; i <= ncomb; i++)
        { a = name;
          b = 1;
          for (c = 0; c < nway; c++)
            { if (c != 0)
                a = stpcpy(a,"_");
              if (b & i)
                a = stpcpy(a,upp[c]);
              else
                a = stpcpy(a,low[c]);
            }
          sprintf(a,".K%d",kmer);

          f = open(name,O_CREAT|O_TRUNC|O_WRONLY,S_IRWXU);
          write(f,&kmer,sizeof(int));
          write(f,&HIST_LOW,sizeof(int));
          write(f,&HIST_HGH,sizeof(int));
          write(f,comb[i-1]+HIST_LOW,sizeof(int64)*((HIST_HGH-HIST_LOW)+1));
          close(f);
        }
    }

    { int c;

      free(name);
      for (c = 0; c < nway; c++)
        { free(low[c]);
          free(upp[c]);
        }
      for (c = 0; c < nway; c++)
        Free_Kmer_Table(T[c]);
      free(comb[0]);
    }
  }

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
