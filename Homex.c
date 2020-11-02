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

#undef  DEBUG_PARTITION
#undef  DEBUG_QUEUE

#include "gene_core.h"

static char *Usage = "-e<int> -g<int>:<int> <source_root>.K<k>";

#define MAX_HOMO_LEN 14

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

  b = (len >> 2);
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
  

Profile *Count_Homopolymer_Errors(Kmer_Table *T)
{ static Profile profile;

  int    kmer  = T->kmer;
  int    tbyte = T->tbyte;
  int    kbyte = T->kbyte;
  int64  nels  = T->nels;
  uint8 *table = T->table;

  Point *counter;
  int    khalf, klong, kbase, kchkl, kextn;

  uint8  suffs[]  = { 0x00, 0xc0, 0xf0, 0xfc };
  uint8  abyte[]  = { 0x00, 0x55, 0xaa, 0xff };

  uint8 *iptr, *nptr;

  int    i;
  uint8 *fing[4];
  uint8 *fend[4];
  uint8 *fbeg[5];

  int    hlen, hsym;
  uint8  suffix[kbyte];

  int   a, b, advn[4];
  int   cn[4];

  fbeg[4] = table;
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

  nptr    = KMER(nels);
  iptr    = table;
  while (iptr < nptr)
    { hlen = khalf-1;
      hsym = SYMBOL(iptr,hlen);
      for (hlen--; hlen >= klong; hlen--)
        if (SYMBOL(iptr,hlen) != hsym)
          break;
      hlen += 1;

      if (hlen < klong)
        { for (iptr += tbyte; iptr < nptr; iptr += tbyte)
            { int x = mypref(iptr,iptr-tbyte,khalf); 
              if (x < khalf)
                break;
            }
          continue;
        }

      hlen = khalf-hlen;
#ifdef DEBUG_PARTITION
      print_seq(iptr,kmer);
      printf("  Len = %d  Sym = %c\n",hlen,dna[hsym]);
      fflush(stdout);
#endif

      { int k;

        mycpy(suffix,iptr,kbyte);
        k = (khalf>>2);
        suffix[k] = (suffix[k] & suffs[khalf&0x3]) | (abyte[hsym] & ~suffs[khalf&0x3]);
        for (k++; k < kbyte; k++)
          suffix[k] = abyte[hsym];
      }
#ifdef DEBUG_PARTITION
      print_seq(suffix,kmer);
      printf("\n");
      fflush(stdout);
#endif

      kbase = khalf + (hlen-1);
      kchkl = khalf + (hlen+2);
      kextn = kmer - kbase;
      for (i = 0; i <= 3; i++)
        fend[i] = NULL;
    
      for (; iptr < nptr; iptr += tbyte)
        { int x = mypref(iptr,suffix,kchkl); 
          if (x < khalf)
            break;
          x -= kbase;
          if (0 <= x && x <= 3)
            { if (fend[x] == NULL)
                fbeg[x] = iptr;
              fend[x] = iptr;
            }
        }

#ifdef DEBUG_PARTITION
      for (i = 0; i <= 3; i++)
        if (fend[i] != NULL)
          printf(" %ld-%ld",(fbeg[i]-table)/tbyte,(fend[i]-table)/tbyte);
      printf(" >> %ld\n",(iptr-table)/tbyte);
      fflush(stdout);
#endif

      if (fend[1] == NULL && fend[2] == NULL)
        continue;

      for (i = 3; i >= 0; i--)
        if (fend[i] == NULL)
          fing[i] = fend[i] = fbeg[i] = table;
        else
          { fing[i] = fbeg[i];
            fend[i] += tbyte;
          }

#ifdef DEBUG_PARTITION
      for (i = 0; i <= 3; i++)
        if (fend[i] == table)
          printf(" ***");
        else
          printf(" %ld-%ld",(fbeg[i]-table)/tbyte,(fend[i]-table)/tbyte);
      printf(" >> %ld\n",(iptr-table)/tbyte);
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
              { int v = mybpcmp(fing[b],fing[i],kbase+b,kbase+i,kextn-i);
                if (v == 0)
                  ADD(i)
                else if (v < 0)
                  SET(i)
              }

#ifdef DEBUG_QUEUE
          for (i = 0; i < a; i++)
            printf(" %d(%d)",advn[i],COUNT_OF(fing[advn[i]]));
#endif

          cn[0] = cn[1] = cn[2] = cn[3] = 0;
          for (i = 0; i < a; i++)
            { b = advn[i];
              cn[b] = COUNT_OF(fing[b]);
              fing[b] += tbyte;
              if (fing[b] == fbeg[b+1])
                fing[b] = fend[b+1];
            }

          if (GOOD_LOW <= cn[1] && cn[1] <= GOOD_HGH && cn[0] <= ERROR && cn[2] <= ERROR)
            { counter[hlen].correct += cn[1];
              counter[hlen].lessone += cn[0]; 
              counter[hlen].plusone += cn[2]; 
#ifdef DEBUG_QUEUE
              printf(" -> %d%c %d:%d:%d\n",hlen,dna[hsym],cn[0],cn[1],cn[2]);
#endif
            }
          else if (GOOD_LOW <= cn[2] && cn[2] <= GOOD_HGH && cn[1] <= ERROR && cn[3] <= ERROR)
            { if (hlen < MAX_HOMO_LEN)
                { counter[hlen+1].correct += cn[2];
                  counter[hlen+1].lessone += cn[1]; 
                  counter[hlen+1].plusone += cn[3]; 
#ifdef DEBUG_QUEUE
                  printf(" -> %d%c %d:%d:%d\n",hlen+1,dna[hsym],cn[1],cn[2],cn[3]);
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
{ Kmer_Table *T;
  Profile    *P;

  (void) print_pack;
  (void) print_seq;
  (void) setup_fmer_table;
  (void) mycmp;
  (void) mycpy;

  { int    i, j, k;
    int    flags[128];
    char  *eptr, *fptr;

    (void) flags;

    ARG_INIT("Haplex");

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

    if (ERROR < 0)
      { fprintf(stderr,"%s: Must give error threshold -e\n",Prog_Name);
        exit (1);
      }
    if (GOOD_LOW < 0)
      { fprintf(stderr,"%s: Must give good count range -g\n",Prog_Name);
        exit (1);
      }
    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -e: Counts <= this value are considered errors.\n");
        fprintf(stderr,"      -g: Counts in this range are considered correct.\n");
        exit (1);
      }
  }

  T = Load_Kmer_Table(argv[1],1);

  P = Count_Homopolymer_Errors(T);

  Free_Kmer_Table(T);

  { int h;

    for (h = 2; h <= MAX_HOMO_LEN; h++)
      { int64 cc = (*P)[0][h].correct + (*P)[3][h].correct;
        int64 cl = (*P)[0][h].lessone + (*P)[3][h].lessone;
        int64 cp = (*P)[0][h].plusone + (*P)[3][h].plusone;

        printf(" %2d at: %10lld %10lld %10lld -> %.1f%%\n",h,cl,cc,cp,(100.*(cl+cp))/(cc+cl+cp));

        cc = (*P)[1][h].correct + (*P)[2][h].correct;
        cl = (*P)[1][h].lessone + (*P)[2][h].lessone;
        cp = (*P)[1][h].plusone + (*P)[2][h].plusone;

        printf(" %2d cg: %10lld %10lld %10lld -> %.1f%%\n",h,cl,cc,cp,(100.*(cl+cp))/(cc+cl+cp));
      }
  }

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
