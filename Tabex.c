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

#include "libfastk.h"
#include "ONElib.h"

static char *Usage = "[-1AC] [-t<int>] <source>[.ktab] [ <address>[-<address>] ]";

static char *One_Schema =
  "1 3 def 2 1                   schema for k-mer table\n"
  ".\n"
  "P 3 kmr                       This is a k-mer table 1-code file\n"
  "D K 4 3 INT 3 INT 3 INT 3 INT k-mer size, prefix length, min. count, & 1st prefix for table\n"
  "O S 1 3 DNA                   concatentation of the suffixes of the k-mers with given prefix\n"
  "D C 1 8 INT_LIST              counts of the suffixes with the given prefix (in lex order)\n"
;

static void Check_Kmer_Stream(Kmer_Stream *S, int64 bidx, int64 eidx)
{ int    hbyte = S->hbyte;
  int    lpre, u;
  uint8 *lsuf;
  char  *seq;

  lsuf = (uint8 *) (seq = Current_Kmer(S,NULL));
  GoTo_Kmer_Index(S,bidx);
  if (S->csuf == NULL)
    return;
  lpre = S->cpre;
  memcpy(lsuf,S->csuf,hbyte);
  for (Next_Kmer_Entry(S); S->cidx < eidx; Next_Kmer_Entry(S))
    { if (S->cpre < lpre || (S->cpre == lpre && memcmp(S->csuf,lsuf,hbyte) <= 0))
        { printf("\nOut of Order %02x %02x\n",S->cpre,lpre);
          for (u = 0; u < S->hbyte; u++)
            printf(" %02x %02x\n",S->csuf[u],lsuf[u]);
          printf(" %9lld: %s = %5d\n",S->cidx,Current_Kmer(S,seq),Current_Count(S));
          GoTo_Kmer_Index(S,S->cidx-1);
          printf(" %9lld: %s = %5d\n",S->cidx,Current_Kmer(S,seq),Current_Count(S));
          Next_Kmer_Entry(S);
          printf(" %9lld: %s = %5d\n",S->cidx,Current_Kmer(S,seq),Current_Count(S));
          break;
        }
      lpre = S->cpre;
      memcpy(lsuf,S->csuf,hbyte);
    }
  free(seq);
  if (S->cidx >= eidx)
    printf("\nTable is OK\n");
  return;
}

static void List_Kmer_Stream(Kmer_Stream *S, int cut, int ASCII, int64 bidx, int64 eidx)
{ char *seq;
  int   c;

  seq = Current_Kmer(S,NULL);
  GoTo_Kmer_Index(S,bidx);
  if (ASCII)
    for (GoTo_Kmer_Index(S,bidx); S->cidx < eidx; Next_Kmer_Entry(S))
      { c = Current_Count(S);
        if (c >= cut)
          printf("%s\t%d\n",Current_Kmer(S,seq),c);
      }
  else
    for (GoTo_Kmer_Index(S,bidx); S->cidx < eidx; Next_Kmer_Entry(S))
      { c = Current_Count(S);
        if (c >= cut)
          printf(" %9lld: %s = %5d\n",S->cidx,Current_Kmer(S,seq),c);
      }
  free(seq);
}

static void One_Kmer_Stream(Kmer_Stream *S, int cut, int64 bidx, int64 eidx, char *command)
{ OneSchema *schema;
  OneFile   *file1;
  char      *suffix, *seq;
  int64     *counts;
  int        gmer, smer, imer, maxlist;

  int        cp, sp;
  int        i, k, c;
  char      *s;

  extern int Max_Bucket(Kmer_Stream *);

  schema = oneSchemaCreateFromText(One_Schema);
  file1  = oneFileOpenWriteNew("-",schema,"kmr",true,1);
  oneAddProvenance(file1,Prog_Name,"1.0","%s >?.kmr",command);

  imer = S->ibyte;
  gmer = 4*imer;
  smer = S->kmer - gmer;

  GoTo_Kmer_Index(S,bidx);

  oneInt(file1,0) = S->kmer;
  oneInt(file1,1) = gmer;
  oneInt(file1,2) = S->minval;
  oneInt(file1,3) = S->cpre;
  oneWriteLine(file1,'K',0,NULL);

  seq     = Current_Kmer(S,NULL);
  maxlist = Max_Bucket(S);
  suffix  = Malloc(maxlist*smer,"Allocating suffix array");
  counts  = Malloc(maxlist*sizeof(int64),"Allocating count array");

  cp = sp = 0;
  for (i = S->cpre; S->cidx < eidx; Next_Kmer_Entry(S))
    { c = Current_Count(S);
      if (S->cpre != i)
        { oneWriteLine(file1,'S',sp,suffix);
          oneWriteLine(file1,'C',cp,counts);
          i = S->cpre;
          cp = sp = 0;
        }
      if (c >= cut)
        { counts[cp++] = c;
          s = Current_Kmer(S,seq)+gmer;
          for (k = 0; k < smer; k++)
            suffix[sp++] = *s++;
        }
    }
  oneWriteLine(file1,'S',sp,suffix);
  oneWriteLine(file1,'C',cp,counts);
  
  free(counts);
  free(suffix);
  free(seq);

  oneFileClose(file1);
  oneSchemaDestroy(schema);
}

static int dna[128] =
  { 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,

    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
  
    0, 1, 0, 1, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    
    0, 1, 0, 1, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
  };  

static int shiftup[128] =
  { 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
          
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
        
    0, 'C', 0, 'G', 0, 0, 0, 'T',
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
          
    0, 'c', 0, 'g', 0, 0, 0, 't',
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
  };
    
static uint8 code[128] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  
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
  s0[len] = s1[len] = s2[len] = 'a';
  
  for (i = 0; i < len; i += 4)
    *t++ = ((code[(int) s0[i]] << 6) | (code[(int) s1[i]] << 4)
         |  (code[(int) s2[i]] << 2) | code[(int) s3[i]] );
          
  s0[len] = c;
  s1[len] = d;
  s2[len] = e;
}       

static int64 Interpret(Kmer_Stream *T, char *x, int beg)
{ int    d, n;
  char  *u;
  uint8 *v;

  if (sscanf(x,"%d%n",&d,&n) == 1)
    { if (x[n] != 0)
        { fprintf(stderr,"%s: Indx %s is not an integer\n",Prog_Name,x);
          exit (1);
        }
      if (d >= T->nels)
        { fprintf(stderr,"%s: Index %s is out of bounds\n",Prog_Name,x);
          exit (1);
        }
      if (beg)
        return ((int64) d);
      else
        return ((int64) d+1);
    }
  for (n = 0; x[n] != '\0'; n++)
    if (!dna[(int) x[n]])
      { fprintf(stderr,"%s: String %s is not dna (acgt)\n",Prog_Name,x);
        exit (1);
      }
  if (n > T->kmer)
    { fprintf(stderr,"%s: String %s is longer than k-mer size (%d)\n",Prog_Name,x,T->kmer);
      exit (1);
    }
  u = Current_Kmer(T,NULL);
  v = (uint8 *) Current_Kmer(T,NULL);
  strcpy(u,x);
  if (!beg)
    { n -= 1;
      while (n >= 0 && (u[n] == 't' || u[n] == 'T'))
        n -= 1;
      if (n < 0)
        return (T->nels);
      else
        u[n] = shiftup[(int) u[n]];
      n += 1;
    }
  while (n < T->kmer)
    u[n++] = 'a';
  compress_norm(u,T->kmer,v);
  GoTo_Kmer_Entry(T,v);
  free(v);
  free(u);
  return (T->cidx);
}

int main(int argc, char *argv[])
{ Kmer_Stream *S;
  char        *command;
  int          CUT;
  int          ASCII;
  int          CHECK;
  int          ONE_CODE;
  int64        bidx, eidx;

  //  Process options and capture command line for provenance

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    (void) flags;

    ARG_INIT("Tabex");

    { int   n, t;
      char *c;

      n = 0;
      for (t = 0; t < argc; t++)
        n += strlen(argv[t])+1;

      command = Malloc(n+1,"Allocating command string");
      if (command == NULL)
        exit (1);

      c = command;
      if (argc >= 1)
        { c += sprintf(c,"%s",argv[0]);
          for (t = 1; t < argc; t++)
            c += sprintf(c," %s",argv[t]);
        }
      *c = '\0';
    }

    CUT = 0;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("1AC")
            break;
          case 't':
            ARG_POSITIVE(CUT,"Cutoff for k-mer table")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    ONE_CODE = flags['1'];
    ASCII    = flags['A'];
    CHECK    = flags['C'];

    if (argc < 2 || argc > 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"          <address> = <int> | <dna:string>\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -t: Trim all k-mers with counts less than threshold\n");
        fprintf(stderr,"      -A: Output tab-delimited ASCII\n");
        fprintf(stderr,"      -C: Check sorting\n");
        fprintf(stderr,"      -1: Produce 1-code as output.\n");
        exit (1);
      }

    if (CHECK && (ONE_CODE || ASCII || CUT > 0))
      { fprintf(stderr,"%s: -C option incompatible with all other options\n",Prog_Name);
        exit (1);
      }
  }

  S = Open_Kmer_Stream(argv[1]);
  if (S == NULL)
    { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,argv[1]);
      exit (1);
    } 

  if (argc == 2)
    { bidx = 0;
      eidx = S->nels;
    }   
  else    
    { char *x = argv[2];
      char *p = index(x,'-');
      if (p != NULL)
        { *p++ = '\0';
          bidx = Interpret(S,x,1);
          eidx = Interpret(S,p,0);
        }
      else
        { bidx = Interpret(S,x,1);
          eidx = Interpret(S,x,0);
        }
    } 

  if (bidx == eidx)
    printf("\nNothing found in range given !\n");

  else if (CHECK)
    Check_Kmer_Stream(S,bidx,eidx);

  else if (ONE_CODE)
    One_Kmer_Stream(S,CUT,bidx,eidx,command);
 
  else
    { if (!ASCII)
        { printf("Opening %d-mer table with ",S->kmer);
          Print_Number(S->nels,0,stdout);
          printf(" entries");
          if (S->minval > 1)
            printf(" occuring %d-or-more times",S->minval);
          printf("\n");
          fflush(stdout);
        }

      List_Kmer_Stream(S,CUT,ASCII,bidx,eidx);
    }

  Free_Kmer_Stream(S);

  free(command);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
