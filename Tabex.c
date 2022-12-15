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

static char *Usage = "[-t<int>] <source_root>[.ktab] (-1 | (LIST|CHECK|(k-mer:string>) ...)";

static char *One_Schema =
  "P 3 kmr                       This is a k-mer table 1-code file\n"
  "D K 3 3 INT 3 INT 3 INT       k-mer size, prefix length, and entry min. count for table\n"
  "O S 1 3 DNA                   concatentation of the suffixes of the k-mers with given prefix\n"
  "D C 1 8 INT_LIST              counts of the suffixes with the given prefix (in lex order)\n";

static int Check_Kmer_Table(Kmer_Table *T)
{ char *curs, *last, *flip;
  int64 i;

  curs = Fetch_Kmer(T,0,NULL);
  last = Fetch_Kmer(T,0,NULL);
  for (i = 1; i < T->nels; i++)
    { if (strcmp(last,Fetch_Kmer(T,i,curs)) >= 0)
        { fprintf(stderr,"\nOut of Order\n");
          fprintf(stderr," %9lld: %s = %5d\n",i-1,last,Fetch_Count(T,i-1));
          fprintf(stderr," %9lld: %s = %5d\n",i,curs,Fetch_Count(T,i));
          fflush(stderr);
          break;
        }
      flip = curs;
      curs = last;
      last = flip;
    }
  free(curs);
  free(last);
  return (i >= T->nels);
}

static void List_Kmer_Table(Kmer_Table *T, int cut, FILE *out)
{ char *seq;
  int64 i;
  int   c;

  seq = Fetch_Kmer(T,0,NULL);
  for (i = 0; i < T->nels; i++)
    { c = Fetch_Count(T,i);
      if (c >= cut)
        fprintf(out," %9lld: %s = %5d\n",i,Fetch_Kmer(T,i,seq),c);
    }
  free(seq);
}

static int Check_Kmer_Stream(Kmer_Stream *S)
{ int    hbyte = S->hbyte;
  int    lpre, u;
  uint8 *lsuf;
  char  *seq;

  lsuf = (uint8 *) (seq = Current_Kmer(S,NULL));
  First_Kmer_Entry(S);
  if (S->csuf == NULL)
    return (1);
  lpre = S->cpre;
  memcpy(lsuf,S->csuf,hbyte);
  for (Next_Kmer_Entry(S); S->csuf != NULL; Next_Kmer_Entry(S))
    { if (S->cpre < lpre || (S->cpre == lpre && memcmp(S->csuf,lsuf,hbyte) < 0))
        { fprintf(stderr,"\nOut of Order %02x %02x\n",S->cpre,lpre);
          for (u = 0; u < S->hbyte; u++)
            printf(" %02x %02x\n",S->csuf[u],lsuf[u]);
          fprintf(stderr," %9lld: %s = %5d\n",S->cidx,Current_Kmer(S,seq),Current_Count(S));
          GoTo_Kmer_Index(S,S->cidx-1);
          fprintf(stderr," %9lld: %s = %5d\n",S->cidx,Current_Kmer(S,seq),Current_Count(S));
          Next_Kmer_Entry(S);
          fprintf(stderr," %9lld: %s = %5d\n",S->cidx,Current_Kmer(S,seq),Current_Count(S));
          fflush(stderr);
          break;
        }
      lpre = S->cpre;
      memcpy(lsuf,S->csuf,hbyte);
    }
  free(seq);
  return (S->csuf == NULL);
}

static void List_Kmer_Stream(Kmer_Stream *S, int cut, FILE *out)
{ char *seq;
  int   c;

  seq = Current_Kmer(S,NULL);
  for (First_Kmer_Entry(S); S->csuf != NULL; Next_Kmer_Entry(S))
    { c = Current_Count(S);
      if (c >= cut)
        fprintf(out," %9lld: %s = %5d\n",S->cidx,Current_Kmer(S,seq),c);
    }
  free(seq);
}

static void One_Kmer_Stream(Kmer_Stream *S, int cut, char *command)
{ OneSchema *schema;
  OneFile   *file1;
  char      *suffix, *seq;
  int64     *counts;
  int        gmer, smer, imer, nbuck, maxlist;

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

  oneInt(file1,0) = S->kmer;
  oneInt(file1,1) = gmer;
  oneInt(file1,2) = S->minval;
  oneWriteLine(file1,'K',0,NULL);

  seq     = Current_Kmer(S,NULL);
  maxlist = Max_Bucket(S);
  suffix  = Malloc(maxlist*smer,"Allocating suffix array");
  counts  = Malloc(maxlist*sizeof(int64),"Allocating count array");

  nbuck = (1 << (2*gmer));
  First_Kmer_Entry(S);
  for (i = 0; i < nbuck; i++)
    { cp = sp = 0;
      while (S->cpre == i)
        { c = Current_Count(S);
          if (c >= cut)
            { counts[cp++] = c;
              s = Current_Kmer(S,seq)+gmer;
              for (k = 0; k < smer; k++)
                suffix[sp++] = *s++;
            }
          Next_Kmer_Entry(S);
        }
      oneWriteLine(file1,'S',sp,suffix);
      oneWriteLine(file1,'C',cp,counts);
    }

  free(counts);
  free(suffix);
  free(seq);

  oneFileClose(file1);
  oneSchemaDestroy(schema);
}

int main(int argc, char *argv[])
{ Kmer_Table  *T;
  Kmer_Stream *S;
  char        *command;
  int          CUT;
  int          STREAM;
  int          ONE_CODE;

  //  Process options and capture command line for provenance

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    (void) flags;

    ARG_INIT("Tabex");

    { int   n, t;
      char *c;

      n = 0;
      for (t = 1; t < argc; t++)
        n += strlen(argv[t])+1;

      command = Malloc(n+1,"Allocating command string");
      if (command == NULL)
        exit (1);

      c = command;
      if (argc >= 1)
        { c += sprintf(c,"%s",argv[1]);
          for (t = 2; t < argc; t++)
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
            ARG_FLAGS("T1")
            break;
          case 't':
            ARG_POSITIVE(CUT,"Cutoff for k-mer table")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    STREAM   = ! flags['T'];   //  This is undocumented and only for developer use.
    ONE_CODE = flags['1'];

    if ((ONE_CODE && argc != 2) || (!ONE_CODE && argc < 3))
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -t: Trim all k-mers with counts less than threshold\n");
        fprintf(stderr,"      -1: Produce 1-code as output.\n");
        exit (1);
      }
  }

  //  Default is to use a Kmer_Stream to realize operations

  if (STREAM)
    { S = Open_Kmer_Stream(argv[1]);
      if (S == NULL)
        { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,argv[1]);
          exit (1);
        } 
    
      if (ONE_CODE)
        One_Kmer_Stream(S,CUT,command);
     
      else
        { int   c;

          printf("Opening %d-mer table with ",S->kmer);
          Print_Number(S->nels,0,stdout);
          printf(" entries");
          if (S->minval > 1)
            printf(" occuring %d-or-more times",S->minval);
          printf("\n");
          fflush(stdout);
    
          for (c = 2; c < argc; c++)
            if (strcmp(argv[c],"LIST") == 0)
              List_Kmer_Stream(S,CUT,stdout);
            else if (strcmp(argv[c],"CHECK") == 0)
              { if (Check_Kmer_Stream(S))
                  printf("The table is OK\n");
              }
            else
              { if ((int) strlen(argv[c]) != S->kmer)
                  printf("%*s: Not a %d-mer\n",S->kmer,argv[c],S->kmer);
                else
                  { if (GoTo_Kmer_String(S,argv[c]))
                      printf("%*s: %5d @ idx = %lld\n",S->kmer,argv[c],Current_Count(S),S->cidx);
                    else
                      printf("%*s: Not found\n",S->kmer,argv[c]);
                  }
              }
        }
    
      Free_Kmer_Stream(S);
    }

  //  But for developers and illustrative purposes we also give a Kmer_Table implementation

  else
    { T = Load_Kmer_Table(argv[1],CUT);
      if (T == NULL)
        { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,argv[1]);
          exit (1);
        } 
    
      printf("Loaded %d-mer table with ",T->kmer);
      Print_Number(T->nels,0,stdout);
      printf(" entries");
      if (T->minval > 1)
        printf(" occuring %d-or-more times",T->minval);
      printf("\n");
      fflush(stdout);
    
      { int   c;
        int64 loc;
    
        for (c = 2; c < argc; c++)
          if (strcmp(argv[c],"LIST") == 0)
            List_Kmer_Table(T,CUT,stdout);
          else if (strcmp(argv[c],"CHECK") == 0)
            { if (Check_Kmer_Table(T))
                printf("The table is OK\n");
            }
          else
            { if ((int) strlen(argv[c]) != T->kmer)
                printf("%*s: Not a %d-mer\n",T->kmer,argv[c],T->kmer);
              else
                { loc = Find_Kmer(T,argv[c]);
                  if (loc < 0)
                    printf("%*s: Not found\n",T->kmer,argv[c]);
                  else
                    printf("%*s: %5d @ idx = %lld\n",T->kmer,argv[c],Fetch_Count(T,loc),loc);
                }
            }
      }
    
      Free_Kmer_Table(T);
    }

  free(command);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
