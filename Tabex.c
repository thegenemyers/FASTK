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

static char *Usage = "[-t<int>] <source_root>[.ktab] (LIST|CHECK|(k-mer:string>) ...";

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

static void List_Kmer_Table(Kmer_Table *T, FILE *out)
{ char *seq;
  int64 i;

  seq = Fetch_Kmer(T,0,NULL);
  for (i = 0; i < T->nels; i++)
    fprintf(out," %9lld: %s = %5d\n",i,Fetch_Kmer(T,i,seq),Fetch_Count(T,i));
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

  seq = Current_Kmer(S,NULL);
  for (First_Kmer_Entry(S); S->csuf != NULL; Next_Kmer_Entry(S))
    if (Current_Count(S) >= cut)
      { fprintf(out," %9lld: %s = %5d\n",S->cidx,Current_Kmer(S,seq),Current_Count(S)); fflush(stdout); }
  free(seq);
}

int main(int argc, char *argv[])
{ Kmer_Table  *T;
  Kmer_Stream *S;
  int          CUT;
  int          STREAM;

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
            ARG_FLAGS("T")
            break;
          case 't':
            ARG_POSITIVE(CUT,"Cutoff for k-mer table")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    STREAM = ! flags['T'];   //  This is undocumented and only for developer use.

    if (argc < 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -t: Trim all k-mers with counts less than threshold\n");
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
    
      printf("Opening %d-mer table with ",S->kmer);
      Print_Number(S->nels,0,stdout);
      printf(" entries");
      if (S->minval > 1)
        printf(" occuring %d-or-more times",S->minval);
      printf("\n");
      fflush(stdout);
    
      { int   c;
        char *seq;
    
        seq = NULL;
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
        free(seq);
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
            List_Kmer_Table(T,stdout);
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

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
