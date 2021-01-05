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

/****************************************************************************************
 *
 *  Test Stub
 *
 *****************************************************************************************/

static int Check_Kmer_Table(Kmer_Table *T)
{ char *curs, *last;
  int64 i;

  curs = NULL;
  last = Fetch_Kmer(T,0,NULL);
  for (i = 1; i < T->nels; i++)
    { curs = Fetch_Kmer(T,i,curs);
      if (strcmp(last,curs) >= 0)
        { fprintf(stderr,"\nOut of Order\n");
          fprintf(stderr," %9lld: %s = %5d\n",i-1,last,Fetch_Count(T,i-1));
          fprintf(stderr," %9lld: %s = %5d\n",i,curs,Fetch_Count(T,i));
          break;
        }
    }
  free(curs);
  free(last);
  return (i >= T->nels);
}

void List_Kmer_Table(Kmer_Table *T, FILE *out)
{ char *seq;
  int64 i;

  seq = NULL;
  for (i = 0; i < T->nels; i++)
    { seq = Fetch_Kmer(T,i,seq);
      fprintf(out," %9lld: %s = %5d\n",i,seq,Fetch_Count(T,i));
    }
}

int main(int argc, char *argv[])
{ Kmer_Table *T;
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

  T = Load_Kmer_Table(argv[1],CUT);
  if (T == NULL)
    { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,argv[1]);
      exit (1);
    } 

  fprintf(stderr,"Loaded %d-mer table with ",T->kmer);
  Print_Number(T->nels,0,stderr);
  fprintf(stderr," entries\n");
  fflush(stderr);

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
                printf("%*s: %5d\n",T->kmer,argv[c],Fetch_Count(T,loc));
            }
        }
  }

  Free_Kmer_Table(T);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
