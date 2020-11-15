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
        List_Kmer_Table(T);
      else if (strcmp(argv[c],"CHECK") == 0)
        Check_Kmer_Table(T);
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
