 /********************************************************************************************
 *
 *  FastK: a rapid disk-based k-mer counter for high-fidelity shotgun data sets.
 *     Uses a novel minimizer-based distribution scheme that permits problems of
 *     arbitrary size, and a two-staged "super-mer then weighted k-mer" sort to acheive
 *     greater speed when error rates are low (1% or less).  Directly produces sequence
 *     profiles.
 *
 *  Author:  Gene Myers
 *  Date  :  October, 2020
 *
 *********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <dirent.h>
#include <sys/stat.h>
#include <math.h>

// flag MOVE determine if Fastmv or Fastcp

#include "gene_core.h"

static char *Usage = "[-in] <source> <dest>";

int main(int argc, char **argv)
{ int   QUERY;
  int   NO_OVERWRITE;
  char *op;

  { int    i, j, k;
    int    flags[128];

#ifdef MOVE
    ARG_INIT("Fastmv")
    op = "mv";
#else
    ARG_INIT("Fastcp")
    op = "cp";
#endif

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("in")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    QUERY        = flags['i'];
    NO_OVERWRITE = flags['n'];

    if (argc != 3)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -i: prompt for each (stub) overwrite.\n");
        exit (1);
      }
  }

  { int   p, len, a, yes, nthreads;
    char *dir, *root;
    char *DIR, *ROOT;
    char *command;
    struct stat B;
    FILE *f;

    len = strlen(argv[1]) + strlen(argv[2]);
    command = Malloc(len+50,"Allocating command buffer");

    dir  = PathTo(argv[1]);
    root = Root(argv[1],"");

    if (stat(argv[2],&B) == 0 && (B.st_mode & S_IFMT) == S_IFDIR)
      { DIR  = Strdup(argv[2],NULL);
        ROOT = Strdup(root,NULL);
      }
    else
      { DIR  = PathTo(argv[2]);
        ROOT = Root(argv[2],"");
      }

    yes = 1;
    if (stat(Catenate(DIR,"/",ROOT,".hist"),&B) == 0)
      { if (NO_OVERWRITE)
          yes = 0;
        else if (QUERY)
          { printf("overwrite %s/%s.hist? ",DIR,ROOT);
            yes = 0;
            while ((a = getc(stdin)) != '\n')
              if (a == 'y' || a == 'Y')
                yes = 1;
              else if (a == 'n' || a == 'N')
                yes = 0;
	  }
      } 
    if (yes)
      { sprintf(command,"%s -f %s/%s.hist %s/%s.hist",op,dir,root,DIR,ROOT);
        system(command);
      }

    yes = 1;
    if (stat(Catenate(DIR,"/",ROOT,".ktab"),&B) == 0)
      { if (NO_OVERWRITE)
          yes = 0;
        else if (QUERY)
          { printf("overwrite %s/%s.ktab? ",DIR,ROOT);
            yes = 0;
            while ((a = getc(stdin)) != '\n')
              if (a == 'y' || a == 'Y')
                yes = 1;
              else if (a == 'n' || a == 'N')
                yes = 0;
          }
      }
    if (yes)
      { f = fopen(Catenate(dir,"/",root,".ktab"),"r");
        if (f != NULL)
          { fread(&nthreads,sizeof(int),1,f);
            fread(&nthreads,sizeof(int),1,f);
            fclose(f);
            for (p = 1; p <= nthreads; p++)
              { sprintf(command,"%s -f %s/.%s.ktab.%d %s/.%s.ktab.%d",op,dir,root,p,DIR,ROOT,p);
                system(command);
              }
            sprintf(command,"%s -f %s/%s.ktab %s/%s.ktab",op,dir,root,DIR,ROOT);
            system(command);
          }
      }

    yes = 1;
    if (stat(Catenate(DIR,"/",ROOT,".prof"),&B) == 0)
      { if (NO_OVERWRITE)
          yes = 0;
        else if (QUERY)
          { printf("overwrite %s/%s.prof? ",DIR,ROOT);
            yes = 0;
            while ((a = getc(stdin)) != '\n')
              if (a == 'y' || a == 'Y')
                yes = 1;
              else if (a == 'n' || a == 'N')
                yes = 0;
          }
      }
    if (yes)
      { f = fopen(Catenate(dir,"/",root,".prof"),"r");
        if (f != NULL)
          { fread(&nthreads,sizeof(int),1,f);
            fread(&nthreads,sizeof(int),1,f);
            fclose(f);
            for (p = 1; p <= nthreads; p++)
              { sprintf(command,"%s -f %s/.%s.pidx.%d %s/.%s.pidx.%d",op,dir,root,p,DIR,ROOT,p);
                system(command);
                sprintf(command,"%s -f %s/.%s.prof.%d %s/.%s.prof.%d",op,dir,root,p,DIR,ROOT,p);
                system(command);
              }
            sprintf(command,"%s -f %s/%s.prof %s/%s.prof",op,dir,root,DIR,ROOT);
            system(command);
          }
      }

    free(root);
    free(dir);
    free(ROOT);
    free(DIR);

    free(command);
  }

  exit (0);
}
