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

  { int   p, len, a, yes, any, nthreads;
    int   alen, rlen;
    int   doh, dot, dop;
    char *dir, *root;
    char *DIR, *ROOT;
    char *command;
    struct stat B;
    FILE *f;

    len = strlen(argv[1]) + strlen(argv[2]);
    command = Malloc(len+50,"Allocating command buffer");

    dir  = PathTo(argv[1]);
    root = Root(argv[1],NULL);

    alen = strlen(argv[1]);
    rlen = strlen(root);
    any  = 0;
    if (alen == rlen)
      doh = dot = dop = 1;
    else
      { doh = (strcmp(argv[1] + rlen,".hist") == 0);
        dot = (strcmp(argv[1] + rlen,".ktab") == 0);
        dop = (strcmp(argv[1] + rlen,".prof") == 0);
        if (doh + dot + dop == 0)
          { free(root);
            root = Strdup(argv[1],NULL);
            doh = dot = dop = 1;
          }
        else
          any = 1;
      }

    if (stat(argv[2],&B) == 0 && (B.st_mode & S_IFMT) == S_IFDIR)
      { DIR  = Strdup(argv[2],NULL);
        ROOT = Strdup(root,NULL);
      }
    else
      { DIR  = PathTo(argv[2]);
        ROOT = Root(argv[2],"");
      }

    yes = doh;
    if (doh && stat(Catenate(DIR,"/",ROOT,".hist"),&B) == 0)
      { any = 1;
        if (NO_OVERWRITE)
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

    yes = dot;
    if (dot && stat(Catenate(DIR,"/",ROOT,".ktab"),&B) == 0)
      { any = 1;
        if (NO_OVERWRITE)
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

    yes = dop;
    if (dop && stat(Catenate(DIR,"/",ROOT,".prof"),&B) == 0)
      { any = 1;
        if (NO_OVERWRITE)
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

    if (any == 0)
      fprintf(stderr,"%s: Warning, no FastK output files with root %s\n",Prog_Name,root);

    free(root);
    free(dir);
    free(ROOT);
    free(DIR);

    free(command);
  }

  exit (0);
}
