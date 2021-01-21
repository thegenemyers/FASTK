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

static char *Usage = "[-inf] <source> <dest>";

int main(int argc, char **argv)
{ int   QUERY;
  int   QUIET;
  int   NO_OVERWRITE;
  char *op;
  struct stat B;

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
            ARG_FLAGS("inf")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    QUERY        = flags['i'];
    QUIET        = flags['f'];
    NO_OVERWRITE = flags['n'];
    if (QUIET)
      QUERY = NO_OVERWRITE = 0;

    if (argc > 3)
      { if ( ! (stat(argv[argc-1],&B) == 0 && (B.st_mode & S_IFMT) == S_IFDIR))
          argc = 2;
      }
    if (argc < 3)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -i: prompt for each (stub) overwrite.\n");
        fprintf(stderr,"      -n: do not overwrite existing files.\n");
        fprintf(stderr,"      -f: force operation quietly\n");
        exit (1);
      }
  }

  { int   p, len, a, yes, any, nthreads;
    int   alen, rlen;
    int   doh, dot, dop;
    char *dir, *root, *fname;
    char *DIR, *ROOT;
    char *command;
    FILE *f;
    int   c;


    len = 0;
    for (c = 1; c < argc; c++)
      if ((int) strlen(argv[c]) > len)
        len = strlen(argv[c]);

    command = Malloc(2*len+100,"Allocating command buffer");

    for (c = 1; c < argc-1; c++)
      {
        dir  = PathTo(argv[c]);
        fname = rindex(argv[c],'/');
        if (fname == NULL)
          fname = argv[c];
        else
          fname = fname + 1;
        root = Root(fname,NULL);

        alen = strlen(fname);
        rlen = strlen(root);
        if (alen == rlen)
          doh = dot = dop = 1;
        else
          { doh = (strcmp(fname + rlen,".hist") == 0);
            dot = (strcmp(fname + rlen,".ktab") == 0);
            dop = (strcmp(fname + rlen,".prof") == 0);
            if (doh + dot + dop == 0)
              { free(root);
                root = Strdup(fname,NULL);
                doh = dot = dop = 1;
              }
          }

        any = 0;
        if (doh && stat(Catenate(dir,"/",root,".hist"),&B) == 0)
          any = 1;
        if (dot && stat(Catenate(dir,"/",root,".ktab"),&B) == 0)
          any = 1;
        if (dop && stat(Catenate(dir,"/",root,".prof"),&B) == 0)
          any = 1;
        if (any == 0 && !QUIET)
          { if (doh+dot+dop == 3)
              fprintf(stderr,"%s: no FastK output files with root %s\n",Prog_Name,root);
            else
              fprintf(stderr,"%s: %s does not exist\n",Prog_Name,argv[c]);
            exit (1);
          }

        if (stat(argv[argc-1],&B) == 0 && (B.st_mode & S_IFMT) == S_IFDIR)
          { DIR  = Strdup(argv[argc-1],NULL);
            ROOT = Strdup(root,NULL);
          }
        else
          { DIR  = PathTo(argv[2]);
            if (doh+dot+dop == 1)
              { if (doh)
                  ROOT = Root(argv[2],".hist");
                else if (dot)
                  ROOT = Root(argv[2],".ktab");
                else // dop
                  ROOT = Root(argv[2],".prof");
              }
            else
              ROOT = Root(argv[2],"");
          }

        yes = doh;
        if (doh && stat(Catenate(DIR,"/",ROOT,".hist"),&B) == 0)
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
        if (yes && stat(Catenate(dir,"/",root,".hist"),&B) == 0)
          { sprintf(command,"%s -f %s/%s.hist %s/%s.hist",op,dir,root,DIR,ROOT);
            system(command);
          }

        yes = dot;
        if (dot && stat(Catenate(DIR,"/",ROOT,".ktab"),&B) == 0)
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

        yes = dop;
        if (dop && stat(Catenate(DIR,"/",ROOT,".prof"),&B) == 0)
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
      }

    free(command);
  }

  exit (0);
}
