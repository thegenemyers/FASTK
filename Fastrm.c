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
#include <sys/stat.h>
#include <math.h>

#include "gene_core.h"

static char *Usage = "[-if] <source> ...";

int main(int argc, char **argv)
{ int QUERY;
  int QUIET;

  { int    i, j, k;
    int    flags[128];

    ARG_INIT("Fastrm")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("if")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    QUERY = flags['i'];
    QUIET = flags['f'];
    if (QUIET)
      QUERY = 0;

    if (argc < 2)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -i: prompt for each (stub) deletion\n");
        fprintf(stderr,"      -f: force operation quietly\n");
        exit (1);
      }
  }

  { int   c, len, yes, any, a;
    int   alen, rlen;
    int   doh, dot, dop;
    char *dir, *root, *fname;
    char *command;
    struct stat B;

    len = 0;
    for (c = 1; c < argc; c++)
      if ((int) strlen(argv[c]) > len)
        len = strlen(argv[c]);

    command = Malloc(3*len+50,"Allocating command buffer");

    for (c = 1; c < argc; c++)
      { dir  = PathTo(argv[c]);
        fname = rindex(argv[c],'/');
        if (fname == NULL)
          fname = argv[c];
        else
          fname += 1;
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

        any  = 0;
        if (doh && stat(Catenate(dir,"/",root,".hist"),&B) == 0)
          { yes = 1;
            any = 1;
            if (QUERY)
              { printf("remove %s/%s.hist? ",dir,root);
                yes = 0;
                while ((a = getc(stdin)) != '\n')
                  if (a == 'y' || a == 'Y')
                    yes = 1;
                  else if (a == 'n' || a == 'N')
                    yes = 0;
              }
            if (yes)
              unlink(Catenate(dir,"/",root,".hist"));
          }

        if (dot && stat(Catenate(dir,"/",root,".ktab"),&B) == 0)
          { yes = 1;
            any = 1;
            if (QUERY)
              { printf("remove %s/%s.ktab & hidden parts? ",dir,root);
                yes = 0;
                while ((a = getc(stdin)) != '\n')
                  if (a == 'y' || a == 'Y')
                    yes = 1;
                  else if (a == 'n' || a == 'N')
                    yes = 0;
              }
            if (yes)
              { sprintf(command,"rm -f %s/%s.ktab %s/.%s.ktab.*",dir,root,dir,root);
                system(command);
              }
          }

        if (dop && stat(Catenate(dir,"/",root,".prof"),&B) == 0)
          { yes = 1;
            any = 1;
            if (QUERY)
              { printf("remove %s/%s.prof & hidden parts? ",dir,root);
                yes = 0;
                while ((a = getc(stdin)) != '\n')
                  if (a == 'y' || a == 'Y')
                    yes = 1;
                  else if (a == 'n' || a == 'N')
                    yes = 0;
              }
            if (yes)
              { sprintf(command,"rm -f %s/%s.prof %s/.%s.pidx.* %s/.%s.prof.*",
                                dir,root,dir,root,dir,root);
                system(command);
              }
          }

        if (any == 0 && !QUIET)
          { if (doh+dot+dop == 3)
              fprintf(stderr,"%s: no FastK output files with root %s\n",Prog_Name,root);
            else
              fprintf(stderr,"%s: %s does not exist\n",Prog_Name,argv[1]);
          }

        free(root);
        free(dir);
      }

    free(command);
  }

  exit (0);
}
