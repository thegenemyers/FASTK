/*********************************************************************************************\
 *
 *  Example code for opening and fetching compressed profiles produced by FastK
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

static char *Usage = "<source_root>[.prof] <read:int> ...";

/****************************************************************************************
 *
 *  Test Stub
 *
 *****************************************************************************************/

int main(int argc, char *argv[])
{ Profile_Index *P;

  { int    i, j, k;
    int    flags[128];

    (void) flags;

    ARG_INIT("Profex");

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
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

  P = Open_Profiles(argv[1]);
  if (P == NULL)
    { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,argv[1]);
      exit (1);
    }

  { int     c, id;
    char   *eptr;
    uint16 *profile;
    int     pmax, plen;

    pmax    = 20000;
    profile = Malloc(pmax*sizeof(uint16),"Profile array");

    for (c = 2; c < argc; c++)
      { id = strtol(argv[c],&eptr,10);
        if (*eptr != '\0' || argv[c][0] == '\0')
          { fprintf(stderr,"%s: argument '%s' is not an integer\n",Prog_Name,argv[c]);
             exit (1);
          }
        if (id <= 0 || id > P->nbase[P->nparts-1])
          { fprintf(stderr,"%s: Id %d is out of range\n",Prog_Name,id);
            exit (1);
          }
        plen = Fetch_Profile(P,(int64) id-1,pmax,profile);
        if (plen > pmax)
          { pmax    = 1.2*plen + 1000;
            profile = Realloc(profile,pmax*sizeof(uint16),"Profile array");
            Fetch_Profile(P,(int64) id-1,pmax,profile);
          }
        printf("\nRead %d:\n",id);
        for (int i = 0; i < plen; i++)
          printf(" %5d: %5d\n",i,profile[i]);
      }
    free(profile);
  }

  Free_Profiles(P);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
