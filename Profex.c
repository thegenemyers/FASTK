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
#include "ONElib.h"

static char *Usage = "[-1Az] <source_root>[.prof] [ <read:int>[-(<read:int>|#)] ... ]";

static char *One_Schema =
  "P 3 prf               This is a 1-code fiel for profiles\n"
  "O P 1 8 INT_LIST      The profile count vector for the next read\n";

/****************************************************************************************
 *
 *  Test Stub
 *
 *****************************************************************************************/

int main(int argc, char *argv[])
{ Profile_Index *P;
  char          *command;
  int            ONE_CODE;
  int            ZFLAG;
  int            ASCII;

  //  Process options and capture command line for provenance

  { int    i, j, k;
    int    flags[128];

    (void) flags;

    ARG_INIT("Profex");

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

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("1Az")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    ONE_CODE = flags['1'];
    ZFLAG    = flags['z'];
    ASCII    = flags['A'];

    if (argc < 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -1: Produce 1-code as output.\n");
        fprintf(stderr,"      -A: tab-delimited ASCII as output.\n");
        fprintf(stderr,"      -z: Compress runs and ignore zeros.\n");
        exit (1);
      }

    if (ONE_CODE)
      { if (ASCII)
          fprintf(stderr,"%s: Warning, -1 overrides the -A flag\n",Prog_Name);
        ASCII = 0;
      }
  }

  P = Open_Profiles(argv[1]);
  if (P == NULL)
    { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,argv[1]);
      exit (1);
    }

  { int     i, c, p;
    int     id1, id2;
    char   *eptr, *fptr;
    uint16 *profile;
    int     pmax, plen;

    OneSchema *schema;
    OneFile   *file1;
    int64     *prof64;

    pmax    = 20000;
    profile = Malloc(pmax*sizeof(uint16),"Profile array");

    if (ONE_CODE)
      { schema = oneSchemaCreateFromText(One_Schema);
        file1  = oneFileOpenWriteNew("-",schema,"prf",1,1);
        oneAddProvenance(file1,"1.0","%s >?.prf",command);
        prof64 = Malloc(pmax*sizeof(int64),"Double Int Profile array");
      }

    for (c = 2; c <= argc; c++)
      { if (c == argc)
          { if (argc > 2)
              break;
            id1 = 1;
            id2 = P->nbase[P->nparts-1];
          }
        else
          { id1 = strtol(argv[c],&eptr,10);
            if (*eptr == '-')
              { if (eptr[1] == '#')
                  { id2  = P->nbase[P->nparts-1];
                    fptr = eptr+2;
                  }
                else
                  id2 = strtol(eptr+1,&fptr,10);
                if (*fptr != '\0')
                  { fprintf(stderr,"%s: argument '%s' is not an integer range\n",Prog_Name,argv[c]);
                    exit (1);
                  }
              }
            else
              { if (*eptr != '\0')
                  { fprintf(stderr,"%s: argument '%s' is not an integer\n",Prog_Name,argv[c]);
                    exit (1);
                  }
                id2 = id1;
              }
          }
        if (id1 > id2)
          { fprintf(stderr,"%s: range %s is empty!\n",Prog_Name,argv[c]);
            exit (1);
          }
        if (id1 <= 0 || id2 > P->nbase[P->nparts-1])
          { if (id1 == id2)
              fprintf(stderr,"%s: Id %d is out of range",Prog_Name,id1);
            else
              fprintf(stderr,"%s: Range %d-%d is out of range",Prog_Name,id1,id2);
            fprintf(stderr," [1,%lld]\n",P->nbase[P->nparts-1]);
            exit (1);
          }

        for (p = id1; p <= id2; p++)
          { plen = Fetch_Profile(P,(int64) p-1,pmax,profile);
            if (plen > pmax)
              { pmax    = 1.2*plen + 1000;
                profile = Realloc(profile,pmax*sizeof(uint16),"Profile array");
                Fetch_Profile(P,(int64) p-1,pmax,profile);
                if (ONE_CODE)
                  prof64 = Realloc(prof64,pmax*sizeof(int64),"Double Int Profile array");
              }
            if (ONE_CODE)
              { for (i = 0; i < plen; i++)
                  prof64[i] = profile[i];
                oneWriteLine(file1,'P',plen,prof64);
              }
            else if (ASCII)
              { printf("\nRead\t%d\n",p); 
                if (ZFLAG)
                  { int last = 0;
                    for (i = 0; i < plen; i++)
                      if (profile[i] != last)
                        { if (last != 0)
                            printf("\t%d\t%d\n",i+P->kmer-1,last);
                          if (profile[i] != 0)
                            printf("%d",i);
                          last = profile[i];
                        }
                    if (last != 0)
                      printf("\t%d\t%d\n",plen+P->kmer-1,last);
                  }
                else
                  for (i = 0; i < plen; i++)
                    printf("%d\n",profile[i]);
              }
            else
              { printf("\nRead %d:\n",p);
                if (ZFLAG)
                  { int last = 0;
                    for (i = 0; i < plen; i++)
                      if (profile[i] != last)
                        { if (last != 0)
                            printf(" - %5d (%d)\n",i+P->kmer-1,last);
                          if (profile[i] != 0)
                            printf(" %5d",i);
                          last = profile[i];
                        }
                    if (last != 0)
                      printf(" - %5d (%d)\n",plen+P->kmer-1,last);
                  }
                else
                  for (i = 0; i < plen; i++)
                    printf(" %5d: %5d\n",i,profile[i]);
              }
          }
      }

    if (ONE_CODE)
      { oneFileClose(file1);
        oneSchemaDestroy(schema);
      }
        
    free(profile);
  }

  Free_Profiles(P);

  free(command);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
