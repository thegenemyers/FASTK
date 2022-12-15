/*******************************************************************************************
 *
 *  Example code for reading and displatying a kmer histogram produced by FastK.
 *
 *  Author:  Gene Myers
 *  Date  :  October 2020
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "libfastk.h"
#include "ONElib.h"

static char *Usage = "[-1] [-kAG] [-h[<int(1)>:]<int(100)>] <source_root>[.hist]";

static char *One_Schema =
  "P 5 khist               a histogram 1-code file\n"
  "D N 1 6 STRING          the name of the FastK .hist file this came from\n"
  "D R 2 3 INT 3 INT       the frequency range [low,hgh] covered\n"
  "O H 1 8 INT_LIST        a (hgh-low)+1 element list of the counts";

int main(int argc, char *argv[])
{ Histogram *H;
  char      *command;
  int    HIST_SET;
  int    HIST_LOW;
  int    HIST_HGH;
  int    UNIQUE;
  int    ASCII;
  int    GSCOPE;
  int    ONE_CODE;

  OneSchema *schema;
  OneFile   *file1;

  //  Process options and capture command line for provenance

  { int    i, j, k;
    int    flags[128];
    char  *eptr, *fptr;

    (void) flags;

    ARG_INIT("Histex")

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

    HIST_SET    = 0;
    HIST_LOW    = 1;
    HIST_HGH    = 100;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("kAG1")
            break;
          case 'h':
            HIST_SET = 1;
            HIST_LOW = strtol(argv[i]+2,&eptr,10);
            if (eptr > argv[i]+2)
              { if (HIST_LOW < 1 || HIST_LOW > 0x7fff)
                  { fprintf(stderr,"%s: Histogram count %d is out of range\n",
                                   Prog_Name,HIST_LOW);
                    exit (1);
                  }
                if (*eptr == ':')
                  { HIST_HGH = strtol(eptr+1,&fptr,10);
                    if (fptr > eptr+1 && *fptr == '\0')
                      { if (HIST_LOW > HIST_HGH)
                          { fprintf(stderr,"%s: Histogram range is invalid\n",Prog_Name);
                            exit (1);
                          }
                        break;
                      }
                  }
                else if (*eptr == '\0')
                  { HIST_HGH = HIST_LOW;
                    HIST_LOW = 1;
                    break;
                  }
              }
            fprintf(stderr,"%s: Syntax of -h option invalid -h[<int(1)>:]<int>\n",Prog_Name);
            exit (1);
        }
      else
        argv[j++] = argv[i];
    argc = j;

    ASCII  = flags['A'];
    UNIQUE = 1-flags['k'];
    GSCOPE = flags['G'];
    ONE_CODE = flags['1'];
    if (HIST_HGH > 0x7fff)
      HIST_HGH = 0x7fff;

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -h: Output histogram of counts in range given\n");
        fprintf(stderr,"      -k: Output histogram of k-mer instance counts (vs. unique k-mers)\n");
        fprintf(stderr,"      -A: Output in simple tab-delimited ASCII format\n");
        fprintf(stderr,"      -G: Output an ASCII format histogram especially for GeneScope.FK\n");
        fprintf(stderr,"      -1: Output in 1-code\n");
        exit (1);
      }

    if (ONE_CODE)
      { if (ASCII)
          fprintf(stderr,"%s: Warning, -1 overrides the -A flag\n",Prog_Name);
        ASCII = 0;
      }

    if (GSCOPE)
      { if (ASCII || !UNIQUE)
          fprintf(stderr,"%s: Warning, -G overrides both -A and -k flags\n",Prog_Name);
        ASCII = UNIQUE = 1;
        if (HIST_SET != 0)
          fprintf(stderr,"%s: Warning: -G forces histogram range to [1,1000]\n",Prog_Name);
        HIST_LOW = 1;
        HIST_HGH = 1000;
        HIST_SET = 1;
      }
  }

  //  Load histogram into "hist"

  H = Load_Histogram(argv[1]);
  if (H == NULL)
    { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,argv[1]);
      exit (1);
    }

  if (HIST_SET)
    { if (HIST_LOW < H->low || HIST_HGH > H->high)
        { fprintf(stderr,"%s: Range of histogram, [%d,%d], does not superset requested range\n",
                         Prog_Name,H->low,H->high); 
          exit (1);
        }
    }
  else
    { if (H->low > HIST_LOW)
        HIST_LOW = H->low;
      if (H->high < HIST_HGH)
        HIST_HGH = H->high;
    }

  Modify_Histogram(H,HIST_LOW,HIST_HGH,UNIQUE);

  if (ONE_CODE)
    { schema = oneSchemaCreateFromText(One_Schema);
      file1  = oneFileOpenWriteNew("-",schema,"khist",true,1);

      oneAddProvenance(file1,Prog_Name,"1.0","%s >?.khist",command);

      oneWriteLine(file1,'N',strlen(argv[1]),argv[1]);

      oneInt(file1,0) = HIST_LOW;
      oneInt(file1,1) = HIST_HGH;
      oneWriteLine(file1,'R',0,NULL);

      oneWriteLine(file1,'H',(HIST_HGH-HIST_LOW)+1,H->hist+HIST_LOW);

      oneFileClose(file1);
      oneSchemaDestroy(schema);
    }

  //  Generate display

  else
    { char       *root;
      int         j;
      int64       ssum, stotal;
      int64      *hist;

      hist = H->hist;
  
      if (ASCII)
        { if (GSCOPE)
            hist[HIST_HGH] = hist[HIST_HGH+2]/HIST_HGH;
          for (j = HIST_LOW; j <= HIST_HGH; j++)
            if (hist[j] > 0)
              printf("%d\t%lld\n",j,hist[j]);
        }
  
      else
        { root = Root(argv[1],NULL);
          if (UNIQUE)
            printf("\nHistogram of unique %d-mers of %s\n",H->kmer,root);
          else
            printf("\nHistogram of %d-mer instances of %s\n",H->kmer,root);
          free(root);
  
          stotal = 0;
      
          for (j = HIST_LOW; j <= HIST_HGH; j++)
            stotal += hist[j];
  
          printf("\n  Input: ");
          Print_Number(stotal,0,stdout);
          if (UNIQUE)
            printf(" unique %d-mers\n",H->kmer);
          else
            printf(" %d-mer instances\n",H->kmer);
  
          if (stotal == 0)
            printf("\n     Empty\n");
  
          else
            { printf("\n     Freq:        Count   Cum. %%\n");
  
              ssum = hist[HIST_HGH];
              if (ssum > 0)
                printf(" >= %5d: %12lld   %5.1f%%\n",HIST_HGH,ssum,(100.*ssum)/stotal);
      
              for (j = HIST_HGH-1; j > HIST_LOW; j--)
                { ssum += hist[j];
                  if (hist[j] > 0)
                    printf("    %5d: %12lld   %5.1f%%\n",j,hist[j],(100.*ssum)/stotal);
                }
      
              if (HIST_HGH > 1 && hist[HIST_LOW] > 0)
                { if (HIST_LOW == 1)
                     printf("    %5d: %12lld   100.0%%\n",1,hist[1]);
                  else
                    printf(" <= %5d: %12lld   100.0%%\n",HIST_LOW,hist[HIST_LOW]);
                }
            }
        }
    }

  Free_Histogram(H);

  free(command);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
