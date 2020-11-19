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

static char *Usage = " [-h[<int(1)>:]<int(100)>] <source_root>[.hist]";

int main(int argc, char *argv[])
{ Histogram *H;
  int    HIST_LOW;
  int    HIST_HGH;

  //  Process arguments

  { int    i, j, k;
    int    flags[128];
    char  *eptr, *fptr;

    (void) flags;

    ARG_INIT("Histex")

    HIST_LOW    = 1;
    HIST_HGH    = 100;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
            break;
          case 'h':
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

    if (HIST_HGH > 0x7fff)
      HIST_HGH = 0x7fff;

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -h: Output histogram of counts in range given\n");
        exit (1);
      }
  }

  //  Load histogram into "cgram"

  H = Load_Histogram(argv[1]);
  if (H == NULL)
    { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,argv[1]);
      exit (1);
    }

  if (HIST_LOW < H->low || HIST_HGH > H->high)
    { fprintf(stderr,"%s: Range of histogram, [%d,%d], does not superset requested range\n",
                     Prog_Name,H->low,H->high); 
      exit (1);
    }

  Subrange_Histogram(H,HIST_LOW,HIST_HGH);

  //  Generate display

  { char       *root;
    int         j;
    int64       ssum, stotal;
    int64      *cgram;

    root = Root(argv[1],NULL);
    printf("\nHistogram of %d-mers of %s\n",H->kmer,root);
    free(root);

    cgram = H->hist;

    stotal = 0;
    for (j = HIST_LOW; j <= HIST_HGH; j++)
      stotal += cgram[j];

    printf("\n  Input: ");
    Print_Number(stotal,0,stdout);
    printf(" %d-mers\n",H->kmer);

    printf("\n     Freq:        Count   Cum. %%\n");
    ssum = 0;
    for (j = HIST_HGH; j > HIST_LOW; j--)
      { ssum += cgram[j];
        if (j == HIST_HGH)
          { if (ssum > 0)
              { printf(" >= %5d: %12lld",j,ssum);
                printf("   %5.1f%%\n",(100.*ssum)/stotal);
              }
          }
        else if (j < HIST_HGH && j > HIST_LOW && cgram[j] > 0)
          { printf("    %5d: %12lld",j,cgram[j]);
            printf("   %5.1f%%\n",(100.*ssum)/stotal);
          }
      }
    if (HIST_LOW > 1)
      printf(" <= %5d: %12lld   100.0%%\n",j,stotal-ssum);
    else
      printf("    %5d: %12lld   100.0%%\n",j,stotal-ssum);
  }

  Free_Histogram(H);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
