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

static char *Usage = " [-uAG] [-h[<int(1)>:]<int(100)>] <source_root>[.hist]";

int main(int argc, char *argv[])
{ Histogram *H;
  int    HIST_SET;
  int    HIST_LOW;
  int    HIST_HGH;
  int    UNIQUE;
  int    ASCII;
  int    GSCOPE;

  //  Process arguments

  { int    i, j, k;
    int    flags[128];
    char  *eptr, *fptr;

    (void) flags;

    ARG_INIT("Histex")

    HIST_SET    = 0;
    HIST_LOW    = 1;
    HIST_HGH    = 100;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("uAG")
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
    UNIQUE = flags['u'];
    GSCOPE = flags['G'];
    if (HIST_HGH > 0x7fff)
      HIST_HGH = 0x7fff;

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -h: Output histogram of counts in range given\n");
        fprintf(stderr,"      -u: Output histogram of unique k-mer counts (vs. instances)\n");
        exit (1);
      }

    if (GSCOPE)
      { if (ASCII || UNIQUE)
          fprintf(stderr,"%s: Warning, -G overrides both -A and -u flags\n",Prog_Name);
        ASCII = UNIQUE = 1;
        if (HIST_LOW != 1)
          fprintf(stderr,"%s: Warning: -G forces histogram range to start at 1\n",Prog_Name);
        HIST_LOW = 1;
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

  //  Generate display

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
            printf("%d %lld\n",j,hist[j]);
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

        printf("\n     Freq:        Count   Cum. %%\n");

        ssum = hist[HIST_HGH];
        if (ssum > 0)
          printf(" >= %5d: %12lld   %5.1f%%\n",HIST_HGH,ssum,(100.*ssum)/stotal);

        for (j = HIST_HGH-1; j > HIST_LOW; j--)
          { ssum += hist[j];
            if (hist[j] > 0)
              printf("    %5d: %12lld   %5.1f%%\n",j,hist[j],(100.*ssum)/stotal);
          }

        if (HIST_HGH > 1)
          { if (HIST_LOW == 1)
              printf("    %5d: %12lld   100.0%%\n",1,hist[1]);
            else
              printf(" <= %5d: %12lld   100.0%%\n",HIST_LOW,hist[HIST_LOW]);
          }
      }
  }

  Free_Histogram(H);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
