 /********************************************************************************************
 *
 *  FastK: a rapid disk-based k-mer counter for high-fidelity shotgun data sets.
 *     Uses a novel minimizer-based distribution scheme that permits problems of
 *     arbitrary size, and a two-staged "super-mer then weighted k-mer" sort to acheive
 *     greater speed when error rates are low (1% or less).  Directly produces read
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
#include <sys/resource.h>
#include <math.h>

#include "gene_core.h"
#include "FastK.h"

#ifdef DEVELOPER

int    DO_STAGE;   //  Which step to perform

#endif

static char *Usage[] = { "[-k<int(40)>] [-h[<int(1)>:]<int>] [-t<int(3)>] [-p] [-bc<int(0)>]",
                         "  [-v] [-T<int(4)>] [-P<dir(/tmp)>] [-M<int(12)>]",
                         "    <data::cram|[bs]am|f[ast][aq][.gz]|db|dam> ..."
                       };

  //  Option Settings

int    VERBOSE;      //  show progress
int    NTHREADS;     //  # of threads to run with
int64  SORT_MEMORY;  //  GB of memory for downstream KMcount sorts
char  *SORT_PATH;    //  where to put external files

int    KMER;         //  desired K-mer length
int    HIST_LOW;     // Start count for histogram
int       HIST_HGH;  // End count for histogram
int    DO_TABLE;     // Zero or table cutoff
int    DO_PROFILE;   // Do or not
int    BC_PREFIX;    // Ignore prefix of each read of this length

  //  Major parameters, sizes of things

int    NPARTS;       //  # of k-mer buckets to use
int    SMER;         //  size of a super-mer for sorts
int64  KMAX;         //  max k-mers in any part
int64  NMAX;         //  max super-mers in any part
int64  RMAX;         //  max run index for any thread

int    MOD_LEN;     //  length of minimizer buffer (power of 2)
int    MOD_MSK;     //  mask for minimzer buffer

int    MAX_SUPER;    //  = (KMER - PAD_LEN) + 1

int    SLEN_BITS;    //  # of bits needed to encode the length of a super-mer (less KMER-1)
uint64 SLEN_BIT_MASK;   //  Bit-mask for super-mer lengths

int RUN_BITS;     //  # of bits encoding largest run index
int RUN_BYTES;    //  # of bytes encoding largest run index

int SLEN_BYTES;   //  # of bytes encoding length of a super-mer
int PLEN_BYTES;   //  # of bytes encoding length of a compressed profile segment
int PROF_BYTES;   //  # of bytes encoding RUN_BYTES+PLEN_BYTES

int SLEN_BYTE_MASK;   //  Byte-mask for super-mer lengths

int KMER_BYTES;   //  # of bytes encoding a KMER
int SMER_BYTES;   //  # of bytes encoding a super-mer
int KMAX_BYTES;   //  # of bytes to encode an int in [0,KMAX)

int KMER_WORD;    //  bytes to hold a k-mer entry
int SMER_WORD;    //  bytes to hold a super-mer entry
int TMER_WORD;    //  bytes to hold a k-mer table entry
int CMER_WORD;    //  bytes to hold a count/index entry

int NUM_READS;    //  number of reads in dataset

int main(int argc, char *argv[])
{ char  *dbrt;

  { int    i, j, k;
    int    flags[128];
    int    memory; 
    char  *eptr, *fptr;

    ARG_INIT("FastK")

    KMER        = 40;
    SORT_MEMORY = 12000000000ll;
    NTHREADS    = 4;
    SORT_PATH   = "/tmp";
    HIST_LOW    = 0;
    HIST_HGH    = 0x7fff;
    DO_TABLE    = 0;
    BC_PREFIX   = 0;
#ifdef DEVELOPER
    DO_STAGE    = 0;
#endif

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vdp")
            break;
          case 'b':
            if (argv[i][2] != 'c')
              { fprintf(stderr,"%s: -%s is not a legal optional argument\n",Prog_Name,argv[i]);
                exit (1);
              }                
            argv[i] += 1;
            ARG_NON_NEGATIVE(BC_PREFIX,"Bar code prefiex")
            argv[i] -= 1;
            break;
          case 'k':
            ARG_POSITIVE(KMER,"K-mer length")
            break;
          case 't':
            if (argv[i][2] == '\0')
              { DO_TABLE = 3;
                break;
              }
            ARG_POSITIVE(DO_TABLE,"Cutoff for k-mer table")
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
          case 'M':
            ARG_POSITIVE(memory,"GB of memory for sorting step")
            SORT_MEMORY = memory * 1000000000ll;
            break;
          case 'P':
            SORT_PATH = argv[i]+2;
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
#ifdef DEVELOPER
          case '1':
            DO_STAGE = 1;
            break;
          case '2':
            DO_STAGE = 2;
            break;
          case '3':
            DO_STAGE = 3;
            break;
          case '4':
            DO_STAGE = 4;
            break;
#endif
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE    = flags['v'];   //  Globally declared in filter.h
    DO_PROFILE = flags['p'];

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, output statistics as proceed.\n");
        fprintf(stderr,"      -T: Use -T threads.\n");
        fprintf(stderr,"      -P: Place block level sorts in directory -P.\n");
        fprintf(stderr,"      -M: Use -M GB of memory in downstream sorting steps of KMcount.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -k: k-mer size.\n");
        fprintf(stderr,"      -h: Output histogram of counts in range given\n");
        fprintf(stderr,"      -t: Produce table of sorted k-mer & counts >= level specified\n");
        fprintf(stderr,"      -p: Produce read count profiles\n");
        fprintf(stderr,"     -bc: Ignore prefix of each read of given length (e.g. bar code)\n");
        exit (1);
      }
  }

  //  Get full path strong for sorting subdirectory (in variable SORT_PATH)

  { char  *cpath, *spath;
    DIR   *dirp;

    if (SORT_PATH[0] != '/')
      { cpath = getcwd(NULL,0);
        if (SORT_PATH[0] == '.')
          { if (SORT_PATH[1] == '/')
              spath = Catenate(cpath,SORT_PATH+1,"","");
            else if (SORT_PATH[1] == '\0')
              spath = cpath;
            else
              { fprintf(stderr,"%s: -P option: . not followed by /\n",Prog_Name);
                exit (1);
              }
          }
        else
          spath = Catenate(cpath,"/",SORT_PATH,"");
        SORT_PATH = Strdup(spath,"Allocating path");
        free(cpath);
      }
    else
      SORT_PATH = Strdup(SORT_PATH,"Allocating path");

    if ((dirp = opendir(SORT_PATH)) == NULL)
      { fprintf(stderr,"%s: -P option: cannot open directory %s\n",Prog_Name,SORT_PATH);
        exit (1);
      }
    closedir(dirp);
  }

  { Input_Partition *io;
    DATA_BLOCK      *block;
    int64            gsize;
    int              rsize, val;

    io = Partition_Input(argc,argv);

    dbrt = First_Root_Name(io);

    if (VERBOSE)
      fprintf(stderr,"\nDetermining minimizer scheme & partition for %s\n",dbrt);

    //  Determine number of buckets and padded minimzer scheme based on first
    //    block of the data set

    block = Get_First_Block(io,1000000000);

    KMER_BYTES = (KMER*2+7) >> 3;

    rsize  = KMER_BYTES + 2;
    gsize  = (block->totlen - KMER*block->nreads)*block->ratio*rsize;
    NPARTS = (gsize-1)/SORT_MEMORY + 1;

    if (VERBOSE)
      { double est = gsize/(1.*rsize);
        if (est >= 5.e8)
          fprintf(stderr,"  Estimate %.3fG %d-mers\n",est/1.e9,KMER);
        else if (est >= 5.e5)
          fprintf(stderr,"  Estimate %.3fM %d-mers\n",est/1.e6,KMER);
        else
          fprintf(stderr,"  Estimate %.3fK %d-mers\n",est/1.e3,KMER);
        if (NPARTS > 1)
          fprintf(stderr,"  Dividing data into %d parts\n",NPARTS);
      }

    MOD_LEN = 1;
    while (MOD_LEN < KMER)
      MOD_LEN <<= 1;
    MOD_LEN <<= 1;
    MOD_MSK = MOD_LEN-1;

    MAX_SUPER = Determine_Scheme(block);

    Free_First_Block(block);

    SMER = MAX_SUPER + KMER - 1;
  
    SLEN_BITS = 0;
    for (val = MAX_SUPER; val > 0; val >>= 1)
      SLEN_BITS += 1;
    SLEN_BIT_MASK = (0x1u << SLEN_BITS)-1;
    SLEN_BYTES    = (SLEN_BITS+7) >> 3;

    SMER_BYTES = (SMER*2+7) >> 3;
    SMER_WORD  = SMER_BYTES + SLEN_BYTES;
    KMER_WORD  = KMER_BYTES + 2;
    PLEN_BYTES = (SLEN_BITS+8) >> 3;
    TMER_WORD  = KMER_BYTES + 2;

    //  Make sure you can open (2 * NPARTS + 2) * NTHREADS + 3 files and then set up data structures
    //    for each such file

    { struct rlimit rlp;
      uint64        nfiles;

      nfiles = (2*NPARTS+2)*NTHREADS+3;
      getrlimit(RLIMIT_NOFILE,&rlp);
      if (nfiles > rlp.rlim_max)
        { fprintf(stderr,"%s: Cannot open %lld files simultaneously\n",Prog_Name,nfiles);
          exit (1);
        }
      rlp.rlim_cur = nfiles;
      setrlimit(RLIMIT_NOFILE,&rlp);
    }

#ifdef DEVELOPER
    if (DO_STAGE == 1)
#endif
      Split_Kmers(io,dbrt);

    Free_Input_Partition(io);
  }

#ifdef DEVELOPER
  if (DO_STAGE == 2)
#endif
    Sorting(".",dbrt);

  if (DO_TABLE > 0)
#ifdef DEVELOPER
    if (DO_STAGE == 3)
      Merge_Tables(".",dbrt);
#else
    Merge_Tables(".",dbrt);
#endif

  if (DO_PROFILE > 0)
#ifdef DEVELOPER
    if (DO_STAGE == 4)
      Merge_Profiles(".",dbrt);
#else
    Merge_Profiles(".",dbrt);
#endif

  free(SORT_PATH);

  Catenate(NULL,NULL,NULL,NULL);  //  frees internal buffers of these routines
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
