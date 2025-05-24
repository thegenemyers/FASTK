/********************************************************************************************
 *
 *  Print a BED-formatted file of the locations of FastK table k-mers in a specified FASTA file
 *
 *  Author  : Nancy Hansen
 *  Date    : January 2025
 *  Modified: Gene Myers, February 2025
 *
 *********************************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>

#define DEBUG

#include "gene_core.h"
#include "libfastk.h"

  //  Usage

static char *Usage =
              "[-vm] [-T<int(4)>] [-P<dir(/tmp)> <kmers>[.ktab] <target>[.\"dna\"] <out:bed>";

static int   VERBOSE;
static int   MERGE;
static int   KMER;

static void write_bed(char *aroot, char *proot, char *pname, FILE *fmap)
{ Profile_Index *RP;
  uint16        *prof;
  int64          pmax, plen;
  int            nscf;
  int            p, x;

  if (VERBOSE)
    fprintf(stderr,"\n  About to open relative kmer profile\n");

  RP = Open_Profiles(pname);
  if (RP == NULL)
    { fprintf(stderr,"\n%s: Cannot open/find FastK profile of %s relative to %s\n",
                     Prog_Name,aroot,proot);
      exit (1);
    }

  pmax = 20000;
  prof = Malloc(pmax*sizeof(uint16),"Profile array");
  if (prof == NULL)
    exit (1);

  nscf = RP->nreads;
  for (p = 0; p < nscf; p++)
    { int beg, end;

      plen = Fetch_Profile(RP,p,pmax,prof);
      if (plen > pmax)
        { pmax = 1.2*plen + 1000;
          prof = Realloc(prof,pmax*sizeof(uint16),"Profile array");
          if (prof == NULL)
            exit (1);
          Fetch_Profile(RP,p,pmax,prof);
        }

      if (MERGE)
        { beg = -1;
          end = -1;
          for (x = 0; x < plen; x++)
            if (prof[x] > 0)
              { if (x > end)
                  { if (beg >= 0)
                      fprintf(fmap,"%d\t%d\t%d\t%s\n",p,beg,end,proot);
                    beg = x;
                  }
                end = x+KMER;
              }
          if (beg >= 0)
            fprintf(fmap,"%d\t%d\t%d\t%s\n",p,beg,end,proot);
        }
      else
        for (x = 0; x <= plen; x++)
          if (prof[x] > 0)
            fprintf(fmap,"%d\t%d\t%d\t%s\n",p,x,x+KMER,proot);
    }

  Free_Profiles(RP);
}


/****************************************************************************************
 *
 *  Main Routine
 *
 *****************************************************************************************/

static int check_table(char *name, int lmer)
{ int   kmer;
  FILE *f;

  f = fopen(name,"r");
  if (f == NULL)
    { fprintf(stderr,"\n%s: Cannot find FastK table %s\n",Prog_Name,name);
      exit (1);
    }
  else
    { fread(&kmer,sizeof(int),1,f);
      if (lmer != 0 && kmer != lmer)
        { fprintf(stderr,"\n%s: Kmer (%d) of table %s != %d\n",Prog_Name,kmer,name,lmer);
          exit (1);
        }
      fclose(f);
      return (kmer);
    }
}

int main(int argc, char *argv[])
{ char  *KTAB, *ASM, *OUT;
  char  *AROOT, *KROOT;
  int    KUNIT;
  char  *KPROF;
  int    NTHREADS;
  FILE  *fmap;
  char  *SORT_PATH;

  char   command[5000];
  
  //  Command line processing

  { int    i, j, k;
    int    flags[128];

    char  *eptr;

    (void) eptr;

    ARG_INIT("KmerMap");

    NTHREADS  = 4;
    SORT_PATH = "/tmp";

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vm")
            break;
          case 'P':
            SORT_PATH = argv[i]+2;
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    MERGE   = flags['m'];

    if (argc != 4)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: verbose output to stderr\n");
        fprintf(stderr,"      -m: merge overlapping k-mer hits\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -T: number of threads to use\n");
        fprintf(stderr,"      -P: Place all temporary files in directory -P.\n");
        exit (1);
      }

    KTAB = argv[1];
    ASM  = argv[2];
    OUT  = argv[3];

    //  Remove any suffixes from argument names

    { char *suffix[] = { ".cram", ".bam", ".sam", ".db", ".dam",
                         ".fastq", ".fasta", ".fq", ".fa",
                         ".fastq.gz", ".fasta.gz", ".fq.gz",  ".fa.gz" };
      char *KPATH;
      int   len;

      for (j = 0; j < 13; j++)
        { len = strlen(ASM) - strlen(suffix[j]);
          if (len > 0 && strcmp(ASM+len,suffix[j]) == 0)
            break;
        }
      if (j >= 13)
        AROOT = strdup(ASM);
      else
        AROOT = Root(ASM,suffix[j]);

      KROOT = Root(KTAB,".ktab");
      KPATH = PathTo(KTAB);

      KMER = check_table(Catenate(KPATH,"/",KROOT,".ktab"),0);

      free(KPATH);

      if (VERBOSE)
        fprintf(stderr,"\n  Kmer size is %d\n",KMER);
    }
  }

  { char templateKP[20] = "._MQY_RP.XXXXXX";

    //  Create relative profile of k-mers against assembly

    if (VERBOSE)
      fprintf(stderr,"\n  Creating a profile of %s against the assembly %s\n",KROOT,AROOT);
  
    KUNIT = mkstemp(templateKP);
    KPROF = templateKP;

    if (VERBOSE)
      { sprintf(command,"FastK -v -T%d -P%s -k%d -p:%s %s -N%s",
                        NTHREADS,SORT_PATH,KMER,KTAB,ASM,KPROF);
        fprintf(stderr,"\n  Running command %s\n",command);
      }
    else
      sprintf(command,"FastK -T%d -P%s -k%d -p:%s %s -N%s",NTHREADS,SORT_PATH,KMER,KTAB,ASM,KPROF);

    SystemX(command);

    // Call write_bed to convert the kmer profile to a bed file

    if (MERGE)
      fmap = fopen(Catenate(OUT,".",AROOT,".kmers.merge.bed"),"w");
    else
      fmap = fopen(Catenate(OUT,".",AROOT,".kmers.bed"),"w");

    write_bed(AROOT, KROOT, KPROF, fmap);

    fclose(fmap);

    sprintf(command,"Fastrm %s",KPROF);
    SystemX(command);
  
    if (VERBOSE)
      printf("\n  Done\n\n");

    close(KUNIT);
  }

  free(KROOT);
  free(AROOT);

  Catenate(NULL,NULL,NULL,NULL);
  free(Prog_Name);

  exit (0);
}
