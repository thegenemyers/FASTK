/*******************************************************************************************
 *
 *  A module adapted from the VGPseq library (written by me).  This input module can
 *  read fasta, fastq, sam, bam, and cram files along with Dazzler db's and dam's.
 *  Given a target # of threads ITHREADS, it first finds partition points in the input file
 *  that split the file into the required # of parts:
 *
 *         Input_Partition *Partition_Input(int argc, char *argv[])
 *
 *  It can then either supply a single large training block:
 *
 *         DATA_BLOCK *Get_First_Block(Input_Partition *parts, int64 numbp)
 *
 *  or in a thread parallel manner read each partition and transmit it in DATA_BLOCKs to a
 *  given call-back routine:
 *
 *         void Scan_All_Input(Input_Partition *parts)
 *
 *  Author:    Gene Myers
 *  Date:      October 2020
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include <time.h>
#include <pthread.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <stdbool.h>

#include "libfastk.h"
#include "FastK.h"

#include "LIBDEFLATE/libdeflate.h"
#include "HTSLIB/htslib/hts.h"
#include "HTSLIB/htslib/hfile.h"
#include "HTSLIB/cram/cram.h"

#include <zlib.h>

#undef    DEBUG_FIND
#undef    DEBUG_IO
#undef    DEBUG_TRAIN
#undef    DEBUG_OUT

#undef    DEBUG_AUTO

#undef    DEBUG_BAM_IO
#undef    DEBUG_BAM_RECORD

#ifdef DEBUG_OUT
#define CALL_BACK  Print_Block
#else
#define CALL_BACK  Distribute_Block
#endif

#define IO_BLOCK 10000000ll

#define DT_BLOCK  1000000ll
#define DT_MINIM   100000ll
#define DT_READS    10000

typedef struct libdeflate_decompressor DEPRESS;

#define INT_MAXLEN 10

#define CRAM   0
#define BAM    1
#define SAM    2
#define FASTQ  3
#define FASTA  4
#define DAZZ   5

static char *Tstring[7] = { "cram", "bam", "sam", "fastq", "fasta", "db", "dam" };

typedef struct
  { char      *path;   //  Full path name
    char      *root;   //  Ascii root name
    char      *pwd;    //  Ascii path name
    int64      fsize;  //  Size of file in bytes
    int        ftype;  //  Type of file
    int        zipd;   //  Is file gzip'd

    int64      zsize;  //  Size of compression index for cram and dazz
    int64     *zoffs;  //  zoffs[i] = offset to compressed block i
    int        DB_all; //  Trim parameters for Dazzler DB's
    int        DB_cut;
  } File_Object;

typedef struct
  { int64  fpos;   //  offset in source file of desired unit (gzip block, cram contaianer)
    uint32 boff;   //  offset within a unit block of desired record
  } Location;

#define SAMPLE 0
#define SPLIT  1

typedef struct
  { File_Object *fobj;   //  array of file objects
    uint8       *buf;    //  block buffer
    int          fid;    //  fid of fobj[bidx].path
    int          bidx;   //  Scan range is [bidx:beg,eidx:end)
    Location     beg; 
    int          eidx;
    Location     end;

    int          action; //  sampling or splitting
    int64        work;   //  total sizes of all files
    DATA_BLOCK   block;  //  block buffer for this thread
    int          thread_id;  //  thread index (in [0,NTHREAD-1]
    void *     (*output_thread)(void *);  //  output routine for file type
    DEPRESS     *decomp; //    decompressor
  } Thread_Arg;

static void do_nothing(Thread_Arg *parm)
{ (void) parm; }



/*******************************************************************************************
 *
 *  Routines to open and close the sequence file types supported
 *
 ********************************************************************************************/

  //  Open and get info about each input file

static int64 *genes_cram_index(char *path, int64 fsize, int64 *zsize);
static void   read_DB_stub(char *path, int *cut, int *all);
static int64 *get_dazz_offsets(FILE *idx, int64 *zsize);

static void Fetch_File(char *arg, File_Object *input, int gz_ok)
{ static char *suffix[] = { ".cram", ".bam", ".sam", ".db", ".dam",
                            ".fastq", ".fasta", ".fq", ".fa",
                            ".fastq.gz", ".fasta.gz", ".fastq", ".fasta",
                            ".fq.gz",  ".fa.gz", ".fq", ".fa" };
  static char *extend[] = { ".cram", ".bam", ".sam", ".db", ".dam",
                            ".fastq", ".fasta", ".fq", ".fa",
                            ".fastq.gz", ".fasta.gz", ".fastq.gz", ".fasta.gz",
                            ".fq.gz",  ".fa.gz", ".fq.gz", ".fa.gz" };

  struct stat stats;
  char  *pwd, *root, *path, *temp;
  int    fid, i;
  int64  fsize, zsize, *zoffs;
  int    ftype, zipd;

  pwd = PathTo(arg);
  for (i = 0; i < 17; i++)
    { root  = Root(arg,suffix[i]);
      fid   = open(Catenate(pwd,"/",root,extend[i]),O_RDONLY);
      if (fid >= 0) break;
      free(root);
    }
  if (fid < 0)
    { fprintf(stderr,"\n%s: Cannot open %s as a .cram|[bs]am|f{ast}[aq][.gz]|db|dam file\n",
                     Prog_Name,arg);
      Clean_Exit(1);
    }

  if (i == 3 || i == 4)
    ftype = DAZZ;
  else if (i >= 5)
    ftype = FASTA - (i%2);
  else
    ftype = i;
  zipd = (i >= 9);

  path = Strdup(Catenate(pwd,"/",root,extend[i]),"Allocating full path name");

  zoffs = NULL;
  if (zipd)
    { if (gz_ok)
        { if (stat(path,&stats) == -1)
            { fprintf(stderr,"\n%s: Cannot get stats for %s\n",Prog_Name,path);
              Clean_Exit(1);
            }
          fsize = stats.st_size;
          zsize = (fsize-1)/IO_BLOCK+1;
        }
      else
        { if (VERBOSE)
            fprintf(stderr,"  Gzipped file %s being temporarily uncompressed\n",arg);
          temp = Strdup(Catenate(SORT_PATH,"/",root,extend[i]),"Allocating full path name");
          temp[strlen(temp)-3] = '\0';
          system(Catenate("gunzip -c ",path," >",temp));
          free(path);
          path = temp;
          zipd  = 0;
          if (stat(path,&stats) == -1)
            { fprintf(stderr,"\n%s: Cannot get stats for %s\n",Prog_Name,path);
              Clean_Exit(1);
            }
          fsize = stats.st_size;
          zsize = (fsize-1)/IO_BLOCK+1;
        }
    }
  else if (ftype == DAZZ)
    { FILE *idx;

      read_DB_stub(path,&(input->DB_cut),&(input->DB_all));

      free(path);
      close(fid);

      path = Strdup(Catenate(pwd,"/.",root,".bps"),"Allocating full path name");
      fid = open(path,O_RDONLY);
      if (fid < 0)
        { fprintf(stderr,"\n%s: Dazzler DB %s does not have a .bps file ! ?\n",Prog_Name,root);
          Clean_Exit(1);
        }
      if (fstat(fid, &stats) == -1)
        { fprintf(stderr,"\n%s: Cannot get stats for %s\n",Prog_Name,arg);
          Clean_Exit(1);
        }
      fsize = stats.st_size;
      idx = fopen(Catenate(pwd,"/.",root,".idx"),"r");
      if (idx == NULL)
        { fprintf(stderr,"\n%s: Dazzler DB %s does not have a .idx file ! ?\n",Prog_Name,root);
          Clean_Exit(1);
        }
      zoffs = get_dazz_offsets(idx,&zsize);
      zoffs[zsize] = fsize;
      fclose(idx);
    }
  else
    { if (fstat(fid, &stats) == -1)
        { fprintf(stderr,"\n%s: Cannot get stats for %s\n",Prog_Name,arg);
          Clean_Exit(1);
        }
      fsize = stats.st_size;
      if (ftype == CRAM)
        zoffs = genes_cram_index(path,fsize,&zsize);
      else
        zsize = (fsize-1)/IO_BLOCK+1;
    }

#ifdef DEBUG_FIND
  fprintf(stderr,"\n%s is a %s file, ftype = %d, zipd = %d, fsize = %lld\n",
                 arg,suffix[i],ftype,zipd,fsize);
#endif

  close(fid);

  input->path  = path;
  input->root  = root;
  input->pwd   = pwd;
  input->fsize = fsize;
  input->ftype = ftype;
  input->zipd  = zipd;
  input->zsize = zsize;
  input->zoffs = zoffs;
}

static void Free_File(File_Object *input, int gz_ok)
{ if (!gz_ok && input->zipd)
    unlink(input->path);
  free(input->zoffs);
  free(input->path);
  free(input->root);
  free(input->pwd);
}

static void Free_Gunzips( File_Object *input, int n)
{ int i;

  for (i = 0; i < n; i++)
    if (input[i].zipd)
      unlink(input[i].path);
}


/*******************************************************************************************
 *
 *  DATA_BLOCK CODE
 *
 ********************************************************************************************/

static int homo_compress(char *seq, int len)
{ int i, x, n;

  n = 0;
  x = seq[n++];
  for (i = 1; i < len; i++)
    if (seq[i] != x)
      seq[n++] = x = seq[i];
  seq[n] = '\0';
  return (n);
} 
  
static int Add_Data_Block(DATA_BLOCK *dset, int len, char *seq)
{ int   n, r;
  int64 o;

  r = 0;
  n = dset->nreads;
  o = dset->boff[n];

  if (dset->rem > 0)
    { len = dset->rem;
      seq = dset->next;
    }
  else if (COMPRESS)
    len = homo_compress(seq,len);

  if (n >= dset->maxrds)
    return (1);
  if (o+len >= dset->maxbps)
    { if (dset->maxbps-o > DT_MINIM)
        { r = 1;
          dset->rem  = len;
          len = (dset->maxbps-o);
          dset->rem -= (len - (KMER-1));
          dset->next = seq + (len - (KMER-1)); 
        }
      else
        return (1);
    }
  else
    dset->rem = 0;
  dset->nreads = ++n;
  memcpy(dset->bases+o,seq,len);
  o += len;
  dset->bases[o] = '\0';
  dset->boff[n] = o+1;
  dset->totlen += len;
  return (r);
}

#if defined(DEBUG_OUT) || defined(DEBUG_TRAIN)

static void Print_Block(DATA_BLOCK *dset, int tid)
{ int64 end;
  int   i;
  char *bases;
  static int tc = 0;

  (void) tid;

  printf("Buffer load rem = %d\n",dset->rem);
  bases = dset->bases;
  for (i = 0; i < dset->nreads; i++)
    { end = dset->boff[i+1]-1;
      printf("%5d %8d:\n",i,tc++);
#undef ONE_LINE
#ifdef ONE_LINE
      printf("    %.80s\n",bases+dset->boff[i]);
#else
      { int64 j;

        for (j = dset->boff[i]; j+80 < end; j += 80)
          printf("    %.80s\n",bases+j);
        if (j < end)
          printf("    %.*s\n",(int) (end-j),bases+j);
       }
#endif
    }
  fflush(stdout);
}

#endif

static void Reset_Data_Block(DATA_BLOCK *dset, int roll)
{ if (roll)
    { memcpy(dset->bases,dset->bases+(dset->boff[dset->nreads]-KMER),KMER-1);
      dset->totlen = KMER-1;
    }
  else
    dset->totlen = 0;
  dset->nreads = 0;
  dset->boff[0] = 0;
}


/*******************************************************************************************
 *
 *  FASTA / FASTQ SPECIFIC CODE
 *
 ********************************************************************************************/

/*******************************************************************************************
 *
 *  fast_nearest: Routine to find next entry given arbitray start point.
 *
 ********************************************************************************************/

  //  Automata states (buffer is processed char at a time)

#define UK    0
#define FK    1
#define FK1   2

#if defined(DEBUG_AUTO)

static char *Name[] =
  { "UK", "FK", "FK1", "QHL" };

#endif

  //  Check blocks [beg,end) starting with the first complete entry and proceeding
  //     into blocks end ... if necessary to complete the last entry started in this
  //     range.

static void fast_nearest(Thread_Arg *data)
{ int          fid    = data->fid;
  int64        beg    = data->beg.fpos;
  uint8       *buf    = data->buf;
  File_Object *inp    = data->fobj + data->bidx;
  int          fastq = (inp->ftype == FASTQ);

  int64 blk, off;
  int   slen;
  int   state;
  int   nl_1, nl_2, pl_1, alive;
  int   b, c;

#ifdef DEBUG_FIND
  fprintf(stderr,"\nFind starting at %lld\n",beg);
#endif

  blk = beg / IO_BLOCK;
  off = beg % IO_BLOCK;

  lseek(fid,blk*IO_BLOCK,SEEK_SET);

  if (fastq)
    state = UK;
  else
    state = FK;
  nl_1 = pl_1 = nl_2 = 1;
  alive = 0;

#ifdef DEBUG_FIND
  fprintf(stderr,"\nFrom block %lld / offset %lld\n",blk,off);
#endif

  while (blk < inp->zsize)
    {
#ifdef DEBUG_FIND
      fprintf(stderr,"  Loading block %lld: @%lld",blk,lseek(fid,0,SEEK_CUR));
#endif
      slen = read(fid,buf,IO_BLOCK);
#ifdef DEBUG_FIND
      fprintf(stderr," %d %d\n",slen,errno);
#endif

      for (b = off; b < slen; b++)
        { c = buf[b];
#ifdef DEBUG_AUTO
          fprintf(stderr,"  %.5s: %c\n",Name[state],c);
#endif
          switch (state)

          { case UK:
              if (c == '@' && alive)
                { data->beg.fpos = blk*IO_BLOCK + b;
                  data->beg.boff = 0;
                  return;
                }
              alive = (c == '\n' && ! (nl_2 || pl_1));
              nl_2 = nl_1;
              nl_1 = (c == '\n');
              pl_1 = (c == '+');
#ifdef DEBUG_AUTO
              fprintf(stderr,"    n2=%d n1=%d p1=%d a=%d\n",nl_2,nl_1,pl_1,alive);
#endif
              break;

            case FK:
              if (c == '\n')
                state = FK1;
              break;

            case FK1:
              if (c == '>')
                { data->beg.fpos = blk*IO_BLOCK + b;
                  data->beg.boff = 0;
                  return;
                }
              else if (c != '\n')
                state = FK;
              break;
           }
        }

      blk += 1;
      off = 0;
    }

  data->beg.fpos = -1;
  data->beg.boff = 0;
  return;
}


/*******************************************************************************************
 *
 *  Second pass thread routines prepare .seq formated output in buffers for each input block.
 *  Note carefully that there are separate blocks for the forward and reverse files if
 *  processing a pair, the reads in the blocks need to be paired.
 *
 ********************************************************************************************/

  //  Automata states (buffer is processed char at a time)

#define QAT     0
#define HSKP    1
#define QSEQ    2
#define QPLS    3
#define QSKP    4
#define ASEQ    5
#define AEOL    6

#if defined(DEBUG_AUTO)

static char *Name2[] =
  { "QAT", "HEAD", "HSKP", "QSEQ", "QPLS", "QSKP", "QQVS", "ASEQ", "AEOL" };

#endif

#define DUMP(scale,notyet,dclose)					  \
if (action == SAMPLE)                                                     \
  { dset->ratio = ((1.*parm->work) / (totread-(notyet))) * scale;	  \
    dset->totlen = dset->boff[dset->nreads] - dset->nreads;		  \
    dclose(fid);                                                          \
    return (NULL);                                                        \
  }                                                                       \
else                                                                      \
  { dset->totlen = dset->boff[dset->nreads] - dset->nreads;		  \
    CALL_BACK(dset,tid);                                                  \
    if (CLOCK)                                                            \
      { cumbps += dset->totlen;                                           \
        if (cumbps >= nxtbps)                                             \
          { fprintf(stderr,"\r  %3d%%",(int) ((100.*cumbps)/estbps));     \
            fflush(stderr);                                               \
            nxtbps = cumbps+pct1;                                         \
          }                                                               \
      }                                                                   \
  }                                                                       \


#define END_SEQ(notyet)								\
{ line[olen++] = lastc = '\0';							\
  dset->boff[++dset->nreads] = olen;						\
  if (olen > omax-DT_MINIM || dset->nreads >= dset->maxrds)			\
    { DUMP(ratio,notyet,close)							\
      Reset_Data_Block(dset,0);							\
      olen = 0;									\
    }										\
}

#define ADD(c)									\
{ if (!COMPRESS || c != lastc)							\
    { if (olen >= omax)								\
        { line[olen++] = '\0';							\
          dset->boff[++dset->nreads] = olen;					\
          dset->rem  = 1;							\
          DUMP(ratio,slen-b,close)						\
          dset->rem  = 0;							\
          Reset_Data_Block(dset,1);						\
          olen = KMER-1;							\
        }									\
      line[olen++] = lastc = c;							\
    }										\
}

  //  Write fast records in relevant partition

static void *fast_output_thread(void *arg)
{ Thread_Arg  *parm   = (Thread_Arg *) arg;
  File_Object *fobj   = parm->fobj;
  uint8       *buf    = parm->buf;
  int          action = parm->action;
  DATA_BLOCK  *dset   = &parm->block;
  int          tid    = parm->thread_id;
  int          fastq  = (fobj->ftype == FASTQ);    //  Is the same for all files (already checked)

  File_Object *inp;
  int          f, fid;
  gzFile       gzid = NULL;
  int64        blk, off;
  int64        epos, eblk, eoff;
  int64        totread;
  double       ratio;

  int   state, lastc;
  int   omax, olen;
  char *line;

  int64 estbps, cumbps, nxtbps, pct1;
  int   CLOCK;

  estbps = nxtbps = pct1 = cumbps = 0;
  if (VERBOSE && tid == 0 && action != SAMPLE)
    { estbps = dset->ratio / ITHREADS;
      nxtbps = pct1 = estbps/100;
      fprintf(stderr,"\n    0%%");
      fflush(stderr);
      CLOCK = 1;
    }
  else
    CLOCK = 0;

  //  Do relevant section of each file assigned to this thread in sequence

  omax = dset->maxbps;
  line = dset->bases;

  totread = 0;
  olen = 0;
  ratio = 1.;

  for (f = parm->bidx; f <= parm->eidx; f++)
    { inp = fobj+f;
      fid = open(inp->path,O_RDONLY);
      if (f < parm->eidx)
        epos = inp->fsize;
      else
        epos = parm->end.fpos;
      if (f > parm->bidx)
        parm->beg.fpos = 0;

#ifdef DEBUG_IO
      fprintf(stderr,"Block: %12lld to %12lld --> %8lld\n",
                     parm->beg.fpos,epos,epos - parm->beg.fpos);
      fflush(stdout);
#endif

      blk   = parm->beg.fpos / IO_BLOCK;
      off   = parm->beg.fpos % IO_BLOCK;
      eblk  = (epos-1) / IO_BLOCK;
      eoff  = (epos-1) % IO_BLOCK + 1;

      if (inp->zipd)
        { gzid = gzdopen(fid,"r");
          eblk = 0x7fffffffffffffffll;
        }
      else
        lseek(fid,blk*IO_BLOCK,SEEK_SET);

      state = QAT;
      lastc = 0;

#ifdef DEBUG_IO
      fprintf(stderr,"\nFrom block %lld / offset %lld\n",blk,off);
#endif

      while (blk <= eblk)
        { int c, b, slen;

#ifdef DEBUG_IO
          fprintf(stderr,"  Loading block %lld: @%lld",blk,lseek(fid,0,SEEK_CUR));
#endif
          if (inp->zipd)
            { slen = gzread(gzid,buf,IO_BLOCK);
              if (blk == 0 && slen != 0 && action == SAMPLE)
                ratio = (1.*slen) / gzoffset(gzid);
            }
          else
            slen = read(fid,buf,IO_BLOCK);

          if (slen == 0)
            break;

          totread += slen;
#ifdef DEBUG_IO
          fprintf(stderr," %d\n",slen);
#endif

          if (blk == eblk && eoff < slen)
            slen = eoff;

          for (b = off; b < slen; b++)
            { c = buf[b];
#ifdef DEBUG_AUTO
              fprintf(stderr,"  %.5s: %c\n",Name2[state],c);
#endif
              switch (state)

              { case QAT:
                  state = HSKP;
                  break;

                case HSKP:
                  if (c == '\n')
                    { if (fastq)
                        state = QSEQ;
                      else
                        state = ASEQ;
                    }
                  break;
                  
                case QSEQ:
                  if (c != '\n')
                    ADD(c)
                  else
                    { END_SEQ(slen-b)
                      state = QPLS;
                    }
                  break;

                case QPLS:
                  if (c == '\n')
                    state = QSKP;
                  break;

                case QSKP:
                  if (c == '\n')
                    state = QAT;
                  break;

                case AEOL:
                  if (c == '>')
                    { END_SEQ(slen-b)
                      state = HSKP;
                    }
                  else if (c != '\n')
                    { ADD(c)
                      state = ASEQ;
                    }
                  break;

                case ASEQ:
                  if (c == '\n')
                    state = AEOL;
                  else
                    ADD(c)
              }
            }
          blk += 1;
          off = 0;
        }
      if (state == AEOL)
        END_SEQ(0)

      if (inp->zipd)
        gzclose_r(gzid);
      else
        close(fid);
    }

  dset->totlen = dset->boff[dset->nreads] - dset->nreads;
  if (action == SAMPLE)
    { dset->ratio = ((1.*parm->work) / totread) * ratio;
      return (NULL);
    }
  else if (dset->nreads > 0)
    CALL_BACK(dset,tid);

  if (CLOCK)
    fprintf(stderr,"\r         \r");

  return (NULL);
}


/*******************************************************************************************
 *
 *  BAM/SAM SPECIFIC CODE
 *
 ********************************************************************************************/

#define BAM_BLOCK  0x10000
#define HEADER_LEN      36
#define SEQ_RUN         40

static int DNA[256] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
    0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
  };

 //  Get value of little endian integer of n-bytes

static inline uint32 getint(uint8 *buf, int n)
{ uint32 val;
  int    k;

  val = 0;
  for (k = n-1; k >= 0; k--)
    val = (val << 8) | buf[k];
  return (val);
}

 //  Next len chars are printable and last is zero?

static inline int valid_name(char *name, int len)
{ int i;

  if (len < 1)
    return (0);
  len -= 1;
  for (i = 0; i < len; i++)
    if ( ! isprint(name[i]))
      return (0);
  return (name[len] == 0);
}

 //  Next len ints are plausible cigar codes

static inline int valid_cigar(int32 *cigar, int len, int range)
{ int i, c;

  for (i = 0; i < len; i++)
    { c = cigar[i];
      if ((c & 0xf) > 8)
        return (0);
      c >>= 4;
      if (c < 0 || c > range)
        return (0);
    }
  return (1);
}


/*******************************************************************************************
 *
 *  Routines to manage a BAM_FILE stream object
 *
 ********************************************************************************************/

  //   Open BAM stream where compressed block start is known.  Compressed BAM blocks are buffered
  //     in a very big IO buffer and the current uncompressed block is in a 64Kbp array.

typedef struct
  { int       fid;              //  file descriptor
    int       last;             //  last block of data in file has been read
    uint8    *buf;              //  IO buffer (of IO_BLOCK bytes, supplied by caller)
    int       blen;             //  # of bytes currently in IO buffer
    int       bptr;             //  start of next BAM block in IO buffer
    uint8     bam[BAM_BLOCK+1]; //  uncompressed bam block
    uint32    bsize;            //  length of compressed bam block
    uint32    ssize;            //  length of uncompressed bam block
    Location  loc;              //  current location in bam file
    DEPRESS  *decomp;
  } BAM_FILE;

  //  Load next len bytes of uncompressed BAM data into array data

static void bam_get(BAM_FILE *file, uint8 *data, int len)
{ int    chk, off;
  uint8 *bam  = file->bam;
  int    boff = file->loc.boff;

  off = 0;
  chk = file->ssize - boff;
  while (len >= chk)
    { int    bptr, blen, bsize;
      uint8 *block, *buf;
      uint32 ssize;
      size_t tsize;

#ifdef DEBUG_BAM_IO
      printf("Move %d bytes to data+%d from bam+%d\n",chk,off,boff);
#endif
      if (data != NULL)
        memcpy(data+off,bam+boff,chk);
      off += chk;
      len -= chk;

      file->loc.fpos += file->bsize;
#ifdef DEBUG_BAM_IO
      printf("File pos %lld\n",file->fpos);
#endif

      if (chk == 0 && len == 0)
        { file->loc.boff = 0;
          return;
        }

      buf   = file->buf;
      bptr  = file->bptr;
      blen  = file->blen;
      block = buf+bptr;
#ifdef DEBUG_BAM_IO
      printf("Block at buffer %d\n",bptr);
#endif
      while (bptr + 18 > blen || bptr + (bsize = getint(block+16,2) + 1) > blen)
        { chk = blen-bptr;
          if (file->last)
            { fprintf(stderr,"\n%s: Corrupted BAM file\n",Prog_Name);
              Clean_Exit(1);
            }
          memmove(buf,block,chk);
          blen = chk + read(file->fid,buf+chk,IO_BLOCK-chk);
#ifdef DEBUG_BAM_IO
          printf("Loaded %d to buf+%d for a total of %d\n",IO_BLOCK-chk,chk,blen);
#endif
          if (blen < IO_BLOCK)
            file->last = 1;
          file->blen = blen;
          bptr = 0;
          block = buf;
        }

      //  Fetch and uncompress next Bam block

      if (libdeflate_gzip_decompress(file->decomp,block,bsize,bam,BAM_BLOCK,&tsize) != 0)
        { fprintf(stderr,"\n%s: Bad gzip block\n",Prog_Name);
          Clean_Exit(1);
        }
      ssize = tsize;
      boff = 0;
#ifdef DEBUG_BAM_IO
      printf("Loaded gzip block of size %d into %d\n",bsize,ssize);
#endif

      file->bsize = bsize;
      file->ssize = ssize;
      file->bptr  = bptr + bsize;
      chk = ssize;
    }

#ifdef DEBUG_BAM_IO
  printf("Xfer %d bytes to data+%d from bam+%d\n",len,off,boff);
#endif
  if (data != NULL)
    memcpy(data+off,bam+boff,len);
  file->loc.boff = boff+len;
}

  //  Startup a bam stream, the location must be valid.

static void bam_start(BAM_FILE *file, int fid, uint8 *buf, Location *loc)
{ file->fid   = fid;
  file->ssize = 0;
  file->bsize = 0;
  file->buf   = buf;
  file->bptr  = 0;
  file->blen  = 0;
  file->last  = 0;
  lseek(fid,loc->fpos,SEEK_SET);
  file->loc.fpos = loc->fpos;
  file->loc.boff = 0;
  bam_get(file,NULL,loc->boff);
}

static int bam_eof(BAM_FILE *file)
{ return (file->loc.boff == file->ssize && file->bptr == file->blen && file->last); }


/*******************************************************************************************
 *
 *  Routines to manage a SAM stream, but as a BAM_FILE (only select fields are used)
 *
 ********************************************************************************************/

  //  Get next line of SAM input if possible

static uint8 *sam_getline(BAM_FILE *file)
{ int    rem;
  int    blen = file->blen;
  int    bptr = file->bptr;
  uint8 *buf  = file->buf;
  uint8  *b, *d;

  b = buf + bptr;
  rem = blen-bptr;
  if (rem == 0)
    d = NULL;
  else
    d = memchr(b,'\n',rem);
  if (d == NULL)
    { if (file->last)
        { fprintf(stderr,"\n%s: Corrupted SAM file",Prog_Name);
          Clean_Exit(1);
        }
      memmove(buf,buf+bptr,rem);
      blen = rem + read(file->fid,buf+rem,IO_BLOCK-rem);
      if (blen < IO_BLOCK)
        file->last = 1;
      file->blen = blen;
      bptr = 0;
      b = buf;
      d = memchr(b,'\n',blen);
      if (d == NULL)
        { if (blen < IO_BLOCK)
            fprintf(stderr,"\n%s: Corrupted SAM file",Prog_Name);
          else
            fprintf(stderr,"\n%s: SAM-line is longer than max %lld\n",Prog_Name,IO_BLOCK);
          Clean_Exit(1);
        }
    }
  d += 1;
  file->bptr = d-buf;
  file->loc.fpos += d-b;
  return (b);
}

  //  Startup a sam stream

static void sam_start(BAM_FILE *file, int fid, uint8 *buf, Location *loc)
{ file->fid   = fid;
  file->buf   = buf;
  file->bptr  = 0;
  file->blen  = 0;
  file->last  = 0;
  lseek(fid,loc->fpos,SEEK_SET);
  file->loc.fpos = loc->fpos;
  file->loc.boff = 0;
}


/*******************************************************************************************
 *
 *  Routines to find bam blocks and valid locations to start scan threads for first pass
 *
 ********************************************************************************************/

  //  Find first record location (skip over header) in parm->fid
  //    Return value is in parm->beg

static void skip_bam_header(Thread_Arg *parm)
{ uint8    *buf = parm->buf;
  int       fid = parm->fid;

  BAM_FILE  _bam, *bam = &_bam;
  Location  zero = { 0ll, 0 };
  uint8     data[4];
  int       i, ntxt, ncnt, nlen;

  //  At start of file so can use BAM stream

  bam->decomp = parm->decomp;
  bam_start(bam,fid,buf,&zero);

#ifdef DEBUG_FIND
  fprintf(stderr,"Header seek\n");
  fflush(stderr);
#endif

  bam_get(bam,data,4);
  if (memcmp(data,"BAM\1",4) != 0)
    { fprintf(stderr,"\n%s: Corrupted BAM header %.4s\n",Prog_Name,data);
      Clean_Exit(1);
    }

  bam_get(bam,data,4);
  ntxt = getint(data,4);
  bam_get(bam,NULL,ntxt);

  bam_get(bam,data,4);
  ncnt = getint(data,4);
  for (i = 0; i < ncnt; i++)
    { bam_get(bam,data,4);
      nlen = getint(data,4);
      bam_get(bam,NULL,nlen+4);
    }

  parm->beg = bam->loc;

#ifdef DEBUG_FIND
  fprintf(stderr,"  Begin @ %lld/%d\n",parm->beg.fpos,parm->beg.boff);
  fflush(stderr);
#endif
}

  //  Find next identifiable entry location forward of parm->fpos in parm->fid
  //    Return value is in parm->beg

static void bam_nearest(Thread_Arg *parm)
{ uint8       *buf  = parm->buf;
  int          fid  = parm->fid;
  DEPRESS     *decomp = parm->decomp;
  int64        fpos = parm->beg.fpos;

  BAM_FILE  _bam, *bam = &_bam;

  uint32 bptr, blen;
  int    last, notfound;

  uint8 *block;
  uint32 bsize, ssize;
  size_t tsize;

#ifdef DEBUG_FIND
  fprintf(stderr,"Searching from %lld\n",fpos);
  fflush(stderr);
#endif

  lseek(fid,fpos,SEEK_SET);
  blen = 0;
  bptr = 0;
  last = 0;

  //  Search until find a gzip block header

  notfound = 1;
  while (notfound)
    { int    j;
      uint32 isize, crc;

      fpos += bptr;      //   Get more data at level of IO blocks
      if (last)
        { fprintf(stderr,"\n%s: Could not find bam block structure!\n",Prog_Name);
          Clean_Exit(1);
        }
      else
        { uint32 x = blen-bptr;
          memmove(buf,buf+bptr,x);
          blen = x + read(fid,buf+x,IO_BLOCK-x);
          if (blen < IO_BLOCK)
            last = 1;
#ifdef DEBUG_FIND
          fprintf(stderr,"Loading %d(last=%d)\n",blen,last);
          fflush(stderr);
#endif
          bptr = 0;
        }

      while (bptr < blen)          //  Search IO block for Gzip block start
        { if (buf[bptr++] != 31)
            continue;
          if ( buf[bptr] != 139)
            continue;
          bptr += 1;
          if ( buf[bptr] != 8)
            continue;
          bptr += 1;
          if ( buf[bptr] != 4)
            continue;
  
#ifdef DEBUG_FIND
          fprintf(stderr,"  Putative header @ %d\n",bptr-3);
          fflush(stderr);
#endif

          if (bptr + 12 > blen)
            { if (last)
                continue;
              bptr -= 3;
              break;
            }
  
          j = bptr+9;
          if (buf[j] != 66)
		  continue;
          j += 1;
          if (buf[j] != 67)
            continue;
          j += 1;
          if (buf[j] != 2)
            continue;
          j += 1;
          if (buf[j] != 0)
            continue;
          j += 1;
    
          bsize = getint(buf+j,2)+1;
          block = buf+(bptr-3);

          if ((bptr-3) + bsize > blen)
            { if (last)
                continue;
              bptr -= 3;
              break;
            }

#ifdef DEBUG_FIND
          fprintf(stderr,"    Putative Extra %d\n",bsize);
          fflush(stderr);
#endif
  
          isize = getint(block+(bsize-4),4);
          crc   = getint(block+(bsize-8),4);
  
          if (libdeflate_gzip_decompress(decomp,block,bsize,bam->bam,BAM_BLOCK,&tsize) != 0)
            continue;
          ssize = tsize;

          if (ssize == isize && crc == libdeflate_crc32(0,bam->bam,ssize))
            { bptr -= 3;
              fpos  += bptr;
              notfound = 0;

#ifdef DEBUG_FIND
              fprintf(stderr,"    First block at %lld (%d)\n",fpos,ssize);
              fflush(stderr);
#endif
	      break;
            }
        }
    }

  //  Have found a gzip/bam block start, now scan blocks until can identify the start
  //    of a sequence entry

  bam->fid      = fid;      //  Kick-start BAM stream object
  bam->last     = last;
  bam->buf      = buf;
  bam->blen     = blen;
  bam->bptr     = bptr+bsize;
  bam->bsize    = bsize;
  bam->ssize    = ssize;
  bam->loc.fpos = fpos;
  bam->loc.boff = 0;
  bam->decomp   = decomp;

  while ( ! bam_eof(bam))
    { int    j, k;
      int    run, out, del;
      int    lname, lcigar, lseq, ldata;

      block = bam->bam;
      ssize = bam->ssize;

      run = HEADER_LEN-1;
      out = 1;
      for (j = HEADER_LEN; j < 10000; j++)
        if (DNA[block[j]])
          { if (out && j >= run+SEQ_RUN)
              {
#ifdef DEBUG_FIND
                fprintf(stderr,"      Possible seq @ %d\n",run+1);
                fflush(stderr);
#endif
                for (k = run-(HEADER_LEN-1); k >= 0; k--)
                  { ldata  = getint(block+k,4);
                    lname  = block[k+12];
                    lcigar = getint(block+(k+16),2);
                    lseq   = getint(block+(k+20),4);
                    if (lname > 0 && lcigar >= 0 && lseq > 0 &&
                        (lseq+1)/2+lseq+lname+(lcigar<<2) < ldata)
                      { del = (k+35+lname+(lcigar<<2)) - run;
                        if (del >= 0 && del < SEQ_RUN/2)
                          { if (valid_name((char *) (block+(k+36)),lname) &&
                                valid_cigar((int32 *) (block+(k+36+lname)),lcigar,lseq))
                              { parm->beg.fpos = bam->loc.fpos;
                                parm->beg.boff = k;
#ifdef DEBUG_FIND
                                fprintf(stderr,"      Found @ %d: '%s':%d\n",k,block+(k+36),lseq);
                                fflush(stderr);
#endif

                                close(fid);

                                return;
                              }
                          }
                      }
                  }
                out = 0;
              }
          }
        else
          { run = j;
            out = 1;
          }

      bam_get(bam,NULL,ssize);
    }

  parm->beg.fpos = -1;
}

  //  Find next identifiable sam entry location forward of parm->fpos in parm->fid
  //    Return location is in parm->beg.  NB: works to skip sam header as well

static void sam_nearest(Thread_Arg *parm)
{ uint8       *buf  = parm->buf;
  int          fid  = parm->fid;

  BAM_FILE  _bam, *bam = &_bam;

  bam->decomp = parm->decomp;
  sam_start(bam,fid,buf,&(parm->beg));

  sam_getline(bam);
  if (parm->beg.fpos == 0)
    { while (buf[bam->bptr] == '@')
        sam_getline(bam);
    }
  parm->beg = bam->loc;
}


/*******************************************************************************************
 *
 *  Routines to scan and parse bam and sam entries
 *
 ********************************************************************************************/

typedef struct
  { int    hlen;
    char  *header;
    uint32 flags;
    int    len;
    char  *seq;
    char  *qvs;
    int    lmax;     //  current max size for seq, arr, and qvs
    int    dmax;     //  current max size for data
    uint8 *data;     //  data buffer
  } samRecord;

static char INT_2_IUPAC[16] = "=acmgrsvtwyhkdbn";

  //  Scan next bam entry and load PacBio info in record 'theR'
 
static int bam_record_scan(BAM_FILE *sf, samRecord *theR)
{ int ldata, lname, lseq, lcigar, aux;

  { uint8  x[36];     //  Process 36 byte header

    bam_get(sf,x,36);

    ldata  = getint(x,4) - 32;
    lname  = getint(x+12,1);
    lcigar = getint(x+16,2);
    lseq   = getint(x+20,4);

    if (ldata < 0 || lseq < 0 || lname < 1)
      { fprintf(stderr,"\n%s: Non-sensical BAM record, file corrupted?\n",Prog_Name);
        Clean_Exit(1);
      }

    aux = lname + ((lseq + 1)>>1) + lseq + (lcigar<<2);
    if (aux > ldata)
      { fprintf(stderr,"\n%s: Non-sensical BAM record, file corrupted?\n",Prog_Name);
        Clean_Exit(1);
      }

    if (lseq > theR->lmax)
      { theR->lmax = 1.2*lseq + 1000;
        theR->seq  = (char *) Realloc(theR->seq,2*theR->lmax,"Reallocating sequence buffer");
        if (theR->seq == NULL)
          Clean_Exit(1);
        theR->qvs  = theR->seq + theR->lmax;
      }

    if (ldata > theR->dmax)
      { theR->dmax = 1.2*ldata + 1000;
        theR->data = (uint8 *) Realloc(theR->data,theR->dmax,"Reallocating data buffer");
        if (theR->data == NULL)
          Clean_Exit(1);
      }

    bam_get(sf,theR->data,ldata);

    if ((getint(x+18,2) & 0x900) != 0)
      { theR->len = 0;
        return (0);
      }

    if (lseq <= 0)
      fprintf(stderr,"%s: WARNING: no sequence for subread !?\n",Prog_Name);

    theR->header = (char *) theR->data;
    theR->hlen   = lname;
    theR->len    = lseq;

    { uint8 *t;
      char  *s;
      int    i, e;

      t = theR->data + (lname + (lcigar<<2));
      s = theR->seq;
      lseq -= 1;
      for (e = i = 0; i < lseq; e++)
        { s[i++] = INT_2_IUPAC[t[e] >> 4];
          s[i++] = INT_2_IUPAC[t[e] & 0xf];
        }
      if (i <= lseq)
        s[i] = INT_2_IUPAC[t[e++] >> 4];
      lseq += 1;

      t += e;
      s  = theR->qvs;
      if (t[0] == 0xff)
        return (0);
      for (i = 0; i < lseq; i++)
        s[i] = t[i] + 33;
    }
  }

  return (1);
}

  //  Scan next bam entry and load PacBio info in record 'theR'

static char  IUPAC_2_DNA[256] =
  { 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',   'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
    'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',   'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
    'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',   'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
    'a', 'c', 'g', 't', 'a', 'a', 'a', 'a',   'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',

    'a', 'a', 'c', 'c', 'a', 'a', 'a', 'g',   'a', 'a', 'a', 'g', 'a', 'a', 'a', 'a',
    'a', 'a', 'a', 'c', 't', 'a', 'a', 'a',   'a', 'c', 'a', 'a', 'a', 'a', 'a', 'a',
    'a', 'a', 'c', 'c', 'a', 'a', 'a', 'g',   'a', 'a', 'a', 'g', 'a', 'a', 'a', 'a',
    'a', 'a', 'a', 'c', 't', 'a', 'a', 'a',   'a', 'c', 'a', 'a', 'a', 'a', 'a', 'a',
  };

#define CHECK(cond, msg)                                \
{ if ((cond))                                           \
    { fprintf(stderr,"\n%s: %s\n",Prog_Name, msg);      \
       Clean_Exit(1); 	                                \
    }                                                   \
}

#define NEXT_ITEM(b,e)                                  \
{ b = e;                                                \
  while (*e != '\0' && *e != '\t')                      \
    e++;                                                \
  CHECK( *e == '\0', "Missing one or more fields")      \
  *e = 0;                                               \
}

static int sam_record_scan(BAM_FILE *sf, samRecord *theR)
{ char      *p;
  int        qlen, flags;

  //  read next line

  theR->data = sam_getline(sf);

  p = (char *) theR->data;

  { char *q, *seq;     //  Load header and sequence from required fields
    int   i;

    NEXT_ITEM(q,p)
    qlen = p-q;
    CHECK( qlen <= 1, "Empty header name")
    CHECK( qlen > 255, "Header is too long")

    theR->header = q;

    flags = strtol(q=p+1,&p,0);
    CHECK( p == q, "Cannot parse flags")

    for (i = 0; i < 7; i++)   // Skip next 8 required fields
      { p = index(p+1,'\t');
        CHECK( p == NULL, "Too few required fields in SAM record, file corrupted?")
      }
    p += 1;

    NEXT_ITEM(q,p)
    qlen = p-q;
    CHECK (*q == '*', "No sequence for read?");

    if (qlen > theR->lmax)
      { theR->lmax = 1.2*qlen + 1000;
        theR->seq  = (char *) Realloc(theR->seq,2*theR->lmax,"Reallocating sequence buffer");
        if (theR->seq == NULL)
          Clean_Exit(1);
        theR->qvs  = theR->seq + theR->lmax;
      }

    if ((flags & 0x900) != 0)
      { theR->len = 0;
        return (0);
      }

    if (qlen <= 0)
      fprintf(stderr,"%s: WARNING: no sequence for subread !?\n",Prog_Name);

    theR->len = qlen;
    seq = theR->seq;
    for (i = 0; i < qlen; i++)
      seq[i] = IUPAC_2_DNA[(int) (*q++)];

    q = ++p;
    p = index(p,'\t');
    CHECK( p == NULL, "No auxilliary tags in SAM record, file corrupted?")

    if (*q == '*')
      return (0);
    qlen = p-q;
    seq = theR->qvs;
    for (i = 0; i < qlen; i++)
      seq[i] = *q++;
  }

  return (1);
}

/*******************************************************************************************
 *
 *  Parallel:  Each thread processes a contiguous stripe across the input files
 *               sending the compressed binary data lines to their assigned OneFile.
 *
 ********************************************************************************************/

  //  Write subread data in samRecord rec to non-NULL file types

static void *bam_output_thread(void *arg)
{ Thread_Arg  *parm  = (Thread_Arg *) arg;
  File_Object *fobj  = parm->fobj;
  uint8       *buf   = parm->buf;
  int          action = parm->action;
  DATA_BLOCK  *dset   = &parm->block;
  int          tid    = parm->thread_id;

  samRecord    _theR, *theR = &_theR;
  BAM_FILE     _bam, *bam = &_bam;

  int64        epos;
  uint32       eoff;
  int          isbam;
  int          f, fid;
  int64        totread, fbeg;

  int64 estbps, cumbps, nxtbps, pct1;
  int   CLOCK;

  estbps = nxtbps = pct1 = cumbps = 0;
  if (VERBOSE && tid == 0 && action != SAMPLE)
    { estbps = dset->ratio / ITHREADS;
      nxtbps = pct1 = estbps/100;
      fprintf(stderr,"\n    0%%");
      fflush(stderr);
      CLOCK = 1;
    }
  else
    CLOCK = 0;

  //  Know the max size of sequence and data from pass 1, so set up accordingly

  theR->data = NULL;
  if (fobj->ftype == BAM)
    { theR->dmax = 50000;
      theR->data = Malloc(theR->dmax,"Allocating sequence array");
      if (theR->data == NULL)
        Clean_Exit(1);
    }
  theR->lmax = 75000;
  theR->seq  = Malloc(2*theR->lmax,"Allocating sequence array");
  if (theR->seq == NULL)
    Clean_Exit(1);
  theR->qvs = theR->seq + theR->lmax;

  bam->decomp = parm->decomp;

  totread = 0;

  for (f = parm->bidx; f <= parm->eidx; f++)
    { fid   = open(fobj[f].path,O_RDONLY);
      isbam = (fobj[f].ftype == BAM);
      if (f < parm->eidx || parm->end.fpos >= fobj[f].fsize)
        { epos  = fobj[f].fsize;
          eoff  = 0;
        }
      else 
        { epos  = parm->end.fpos;
          eoff  = parm->end.boff;
        }
      if (f > parm->bidx || parm->beg.fpos == 0)
        { parm->beg.fpos = 0;
          parm->fid      = fid;
          if (isbam)
            skip_bam_header(parm);
          else
            sam_nearest(parm);
        }

      if (isbam)
        bam_start(bam,fid,buf,&(parm->beg));
      else
        sam_start(bam,fid,buf,&(parm->beg));

      fbeg = lseek(fid,0,SEEK_CUR);

#ifdef DEBUG_IO
      printf("Block: %12lld / %5d to %12lld / %5d --> %8lld\n",bam->loc.fpos,bam->loc.boff,
                                                               epos,eoff,epos - bam->loc.fpos);
      fflush(stdout);
#endif

      while (1)
        { if (bam->loc.fpos >= epos && bam->loc.boff >= eoff)
            break;

          if (isbam)
            bam_record_scan(bam,theR);
          else
            sam_record_scan(bam,theR);

          if (theR->len <= 0)
            continue;

#ifdef DEBUG_BAM_RECORD
          fprintf(stderr,"S = '%s'\n",theR->seq);
          if (hasQV)
            fprintf(stderr,"Q = '%.*s'\n",theR->len,theR->qvs);
#endif

          while (Add_Data_Block(dset,theR->len,theR->seq))
            { if (action == SAMPLE)
                if (isbam)
                  { int unused = (bam->blen - (bam->bptr + bam->bsize))
                               - bam->loc.boff * ((1.*bam->bsize) / bam->ssize);
                    dset->ratio = (1.*parm->work)
                                / ((totread+lseek(fid,0,SEEK_CUR))-(unused+fbeg+dset->rem));
                    close(fid);
                    return (NULL);
                  }
                else
                  { dset->ratio = (1.*parm->work) / (totread+bam->loc.fpos-dset->rem);
                    close(fid);
                    return (NULL);
                  }
              else
                { CALL_BACK(dset,tid);
                  if (CLOCK)
                    { cumbps += dset->totlen;
                      if (cumbps >= nxtbps)
                        { fprintf(stderr,"\r  %3d%%",(int) ((100.*cumbps)/estbps));
                          fflush(stderr);
                          nxtbps = cumbps+pct1;
                        }
                    }
                  Reset_Data_Block(dset,0);
                }
            }
        }

      totread += lseek(fid,0,SEEK_CUR) - fbeg;
      close(fid);
    }

  if (action == SAMPLE)
    { dset->ratio = (1.*parm->work) / totread;
      return (NULL);
    }

  if (dset->nreads > 0)
    CALL_BACK(dset,tid);

  if (CLOCK)
    fprintf(stderr,"\r         \r");

  free(theR->seq);
  if (fobj->ftype == BAM)
    free(theR->data);

  return (NULL);
}


/*******************************************************************************************
 *
 *  CRAM SPECIFIC CODE
 *
 ********************************************************************************************/

/*******************************************************************************************
 *
 *  cram_nearest: Routine to find next entry given arbitray start point.
 *
 ********************************************************************************************/

#include "htslib/hts.h"
#include "htslib/hfile.h"
#include "cram/cram.h"

int ITF_LEN[16] = { 0, 0, 0, 0,
                    0, 0, 0, 0,
                    1, 1, 1, 1,
                    2, 2, 3, 4 };

static void scan_itf8(hFILE *fp)
{ char buf[4];
  int  nb = ITF_LEN[hgetc(fp)>>4];
  if (nb > 0)
    hread(fp,buf,nb);
}

const int LTF_LEN[256] = {
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,

    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,

    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,

    2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,
    2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,
    3, 3, 3, 3,  3, 3, 3, 3,  3, 3, 3, 3,  3, 3, 3, 3,
    4, 4, 4, 4,  4, 4, 4, 4,  5, 5, 5, 5,  6, 6, 7, 8 };

static void scan_ltf8(hFILE *fp)
{ char buf[8];
  int  nb = LTF_LEN[hgetc(fp)];
  if (nb > 0)
    hread(fp,buf,nb);
}

static int int32_decode(hFILE *fp, int32 *val)
{ char *v = (char *) val;

  if (hread(fp,v,4) != 4)
    return (1);
#if __ORDER_LITTLE_ENDIAN__ != __BYTE_ORDER__
  { char  x;
    x = v[0];
    v[0] = v[3];
    v[3] = x;
    x = v[1];
    v[1] = v[2];
    v[2] = x;
  }
#endif
  return (0);
}

// reading cram block, header is a block so wrapped.

static int64 scan_container(cram_fd *fd)
{ int     i, len;
  int32   nslice;
  hFILE  *fp = fd->fp;

  if (int32_decode(fp,&len))                    //  total length
    return (-1);
  scan_itf8(fp);                                //  ref seq id
  scan_itf8(fp);                                //  start pos on ref
  scan_itf8(fp);                                //  align span
  scan_itf8(fp);                                //  # of records
  if (CRAM_MAJOR_VERS(fd->version) > 1)
    { if (CRAM_MAJOR_VERS(fd->version) >= 3)    //  record counter
        scan_ltf8(fp);
      else
        scan_itf8(fp);
      scan_ltf8(fp);                            //  # bps
    }
  scan_itf8(fp);                                //  # of blocks
  itf8_decode(fd,&nslice);                      //  # of slices
  for (i = 0; i < nslice; i++)                  //  landmarks
    scan_itf8(fp);
  hread(fp,&nslice,4);                          //  crc code
  hseek(fp,len,SEEK_CUR);                       //  skip contents of container
  return (htell(fp));
}

static int64 *genes_cram_index(char *path, int64 fsize, int64 *zsize)
{ cram_fd *fd;
  int64    cash[10], *zoff;
  int64    s, e;
  int      i, j;
 
  fd = cram_open(path,"r");

  cash[0] = htell(fd->fp);
  for (i = 1; i < 10; i++)
    { cash[i] = s = scan_container(fd);
      if (s == fsize)
        { zoff = Malloc(sizeof(int64)*i,"Allocating cram index"); 
          if (zoff == NULL)
            Clean_Exit(1);
          for (j = 0; j < i; i++)
            zoff[j] = cash[i];
          *zsize = i-1;
          return (zoff);
        }
    }

  e = ((fsize-(cash[0]+38))/(cash[9]-cash[0])+1)*9 + 100;
  zoff = Malloc(sizeof(int64)*e,"Allocating cram index"); 
  if (zoff == NULL)
    Clean_Exit(1);
  for (j = 0; j < i; j++)
    zoff[j] = cash[j];

  while (s != fsize)
    { zoff[i++] = s = scan_container(fd);
      if (i >= e)
        { e = ((fsize-(zoff[0]+38.))/(zoff[i-1]-zoff[0]))*(i-1) + 100;
          zoff = Realloc(zoff,sizeof(int64)*e,"Allocating cram index"); 
          if (zoff == NULL)
            Clean_Exit(1);
        }
    }

  cram_close(fd);

  *zsize = i-2;
  return (zoff);
}       

static void cram_nearest(Thread_Arg *data)
{ File_Object *inp   = data->fobj + data->bidx;
  int64       *zoffs = inp->zoffs;
  int64        zsize = inp->zsize;
  int64        pos   = data->beg.fpos;
  int i;

  for (i = 0; i < zsize; i++)
    if (zoffs[i] >= pos)
      break;
  if (i == zsize)
    data->beg.fpos = -1;
  else
    data->beg.fpos = zoffs[i];
}

static void *cram_output_thread(void *arg)
{ Thread_Arg  *parm   = (Thread_Arg *) arg;
  File_Object *fobj   = parm->fobj;
  int          action = parm->action;
  DATA_BLOCK  *dset   = &parm->block;
  int          tid    = parm->thread_id;

  File_Object *inp;
  int          f;
  cram_fd     *fid;
  int64        bpos, epos;
  int64        totread;
  char        *line;
  int          o, omax;

  int64 estbps, cumbps, nxtbps, pct1;
  int   CLOCK;

  cumbps = nxtbps = estbps = pct1 = 0;
  if (VERBOSE && tid == 0 && action != SAMPLE)
    { estbps = dset->ratio / ITHREADS;
      nxtbps = pct1 = estbps/100;
      fprintf(stderr,"\n    0%%");
      fflush(stderr);
      CLOCK = 1;
    }
  else
    CLOCK = 0;

  totread = 0;
  omax    = dset->maxbps;
  line    = dset->bases;

  for (f = parm->bidx; f <= parm->eidx; f++)
    { inp = fobj+f;
      fid = cram_open(inp->path,"r");
      if (f < parm->eidx || parm->end.fpos >= inp->zoffs[inp->zsize])
        epos  = inp->zoffs[inp->zsize];
      else
        epos  = parm->end.fpos;
      if (f > parm->bidx || parm->beg.fpos < inp->zoffs[0])
        bpos  = inp->zoffs[0];
      else
        bpos  = parm->beg.fpos;
      hseek(fid->fp,bpos,SEEK_SET);

#ifdef DEBUG_IO
      fprintf(stderr,"Block: %12lld to %12lld --> %8lld\n",bpos,epos,epos - bpos);
      fflush(stderr);
#endif

      o = 0;
      while (1)
        { cram_record *rec;
          char        *seq;
          int          len, ovl;

          rec = cram_get_seq(fid);
          if (rec == NULL)
            break;

          if (htell(fid->fp) > epos)
            break;

          seq = (char *) (rec->s->seqs_blk->data+rec->seq);
          if (COMPRESS)
            len = homo_compress(seq,rec->len);
          else
            len = rec->len;

          while (o+len > omax)
            { ovl = omax-o;
              memcpy(line+o,seq,ovl);
              line[omax] = '\0';
              dset->boff[++dset->nreads] = omax+1;
              dset->rem = 1;
              DUMP(1.,(bpos-htell(fid->fp))+(len-ovl),cram_close)
              dset->rem = 0;
              Reset_Data_Block(dset,1);
              o = KMER-1;
              len -= ovl;
              seq += ovl; 
            }
          memcpy(line+o,seq,len);
          o += len;
          line[o++] = '\0';
          dset->boff[++dset->nreads] = o;
          if (o > omax-DT_MINIM || dset->nreads >= dset->maxrds)
            { DUMP(1.,bpos-htell(fid->fp),cram_close)
              Reset_Data_Block(dset,0);
              o = 0;
            }
        }

      totread += epos-bpos;
      cram_close(fid);
    }

  dset->totlen = dset->boff[dset->nreads] - dset->nreads;
  if (action == SAMPLE)
    { dset->ratio = (1.*parm->work) / totread;
      return (NULL);
    }
  else if (dset->nreads > 0)
    CALL_BACK(dset,tid);

  if (CLOCK)
    fprintf(stderr,"\r         \r");

  return (NULL);
}



/*******************************************************************************************
 *
 *  DAZZLER SPECIFIC CODE
 *
 ********************************************************************************************/

#define DB_QV   0x03ff   //  Mask for 3-digit quality value
#define DB_CCS  0x0400   //  This is the second or later of a group of subreads from a given insert
#define DB_BEST 0x0800   //  This is the "best" subread of a given insert (may be the only 1)

#define DB_ARROW 0x2     //  DB is an arrow DB
#define DB_ALL   0x1     //  all wells are in the trimmed DB

typedef struct
  { int     origin; //  Well # (DB), Contig # (DAM)
    int     rlen;   //  Length of the sequence (Last pulse = fpulse + rlen)
    int     fpulse; //  First pulse (DB), left index of contig in scaffold (DAM)
    int64   boff;   //  Offset (in bytes) of compressed read in 'bases' file, or offset of
                    //    uncompressed bases in memory block
    int64   coff;   //  Offset (in bytes) of compressed quiva streams in '.qvs' file (DB),
                    //  Offset (in bytes) of scaffold header string in '.hdr' file (DAM)
                    //  4 compressed shorts containing snr info if an arrow DB.
    int     flags;  //  QV of read + flags above (DB only)
  } DAZZ_READ;

typedef struct
  { int         ureads;     //  Total number of reads in untrimmed DB
    int         treads;     //  Total number of reads in trimmed DB
    int         cutoff;     //  Minimum read length in block (-1 if not yet set)
    int         allarr;     //  DB_ALL | DB_ARROW
    float       freq[4];    //  frequency of A, C, G, T, respectively

    //  Set with respect to "active" part of DB (all vs block, untrimmed vs trimmed)

    int         maxlen;     //  length of maximum read (initially over all DB)
    int64       totlen;     //  total # of bases (initially over all DB)

    int         nreads;     //  # of reads in actively loaded portion of DB
    int         trimmed;    //  DB has been trimmed by cutoff/all
    int         part;       //  DB block (if > 0), total DB (if == 0)
    int         ufirst;     //  Index of first read in block (without trimming)
    int         tfirst;     //  Index of first read in block (with trimming)

       //  In order to avoid forcing users to have to rebuild all thier DBs to accommodate
       //    the addition of fields for the size of the actively loaded trimmed and untrimmed
       //    blocks, an additional read record is allocated in "reads" when a DB is loaded into
       //    memory (reads[-1]) and the two desired fields are crammed into the first two
       //    integer spaces of the record.

    char       *path;       //  Root name of DB for .bps, .qvs, and tracks
    int         loaded;     //  Are reads loaded in memory?
    void       *bases;      //  file pointer for bases file (to fetch reads from),
                            //    or memory pointer to uncompressed block of all sequences.
    DAZZ_READ  *reads;      //  Array [-1..nreads] of DAZZ_READ
    void       *tracks;     //  Linked list of loaded tracks
  } DAZZ_DB;

static void read_DB_stub(char *path, int *cut, int *all)
{ FILE *dbfile;
  char  buf1[10100];
  char  buf2[10100];
  int   nread;

  int   i;
  int   nfiles;
  int   nblocks;
  int64 size;

  dbfile = fopen(path,"r");
  if (dbfile == NULL)
    { fprintf(stderr,"\n%s: Cannot open stub file %s\n",Prog_Name,path);
      Clean_Exit(1);
    }

  if (fscanf(dbfile,"files = %9d\n",&nfiles) != 1)
    goto stub_trash;

  for (i = 0; i < nfiles; i++)
    if (fscanf(dbfile,"  %9d %s %s\n",&nread,buf1,buf2) != 3)
      goto stub_trash;

  if (fscanf(dbfile,"blocks = %9d\n",&nblocks) != 1)
    { fprintf(stderr,"\n%s: DB %s has not been split or its stub file is junk\n",Prog_Name,path);
      fclose(dbfile);
      Clean_Exit(1);
    }

  if (fscanf(dbfile,"size = %11lld cutoff = %9d all = %1d\n",&size,cut,all) != 3)
    goto stub_trash;

  fclose(dbfile);
  return;

stub_trash:
  fprintf(stderr,"\n%s: Stub file %s is junk\n",Prog_Name,path);
  fclose(dbfile);
  Clean_Exit(1);
}

static int64 *get_dazz_offsets(FILE *idx, int64 *zsize)
{ DAZZ_DB   header;
  DAZZ_READ duplo[2];
  int64    *index;
  int       i, p, nreads;

  fread(&header,sizeof(DAZZ_DB),1,idx);
  nreads = header.ureads;
  index = (int64 *) Malloc(sizeof(int64)*(nreads+3)/2,"Allocating DB index");
  p = 0;
  for (i = 0; i < nreads; i += 2)
    { fread(duplo,sizeof(DAZZ_READ),2,idx);
      index[p++] = duplo[0].boff;
    }
  *zsize = p;
  return (index); 
}

static int64 get_dazz_lengths(FILE *idx, int64 *index, int DB_cut, int DB_all)
{ DAZZ_DB   header;
  DAZZ_READ read;
  int      *in = (int *) index;
  int       cutoff, allflag;
  int       i, nreads;

  cutoff = DB_cut;
  if (DB_all)
    allflag = 0;
  else
    allflag = DB_BEST;

  fread(&header,sizeof(DAZZ_DB),1,idx);
  nreads = header.ureads;
  for (i = 0; i < nreads; i++)
    { fread(&read,sizeof(DAZZ_READ),1,idx);
      if (read.rlen >= cutoff && (read.flags & DB_BEST) >= allflag)
        in[i] = read.rlen;
      else
        in[i] = -read.rlen;
    }
  return (nreads); 
}

static void dazz_nearest(Thread_Arg *data)
{ File_Object *inp   = data->fobj + data->bidx;
  int64       *zoffs = inp->zoffs;
  int64        zsize = inp->zsize;
  int64        pos   = data->beg.fpos;
  int i;

  i = ((1.*pos) / inp->fsize) * zsize;
  if (zoffs[i] >= data->beg.fpos)
    { while (i >= 0)
        { if (zoffs[i] < data->beg.fpos)
            break;
          i -= 1;
        }
      i += 1;
    }
  else
    while (i < zsize)
      { if (zoffs[i] >= data->beg.fpos)
          break;
        i += 1;
      }
  if (i == zsize)
    data->beg.fpos = -1;
  else
    { data->beg.fpos = zoffs[i];
      data->beg.boff = 2*i;
    }
}

static void uncompress_read(int len, char *s)
{ static char letter[4] = { 'a', 'c', 'g', 't' };
  int   i, tlen, byte;
  char *s0, *s1, *s2, *s3;
  char *t;

  s0 = s;
  s1 = s0+1;
  s2 = s1+1;
  s3 = s2+1;

  tlen = (len-1)/4;

  t = s+tlen;
  for (i = tlen*4; i >= 0; i -= 4)
    { byte = *t--;
      s0[i] = letter[(byte >> 6) & 0x3];
      s1[i] = letter[(byte >> 4) & 0x3];
      s2[i] = letter[(byte >> 2) & 0x3];
      s3[i] = letter[byte & 0x3];
    }
  s[len] = '\0';
}

static void *dazz_output_thread(void *arg)
{ Thread_Arg  *parm   = (Thread_Arg *) arg;
  File_Object *fobj   = parm->fobj;
  int          action = parm->action;
  DATA_BLOCK  *dset   = &parm->block;
  int          tid    = parm->thread_id;

  File_Object *inp;
  int          f;
  FILE        *fid;
  int64        bpos, epos;
  int64        totread;

  int  *zoffs;
  int   r, len;
  int   covl, ovl, new;
  int   o, omax;
  int  *boff;
  char *line;

  int64 estbps, cumbps, nxtbps, pct1;
  int   CLOCK;

  nxtbps = pct1 = estbps = cumbps = 0;
  if (VERBOSE && tid == 0 && action != SAMPLE)
    { estbps = dset->ratio / ITHREADS;
      nxtbps = pct1 = estbps/100;
      fprintf(stderr,"\n    0%%");
      fflush(stderr);
      CLOCK = 1;
    }
  else
    CLOCK = 0;

  totread = 0;
  omax    = dset->maxbps;
  line    = dset->bases;
  boff    = dset->boff;

  for (f = parm->bidx; f <= parm->eidx; f++)
    { inp = fobj+f;
      fid = fopen(inp->path,"r");
      if (f < parm->eidx)
        epos  = inp->zsize;
      else
        epos  = parm->end.boff;
      if (f > parm->bidx)
        { bpos  = 0;
          r     = 0;
	}
      else
        { bpos  = parm->beg.fpos;
          r     = parm->beg.boff;
        }
      fseek(fid,bpos,SEEK_SET);
      zoffs = (int *) inp->zoffs;

#ifdef DEBUG_IO
      fprintf(stderr,"Block: %12lld: %d to %lld\n",bpos,r,epos);
      fflush(stderr);
#endif

      o = 0;
      while (r < epos)
        { len = zoffs[r++];
          if (len < 0)
            { fseeko(fid,(3-len)>>2,SEEK_CUR);
              continue;
            }

          covl = 0;
          while (o+len > omax)
            { int x = 0;

              ovl = ((omax-o) >> 2) << 2;
              fread(line+o,ovl>>2,1,fid);
              uncompress_read(ovl,line+o);
              len -= ovl;
              if (COMPRESS)
                { if (covl > 0)
                    new = homo_compress(line+(o-1),ovl+1) - 1;
                  else
                    new = homo_compress(line+o,ovl);
                  o += new;
                  covl += new;
                  if (covl <= KMER)
                    continue;
                  x = line[--o];
                }
              else
                o += ovl;
              line[o++] = '\0';
              boff[++dset->nreads] = o;
              dset->rem = 1;
              DUMP(1.,(bpos-ftello(fid))+len,fclose)
              dset->rem = 0;
              Reset_Data_Block(dset,1);
              o = KMER-1;
              if (COMPRESS)
                { line[o++] = x;
                  covl = KMER;
                }
            }
          fread(line+o,(len+3)>>2,1,fid);
          uncompress_read(len,line+o);
          if (COMPRESS)
            { if (covl > 0)
                len = homo_compress(line+(o-1),len+1) - 1;
              else
                len = homo_compress(line+o,len);
            }
          o += len;
          line[o++] = '\0';
          boff[++dset->nreads] = o;
          if (o > omax-DT_MINIM || dset->nreads >= dset->maxrds)
            { DUMP(1.,bpos-ftello(fid),fclose)
              Reset_Data_Block(dset,0);
              o = 0;
            }
        }

      totread += ftello(fid)-bpos;
      fclose(fid);
    }

  dset->totlen = dset->boff[dset->nreads] - dset->nreads;
  if (action == SAMPLE)
    { dset->ratio = (1.*parm->work) / totread;
      return (NULL);
    }
  else if (dset->nreads > 0)
    CALL_BACK(dset,tid);

  if (CLOCK)
    fprintf(stderr,"\r         \r");

  return (NULL);
}


/****************************************************************************************
 *
 *  The top-level interface
 *
 ****************************************************************************************/

  //  Find the ITHREADS partition points in files in argv[0..argc) and return in an
  //    Input_Partition data structure

Input_Partition *Partition_Input(int argc, char *argv[])
{ int         nfiles;
  int         ftype;
  int         need_decon;
  int         need_buf;
  void *    (*output_thread)(void *);
  void      (*scan_header)(Thread_Arg *);
  void      (*find_nearest)(Thread_Arg *);

  File_Object *fobj;
  Thread_Arg  *parm;

  nfiles = argc-1;

  parm = (Thread_Arg *) Malloc(sizeof(Thread_Arg)*NTHREADS,"Allocating input threads");
  fobj = (File_Object *) Malloc (sizeof(File_Object)*nfiles,"Allocating file records"); 
  if (parm == NULL || fobj == NULL)
    Clean_Exit(1);

  //  Find partition points dividing data in all files into NTHREADS roughly equal parts
  //    and then in parallel threads produce the output for each part.

  { int    f, i, t;
    int    gz_ok;
    int64  b, w;
    int64  work;
    uint8 *bf;
    int    file_split;

    //  Get name and size of each file in 'fobj[]', determine type, etc.

    file_split = 0;
    need_decon = 0;
    need_buf   = 0;
    ftype      = FASTA;
    output_thread = fast_output_thread;
    scan_header   = do_nothing;
    find_nearest  = fast_nearest;

    gz_ok = (nfiles >= 3 || nfiles >= NTHREADS/2);
    work = 0;
    for (f = 0; f < nfiles; f++)
      { Fetch_File(argv[f+1],fobj+f,gz_ok);

        if (f > 0)
          { if (fobj[f].ftype != ftype)
              { fprintf(stderr,"\n%s: All files must be of the same type\n",Prog_Name);
                if (gz_ok)
                  Free_Gunzips(fobj,f);
                Clean_Exit(1);
              }
          }
        else
          { ftype = fobj[0].ftype;
            if (ftype == FASTA || ftype == FASTQ)
              { output_thread = fast_output_thread;
                  scan_header   = do_nothing;
                find_nearest  = fast_nearest;
              }
            else if (ftype == BAM)
              { output_thread = bam_output_thread;
                scan_header   = skip_bam_header;
                find_nearest  = bam_nearest;
              }
            else if (ftype == SAM)
              { output_thread = bam_output_thread;
                scan_header   = sam_nearest;
                find_nearest  = sam_nearest;
              }
            else if (ftype == CRAM)
              { output_thread = cram_output_thread;
                scan_header   = do_nothing;
                find_nearest  = cram_nearest;
              }
            else // ftype == DAM or DB
              { output_thread = dazz_output_thread;
                scan_header   = do_nothing;
                find_nearest  = dazz_nearest;
              }
          }
        if ( ! (ftype == CRAM || ftype == DAZZ))
          need_buf = 1;
        if (ftype == BAM)
          need_decon = 1;
        if (fobj[f].zipd)
          file_split = 1;

        work += fobj[f].fsize;
      }
    parm[0].work = work;

    //  If gzip'd files then divide whole files

    if (file_split)
      { if (nfiles <= 1.5*NTHREADS)
          ITHREADS = nfiles;
        else
          ITHREADS = NTHREADS;

        if (VERBOSE)
          { fprintf(stderr,"\nUsing %d threads to read %d %s files some or all",
                           ITHREADS,nfiles,Tstring[fobj->ftype]);
            fprintf(stderr," of which are compressed\n");
            fflush(stderr);
          }

        //  Allocate IO buffer space

        if (need_buf)
          { bf = Malloc(ITHREADS*IO_BLOCK,"Allocating IO_Buffer\n");
            if (bf == NULL)
              { if (gz_ok)
                  Free_Gunzips(fobj,nfiles);
                Clean_Exit(1);
              }
          }
        else
          bf = NULL;

        for (t = 0; t < ITHREADS; t++)
          { parm[t].fobj          = fobj;
            parm[t].output_thread = output_thread;
            parm[t].thread_id     = t;
            parm[t].action        = SPLIT;
            parm[t].bidx          = (t*nfiles)/ITHREADS;
            parm[t].beg.fpos      = 0;
            parm[t].beg.boff      = 0;
            parm[t].eidx          = ((t+1)*nfiles)/ITHREADS-1;
            parm[t].end.fpos      = fobj[parm[t].eidx].fsize;
            parm[t].end.boff      = 0;
            parm[t].buf           = bf + t*IO_BLOCK;
            parm[t].decomp        = NULL;
          }
      }

    //  else split files at achieve load balancing

    else
      {

    //  Allocate work evenly amongst threads, setting up search start
    //    point for each thread.  Also find the beginning of data in
    //    each file that a thread will start in (place in end.fpos)

    if (work/NTHREADS < .02*IO_BLOCK)
      { ITHREADS = work/(.02*IO_BLOCK);
        if (ITHREADS <= 0)
          ITHREADS = 1;
      }
    else
      ITHREADS = NTHREADS;

    //  Allocate IO buffer space

    if (need_buf)
      { bf = Malloc(ITHREADS*IO_BLOCK,"Allocating IO_Buffer\n");
        if (bf == NULL)
          { if (gz_ok)
              Free_Gunzips(fobj,nfiles);
            Clean_Exit(1);
          }
      }
    else
      bf = NULL;

    if (VERBOSE)
      { if (nfiles > 1)
          fprintf(stderr,"\nPartitioning %d %s.%s files into %d parts\n",
                         nfiles,fobj->zipd?"compressed ":"",Tstring[fobj->ftype],NTHREADS);
        else
          fprintf(stderr,"\nPartitioning %d %s.%s file into %d parts\n",
                         nfiles,fobj->zipd?"compressed ":"",Tstring[fobj->ftype],NTHREADS);
        fflush(stderr);
      }

    f = 0;
    t = -1;
    w = fobj[f].fsize;
    for (i = 0; i < ITHREADS; i++)
      { while (w < (i*work)/ITHREADS - .01*IO_BLOCK)
          { f += 1;
            w += fobj[f].fsize;
          }
        b = (i*work)/ITHREADS - (w-fobj[f].fsize); 
        if (b < 0)
          b = 0;

        if (t >= 0 && (f < parm[t].bidx || (f == parm[t].bidx && b <= parm[t].beg.fpos)))
          continue;
        t += 1;

        parm[t].fobj          = fobj;
        parm[t].output_thread = output_thread;
        parm[t].thread_id     = t;
        parm[t].action        = SPLIT;
        if (need_buf)
          parm[t].buf = bf + t*IO_BLOCK;
        else
          parm[t].buf = NULL;
        if (need_decon)
          parm[t].decomp = libdeflate_alloc_decompressor();
        else
          parm[t].decomp = NULL;

        if (b != 0)
          { parm[t].fid = open(fobj[f].path,O_RDONLY);

            parm[t].beg.fpos = 0;
            parm[t].beg.boff = 0;
            scan_header(parm+t);
            parm[t].end = parm[t].beg;

            parm[t].beg.fpos = b;
            parm[t].bidx = f;

            find_nearest(parm+t);
            close(parm[t].fid);

            if (parm[t].beg.fpos < 0)
              { parm[t].beg.fpos = 0;
                parm[t].bidx += 1;
                if (parm[t].bidx >= nfiles)
                  { t -= 1;
                    break;
                  }
              }
            else if (parm[t].beg.fpos <= parm[t].end.fpos)
              parm[t].beg.fpos = 0;
          }
        else
          { parm[t].bidx = f;
            parm[t].beg.fpos = b;
            parm[t].beg.boff = 0;
          }

#ifdef DEBUG_FIND
        fprintf(stderr," %2d: %1d %10lld (%10lld/%d) -> %d %10lld (%d)\n",
                       t,f,b,parm[t].end.fpos,parm[t].end.boff,
                             parm[t].bidx,parm[t].beg.fpos,parm[t].beg.boff);
        fflush(stderr);
#endif
      }

    t += 1;
    if (t < ITHREADS)
      { if (need_buf)
          { bf = Realloc(bf,t*IO_BLOCK,"Allocating IO_Buffer\n");
            if (bf == NULL)
              { if (gz_ok)
                  Free_Gunzips(fobj,nfiles);
                Clean_Exit(1);
              }
            for (i = 0; i < t; i++)
              if (need_buf)
                parm[i].buf = bf + i*IO_BLOCK;
          }
        ITHREADS = t;
      }

    //  If cannot use all threads report it

    if (VERBOSE && NTHREADS != ITHREADS)
      { if (work/NTHREADS < .02*IO_BLOCK)
          fprintf(stderr,"  File%s so small will use only %d thread%s\n",
                         nfiles>1?"s are":" is",ITHREADS,ITHREADS>1?"s ":" ");
        else
          fprintf(stderr,"  File%scould only be divided to use %d thread%s\n",
                         nfiles>1?"s ":" ",ITHREADS,ITHREADS>1?"s ":" ");
      }

    //  Develop end points of each threads work using the start point of the next thread

    for (i = 1; i < ITHREADS; i++)
      if (parm[i].beg.fpos == 0)
        { parm[i-1].end.fpos = fobj[parm[i].bidx-1].fsize;
          parm[i-1].end.boff = 0;
          parm[i-1].eidx     = parm[i].bidx-1;
        }
      else
        { parm[i-1].end.fpos = parm[i].beg.fpos;
          parm[i-1].end.boff = parm[i].beg.boff;
          parm[i-1].eidx     = parm[i].bidx;
        }
    parm[ITHREADS-1].end.fpos = fobj[nfiles-1].fsize;
    parm[ITHREADS-1].end.boff = 0;
    parm[ITHREADS-1].eidx = nfiles-1;

    if (ftype == DAZZ)
      { int   f;
        FILE *idx;

        for (f = 0; f < nfiles; f++)
          { strcpy(fobj[f].path+(strlen(fobj[f].path)-3),"idx");
            idx = fopen(fobj[f].path,"r");
            fobj[f].zsize = get_dazz_lengths(idx,fobj[f].zoffs,fobj[f].DB_cut,fobj[f].DB_all);
            fclose(idx);
            strcpy(fobj[f].path+(strlen(fobj[f].path)-3),"bps");
          }
        parm[ITHREADS-1].end.boff = fobj[nfiles-1].zsize;
        parm[0].beg.boff = 0;
      }

    }  //  end load balancing code
  }

#if defined(DEBUG_FIND) || defined(DEBUG_IO)
  { int i;

    fprintf(stderr,"\nPartition:\n");
    for (i = 0; i < ITHREADS; i++)
      { fprintf(stderr," %2d: %2d / %12lld / %5d",
                       i,parm[i].bidx,parm[i].beg.fpos,parm[i].beg.boff);
        fprintf(stderr,"  -  %2d / %12lld / %5d\n",
                         parm[i].eidx,parm[i].end.fpos,parm[i].end.boff);
      }
    fflush(stderr);
  }
#endif

  return ((Input_Partition *) parm);
}

  //  Get a single block from the start of the input

static Thread_Arg cust;

DATA_BLOCK *Get_First_Block(Input_Partition *parts, int64 numbp)
{ Thread_Arg *parm = (Thread_Arg *) parts;

  cust = parm[0];
  cust.action = SAMPLE;

  cust.block.maxrds = numbp/150;
  cust.block.maxbps = numbp + cust.block.maxrds;
  cust.block.bases  = Malloc(sizeof(char)*(cust.block.maxbps+1),"Allocating first data block");
  cust.block.boff   = Malloc(sizeof(int)*(cust.block.maxrds+1),"Allocating first data block");

  cust.eidx   = parm[ITHREADS-1].eidx;
  cust.end    = parm[ITHREADS-1].end;

  Reset_Data_Block(&cust.block,0);
  cust.block.rem = 0;
  cust.output_thread(&cust);

#ifdef DEBUG_TRAIN
  Print_Block(&cust.block,0);
  Clean_Exit(1);
#endif

  return (&cust.block);
}

  //  Free the first block returned by Get_First_Block

void Free_First_Block(DATA_BLOCK *block)
{ free(block->bases);
  free(block->boff);
}

  //  Return root or pwd of the name of the first file

char *First_Root(Input_Partition *io)
{ Thread_Arg *parm = (Thread_Arg *) io;
  char       *root;

  root = Strdup(parm[0].fobj[0].root,"Allocating root name for FastK");
  return (root);
}

char *First_Pwd(Input_Partition *io)
{ Thread_Arg *parm = (Thread_Arg *) io;
  char       *pwd;

  pwd = Strdup(parm[0].fobj[0].pwd,"Allocating pwd for FastK");
  return (pwd);
}

   //  Distribute k-mers

void Scan_All_Input(Input_Partition *parts)
{ Thread_Arg *parm = (Thread_Arg *) parts;
#if !defined(DEBUG_IO) && !defined(DEBUG_OUT)
  pthread_t   threads[ITHREADS];
#endif
  char  *bases;
  int   *boff;
  int    i;

  parm[0].block.ratio = cust.block.ratio * cust.block.totlen;

  bases = Malloc(sizeof(char)*(DT_BLOCK+3)*ITHREADS,"Allocating data blocks");
  boff  = Malloc(sizeof(int)*(DT_READS+1)*ITHREADS,"Allocating data blocks");
  for (i = 0; i < ITHREADS; i++)
    { parm[i].block.bases  = bases + (DT_BLOCK+3)*i;
      parm[i].block.boff   = boff  + (DT_READS+1)*i;
      parm[i].block.maxbps = DT_BLOCK;
      parm[i].block.maxrds = DT_READS;
      Reset_Data_Block(&parm[i].block,0);
      parm[i].block.rem    = 0;
    }

#if defined(DEBUG_IO) || defined(DEBUG_OUT)
  for (i = 0; i < ITHREADS; i++)
    { fprintf(stderr,"Thread %d\n",i);
      parm[0].output_thread(parm+i);
    }
#else
  for (i = 0; i < ITHREADS; i++)
    pthread_create(threads+i,NULL,parm[0].output_thread,parm+i);

  for (i = 0; i < ITHREADS; i++)
    pthread_join(threads[i],NULL);
#endif

#ifdef DEBUG_OUT
  exit (0);
#endif
  free(bases);
  free(boff);
}

  //  Free an Input_Partition data structure

void Free_Input_Partition(Input_Partition *parts)
{ Thread_Arg *parm = (Thread_Arg *) parts;
  int nfiles, gz_ok;
  int i, f;

  nfiles = parm[ITHREADS-1].eidx + 1;
  gz_ok = (nfiles >= 3 || nfiles >= NTHREADS/2);

  if (parm[0].decomp != NULL)
    for (i = 0; i < ITHREADS; i++)
      libdeflate_free_decompressor(parm[i].decomp);
  free(parm[0].buf);
  for (f = 0; f <= parm[ITHREADS-1].eidx; f++)
    Free_File(parm[0].fobj+f,gz_ok);
  free(parm[0].fobj);
  free(parm);
}
