/*******************************************************************************************
 *
 *  Phase 4 of FastK: Given NPARTS sorted super-mer profiles for NTHREADS in NPANELS
 *    parts per thread, merge the super-mer NPARTS x NPANELS files for each thread
 *    (in subdirectory SORT_PATH) into single compressed read-profile files (in
 *    subdirectory "path").
 *
 *  Author:  Gene Myers
 *  Date  :  Oectober 2020
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <fcntl.h>
#include <pthread.h>

#include <errno.h>

#include "gene_core.h"
#include "FastK.h"

#undef    DEBUG
#undef    DEBUG_SEGS
#undef    SHOW_RUN

#define THREAD  pthread_t

  //  IO block data structure and block fetcher

static int64 BUFLEN_UINT8;
static int64 BUFLEN_INT64;
static int64 BUFLEN_IBYTE;

typedef struct
  { int     stream;
    int     panel;
    uint8  *block;
    uint8  *ptr;
    uint8  *top;
    int64   rid;
    int     lst;
  } IO_block;

static void reload(IO_block *in)
{ int64 del;

  del = in->top - in->ptr;
  if (del > 0)
    memmove(in->block, in->ptr, del);
  in->ptr  = in->block;
  in->top  = in->block + del;
  in->top += read(in->stream,in->top,BUFLEN_UINT8-del);
}

  //  Thread to merge files of super-mer profiles

static int PAN_SIZE;   //  # of profile runs to output as a block

typedef struct
  { int    len;   //  length of super-mer profile (including length bytes)
    int    lst;   //  ends a read
    uint8 *frag;  //  encoding of super-mer profile
  } Entry;

typedef struct
  { char     *root;    //  Path & prefix of each part file
    IO_block *io;      //  input & output buffers
    Entry    *chord;   //  super-mer profiles vector
    int       wch;     //  Number of this thread
    int       afile;   //  A-file output
    int       dfile;   //  D-file output
    int64     nreads;
  } Track_Arg;

static int64 totin;  //  Total bytes for part 0, thread 0

static void *merge_profile_thread(void *arg)
{ Track_Arg   *data  = (Track_Arg *) arg;
  int          afile = data->afile;
  int          dfile = data->dfile;
  IO_block    *io    = data->io;
  Entry       *chord = data->chord;
  int          maxS  = 2*MAX_SUPER;

  uint8 *dbuf;
  int64 *abuf, *aptr, *atop;
  int64  offset;
  int    nreads;
  char  *fname;
  int64  panel, nanel;
  int    naval, loaded, wlast;
  uint16 lcont, d;
  uint8 *db = (uint8 *)  &d;
  int    n;

  int64  pct1, partin, nextin;
  int    CLOCK;

  int    len = 0;
  int64  rid = 0;
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
  uint8  *lb = ((uint8 *) &len) - 1;
  uint8  *rb = ((uint8 *) &rid) - 1;
#else
  uint8  *lb = ((uint8 *) &len) + (sizeof(int)-PLEN_BYTES);
  uint8  *rb = ((uint8 *) &rid)) + (sizeof(int64)-RUN_BYTES);
#endif

#ifdef DEBUG
  printf("THREAD %d\n",data->wch);
#endif

  fname = Malloc(strlen(data->root) + 100,"File name buffer");
  if (fname == NULL)
    exit (1);

  //  Set up clock if required

  nextin = pct1 = partin = 0;
  if (VERBOSE && data->wch == 0)
    { nextin = pct1 = totin/100;
      fprintf(stderr,"\n    0%%");
      fflush(stderr);
      CLOCK = 1;
    }
  else
    CLOCK = 0;

  //  Setup A-file and D-file buffers

  abuf = (int64 *) io[NPARTS].block;
  aptr = abuf;
  atop = abuf + BUFLEN_INT64;

  dbuf = chord[0].frag;

  //  Read the first run-id for each part

  offset = 0;
  nreads = 0;

  wlast = 1;
  lcont = 0;
  naval = NPARTS;
  panel = 0x7fffffffffffffffll;
  for (n = 0; n < NPARTS; n++)
    { IO_block *src;
      uint8    *sptr;
      int       f;
  
      src        = io+n;
      sptr       = src->block;
      src->ptr   = src->top = sptr;
      src->panel = -1;

      //  Open panels until get to a non-empty one or all are exhausted

      while (sptr >= src->top)
        { if (src->panel >= 0)
            { close(src->stream);
#ifndef DEVELOPER
              sprintf(fname,"%s.%d.P%d.%d",data->root,n,data->wch,src->panel);
              unlink(fname);
#endif
#ifdef DEBUG
              printf("Closing Panel %d:%d.%d\n",n,data->wch,src->panel);
#endif
            }

          src->panel += 1;
          if (src->panel >= NPANELS)
            break;

#ifdef DEBUG
          printf("Starting Panel %d:%d.%d\n",n,data->wch,src->panel);
          fflush(stdout);
#endif
          sprintf(fname,"%s.%d.P%d.%d",data->root,n,data->wch,src->panel);
          f = open(fname,O_RDONLY);
          if (f == -1)
            { fprintf(stderr,"%s: Cannot open external file %s in %s\n",
                             Prog_Name,fname,SORT_PATH);
              exit (1);
            }

#ifdef DEVELOPER
          if (src->panel == 0 && data->wch == 0)
            lseek(f,3*sizeof(int),SEEK_SET);
#endif
          src->top    = sptr + read(f,sptr,BUFLEN_UINT8);
        }

      //  If panels exhausted then part is exhausted (has no data)

      if (src->panel >= NPANELS)
        { naval -= 1;
          src->rid = 0x7fffffffffffffffll;
          continue;
        }

      //  Scan run-id

      if (*sptr & 0x80)
        { src->lst = 1;
          *sptr &= 0x7f;
        }
      else
        src->lst = 0;
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
      for (f = RUN_BYTES; f > 0; f--)
        rb[f] = *sptr++;
#else
      for (f = 0; f < RUN_BYTES; f++)
        rb[f] = *sptr++;
#endif
      if (rid < panel)
        panel = rid;
      src->rid = rid;
      src->ptr = sptr;
    }

  //  For each run of PAN_SIZE consecutive super-mer profiles do

  for ( ; naval > 0; panel = nanel)
    { nanel  = panel + PAN_SIZE;
      loaded = 0;

#ifdef DEBUG_SEGS
      printf("Panning %lld - %lld\n",panel,nanel);
      for (n = 0; n < PAN_SIZE; n++)
        chord[n].len = 0;
#endif

      //  Load profiles in current run range [panel,nanel) from the parts part so as to fill range

      for (n = 0; n < NPARTS; n++)
        { IO_block *src;
          uint8    *sptr;
          int       f;

          src  = io+n;
          sptr = src->ptr;
          rid  = src->rid;
          while (rid < nanel)    //  While next profile is in range
            { loaded += 1;                      //  Load it
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
              for (f = PLEN_BYTES; f > 0; f--)
                lb[f] = *sptr++;
#else
              for (f = 0; f < PLEN_BYTES; f++)
                lb[f] = *sptr++;
#endif
              if (len > 2)
                len += 2;
#ifdef DEBUG_SEGS
              printf("   %lld: %d(%d) %d\n",rid,len,src->lst,n);
#endif
              rid -= panel;
              chord[rid].len = len;
              chord[rid].lst = src->lst;
              memcpy(chord[rid].frag,sptr,len);
              sptr += len;

              //  Fetch more data from panel if needed

              if (sptr + maxS >= src->top && src->top - src->block >= BUFLEN_UINT8)
                { src->ptr = sptr;
                  if (CLOCK && n == 0)
                    { partin += src->ptr - src->block;
                      reload(src);
                    }
                  else
                    reload(src);
                  sptr = src->ptr;
                }

              //  If panel EOF then open next panel(s) to get next run-id

              while (sptr >= src->top)
                { int f;

                  if (CLOCK && n == 0 && src->top - src->block < BUFLEN_UINT8)
                    partin += src->top - src->block;
                  close(src->stream);
#ifndef DEVELOPER
                  sprintf(fname,"%s.%d.P%d.%d",data->root,n,data->wch,src->panel);
                  unlink(fname);
#endif
                  src->panel += 1;
                  if (src->panel >= NPANELS)
                    break;

#ifdef DEBUG
                  printf("Starting Panel %d:%d.%d\n",n,data->wch,src->panel);
                  fflush(stdout);
#endif
                  sprintf(fname,"%s.%d.P%d.%d",data->root,n,data->wch,src->panel);
                  f = open(fname,O_RDONLY);
                  if (f == -1)
                    { fprintf(stderr,"%s: A Cannot open external file %s in %s\n",
                                     Prog_Name,fname,SORT_PATH);
                      exit (1);
                    }
                  sptr = src->block;
                  src->stream = f;
                  src->top    = sptr + read(f,sptr,BUFLEN_UINT8);
                }

              //  If all panels exhausted then thread is exhausted

              if (src->panel >= NPANELS)
                { naval -= 1;
                  rid = 0x7fffffffffffffffll;
                  break;
                }

              //  Scan run-id
             
              if (*sptr & 0x80)
                { src->lst = 1;
                  *sptr &= 0x7f;
                }
              else
                src->lst = 0;
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
              for (f = RUN_BYTES; f > 0; f--)
                rb[f] = *sptr++;
#else
              for (f = 0; f < RUN_BYTES; f++)
                rb[f] = *sptr++;
#endif
#ifdef DEBUG_SEGS
              printf("   %lld(%d) <%d> [%d]\n",rid,src->lst,n,data->wch);
              fflush(stdout);
#endif
            }
          src->rid = rid;
          src->ptr = sptr;
        }

#ifdef DEBUG_SEGS
      if (loaded < PAN_SIZE)
        { printf("Loaded %d (%d)\n",loaded,data->wch);
          for (n = 0; n < PAN_SIZE; n++)
            if (chord[n].len > 0)
              printf("  %d -- %lld\n",n,panel+n);
        }
      fflush(stdout);
#endif

      //  Compress and join the profile fragments in the current range and output

      { uint8 *o, *ptr;

        o = dbuf;
        for (n = 0; n < loaded; n++)
          { len = chord[n].len;
            ptr = chord[n].frag;

            //  Get 1st count of profile n

            d = *((uint16 *) ptr);
            ptr += 2;
#ifdef SHOW_RUN
            if (wlast)
              printf("READ\n  %5d:: %3d: {%hu}",n,len,d);
            else
              printf("  %5d:: %3d: {%hd}",n,len,d);
#endif

            //  If start of read then encode as absolute

            if (wlast)
              { if (d < 128)
                  *o++ = d;
                else
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
                  { *o++ = db[1] | 0x80;
                    *o++ = db[0];
#else
                  { *o++ = db[0] | 0x80;
                    *o++ = db[1];
#endif
#ifdef SHOW_RUN
                    printf("+");
#endif
                  }
                lcont = d;
              }

            //  Otherwise encode 1st forward difference

            else
              { d -= lcont;
                if (d == 0)
                  *o++ = 0x01;
                else if (d > 0xffe1 || d < 32)
                  *o++ = 0x40 | (d & 0x3f);
                else
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
                  { *o++ = db[1] | 0x80;
                    *o++ = db[0];
#else
                  { *o++ = db[0] | 0x80;
                    *o++ = db[1];
#endif
#ifdef SHOW_RUN
                    printf("+");
#endif
                  }
                lcont += d;
              } 

            //  Transfer remainder of profile keeping track of absolute count
          
            if (len > 2)
              { len -= 4;
                memmove(o,ptr,len);
                lcont = *((uint16 *) (ptr+len));
#ifdef SHOW_RUN
                for (int u = 0; u < len; u++)
                  printf(" %02x",ptr[u]);
                printf(" {%d}",lcont);
#endif
                o += len;
              }

#ifdef SHOW_RUN
            printf("\n");
#endif

            //  If last profile of a read, then add to A-file index

            wlast = chord[n].lst;
            if (wlast)
              { if (aptr >= atop)
                  { write(afile,abuf,BUFLEN_IBYTE);
                    aptr = abuf;
                  }
                *aptr++ = offset + (o-dbuf);
                nreads += 1;
              }
          }

        //  Write the merged, compressed profiles for the range

        write(dfile,dbuf,o-dbuf);
        offset += o-dbuf;

        //  Check clock

        if (CLOCK && partin + (io->ptr - io->block) > nextin)
          { nextin = partin + (io->ptr - io->block);
            fprintf(stderr,"\r  %3d%%",(int) ((100.*nextin)/totin));
            fflush(stderr);
            nextin += pct1;
          }
      }
    }

  //  Flush the A-file buffer

  if (aptr > abuf)
    write(afile,abuf,(aptr-abuf)*sizeof(int64));

  data->nreads = nreads;

  free(fname);

  if (CLOCK)
    fprintf(stderr,"\r         \r");

  return (NULL);
}

void Merge_Profiles(char *dpwd, char *dbrt)
{ char       *fname;
  IO_block   *io;
  uint8      *blocks;
  Entry      *chord;
  uint8      *_chord;

  if (VERBOSE)
    { fprintf(stderr,"\nPhase 4 (-p option): Merging Profile Fragments\n");
      fflush(stderr);
    }

  //  Chek that input files are at sort path and determine the number of panels,
  //     threads, and partitions

  fname = Malloc(strlen(dpwd) + strlen(SORT_PATH) + strlen(dbrt) + 100,"File name buffer");
  if (fname == NULL)
    exit (1);

  //  Get rid of any previous results for this DB in this directory with this KMER

  sprintf(fname,"rm -f %s/%s.K%d.A*",dpwd,dpwd,KMER);
  system(fname);
  sprintf(fname,"rm -f %s/%s.K%d.P*",dpwd,dpwd,KMER);
  system(fname);

  //  Allocate all working data structures

  { int         f;
    struct stat info;

    sprintf(fname,"%s/%s.0.P0.0",SORT_PATH,dbrt);
    f = open(fname,O_RDONLY);
    if (f == -1)
      { fprintf(stderr,"%s: A Cannot open external file %s in %s\n",
                       Prog_Name,fname,SORT_PATH);
        exit (1);
      }
  
#ifdef DEVELOPER
     read(f,&RUN_BYTES,sizeof(int));
     read(f,&PLEN_BYTES,sizeof(int));
     read(f,&MAX_SUPER,sizeof(int));
#endif

     if (VERBOSE)
       { fstat(f,&info);
         totin = info.st_size;

         close(f);

         for (int i = 1; i < NPANELS; i++)
           { sprintf(fname,"%s/%s.0.P0.%d",SORT_PATH,dbrt,i);
             stat(fname,&info);
             totin += info.st_size;
           }
       }
     else
       close(f);
   }

  BUFLEN_UINT8 = SORT_MEMORY/((NPARTS+1)*NTHREADS);
  if (BUFLEN_UINT8 > 0x7fffffffll)
    BUFLEN_UINT8 = 0x7ffffff8ll;
  if (BUFLEN_UINT8 < 2*MAX_SUPER)
    BUFLEN_UINT8 = 2*MAX_SUPER;

  BUFLEN_INT64 = BUFLEN_UINT8 / sizeof(int64);
  BUFLEN_IBYTE = BUFLEN_INT64 * sizeof(int64);
  PAN_SIZE     = 1024*NPARTS;

  io     = (IO_block *) Malloc(sizeof(IO_block)*(NPARTS+1)*NTHREADS,"Allocating IO buffers");
  blocks = (uint8 *) Malloc(BUFLEN_UINT8*(NPARTS+1)*NTHREADS,"Allocating IO buffers");
  chord  = (Entry *) Malloc(PAN_SIZE*sizeof(Entry)*NTHREADS,"Allocating IO buffers");
  _chord = (uint8 *) Malloc(PAN_SIZE*2*MAX_SUPER*NTHREADS,"Allocating IO buffers");
  if (io == NULL || blocks == NULL || chord == NULL || _chord == NULL)
    exit (1);

  //  Open up A- and D-files, assign blocks for the inputs, and setup thread params 

  { Track_Arg   parmk[NTHREADS];
#ifndef DEBUG
    THREAD      threads[NTHREADS];
#endif
    char *root;
    int   t, n, p;
    
    root = Malloc(strlen(SORT_PATH) + strlen(dbrt) + 10,"File name buffer");
    if (root == NULL)
      exit (1);
    sprintf(root,"%s/%s",SORT_PATH,dbrt);

    p = 0;
    for (t = 0; t < NTHREADS; t++)
      { int   f, g;
        int64 zero = 0;
            
        sprintf(fname,"%s/%s.K%d.A%d",dpwd,dbrt,KMER,t+1);
        f = open(fname,O_WRONLY|O_CREAT|O_TRUNC,S_IRWXU|S_IRWXG|S_IRWXO);
        if (f == -1)
          { fprintf(stderr,"%s: Cannot open external file %s in %s\n",
                           Prog_Name,fname,SORT_PATH);
            exit (1);
          }
        sprintf(fname,"%s/%s.K%d.P%d",dpwd,dbrt,KMER,t+1);
        g = open(fname,O_WRONLY|O_CREAT|O_TRUNC,S_IRWXU|S_IRWXG|S_IRWXO);
        if (g == -1)
          { fprintf(stderr,"%s: Cannot open external file %s in %s\n",
                           Prog_Name,fname,SORT_PATH);
            exit (1);
          }
        
        parmk[t].root  = root;
        parmk[t].wch   = t;
        parmk[t].afile = f;
        parmk[t].dfile = g;
        parmk[t].io    = io + p;
        parmk[t].chord = chord + PAN_SIZE*t;


        for (n = 0; n <= NPARTS; n++, p++)
          io[p].block = blocks + p*BUFLEN_UINT8;

        write(f,&KMER,sizeof(int));
        write(f,&zero,sizeof(int64));
        write(f,&zero,sizeof(int64));
      }

    //  Setup the fragment buffers for each range chord

    for (n = 0; n < PAN_SIZE*NTHREADS; n++)
      chord[n].frag = _chord + 2*MAX_SUPER*n;

    //  In parallel process each of the NTHREADS partitions

#ifdef DEBUG
    for (t = 0; t < NTHREADS; t++)
      merge_profile_thread(parmk+t);
#else
    for (t = 1; t < NTHREADS; t++)
      pthread_create(threads+t,NULL,merge_profile_thread,parmk+t);
    merge_profile_thread(parmk);
    for (t = 1; t < NTHREADS; t++)
      pthread_join(threads[t],NULL);
#endif

    //  Rewind and set the header of each A-file

    { int64 nreads;
      int   f;

      nreads = 0;
      for (t = 0; t < NTHREADS; t++)
        { f = parmk[t].afile;

          lseek(f,sizeof(int),SEEK_SET);
          write(f,&nreads,sizeof(int64));
          write(f,&(parmk[t].nreads),sizeof(int64));
          nreads += parmk[t].nreads;

          close(f);
          close(parmk[t].dfile);
        }
    }

    //  Release working data

    free(root);
    free(_chord);
    free(chord);
    free(blocks);
    free(io);
    free(fname);
  }
}
