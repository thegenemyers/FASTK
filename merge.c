/*******************************************************************************************
 *
 *  Phase 4 of FastK: Given NPARTS sorted super-mer profiles for ITHREADS in NPANELS
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

#include "libfastk.h"
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
    char     *aname;
    int       dfile;   //  D-file output
    char     *dname;
    int64     nreads;  //  # of reads seen by this thread
    int64     nbase;   //  First rid for thread
  } Track_Arg;

static int64 totin;  //  Total bytes for part 0, thread 0

static void *merge_profile_thread(void *arg)
{ Track_Arg   *data  = (Track_Arg *) arg;
  int          afile = data->afile;
  int          dfile = data->dfile;
  IO_block    *io    = data->io;
  Entry       *chord = data->chord;
  int64        nbase = data->nbase;
  int          maxS  = 2*MAX_SUPER + 2 + PLEN_BYTES + RUN_BYTES;

  uint8 *dbuf;
  int64 *abuf, *aptr, *atop;
  int64  offset;
  int    nreads;
  char  *fname;
  int64  panel, nanel;
  int    naval, wlast;
  uint16 lcont, d;
  uint8 *db = (uint8 *)  &d;
  int    n;
  int64  nidx = 0;

  FILE *nfile;  //  Invalid k-mers file
  int64 iridx;
  int   ileng, ilast = 0;

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
    Clean_Exit(1);

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

  //  Open invalid interval file

  sprintf(fname,"%s.NS.T%d",data->root,data->wch);
  nfile = fopen(fname,"r");
  if (nfile == NULL)
    { fprintf(stderr,"\n%s: Cannot open external file %s in %s\n",
                     Prog_Name,fname,SORT_PATH);
      Clean_Exit(1);
    }
  if (fread(&iridx,sizeof(int64),1,nfile) < 1)
    iridx = 0x7fffffffffffffffll;
  else
    { if ((iridx & 0x8000000000000000ll) != 0)
        { ilast = 1;
          iridx &= 0x7fffffffffffffffll;
        }
      else
        ilast = 0;
      fread(&ileng,sizeof(int),1,nfile);
    }

  //  Read the first run-id for each part

  offset = 0;
  nreads = 0;

  wlast = 1;
  lcont = 0;
  naval = NPARTS;
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
            { fprintf(stderr,"\n%s: Cannot open external file %s in %s\n",
                             Prog_Name,fname,SORT_PATH);
              Clean_Exit(1);
            }

#ifdef DEVELOPER
          if (src->panel == 0)
            { if (data->wch == 0)
                lseek(f,3*sizeof(int),SEEK_SET);
              read(f,&nidx,sizeof(int64));
            }
#else
          nidx = NUM_RID[data->wch];
#endif
          src->stream = f;
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
      src->rid = rid;
      src->ptr = sptr;
    }

  //  For each run of PAN_SIZE consecutive super-mer profiles do

  nidx += nbase;
  for (panel = nbase; naval > 0; panel = nanel)
    { nanel  = panel + PAN_SIZE;
      if (nanel > nidx)
        nanel = nidx;

#ifdef DEBUG_SEGS
      printf("Panning %lld - %lld\n",panel,nanel);
#endif
      for (n = 0; n < PAN_SIZE; n++)
        chord[n].len = 0;

      //  Load profiles in current run range [panel,nanel) from the parts part so as to fill range

      for (n = 0; n < NPARTS; n++)
        { IO_block *src;
          uint8    *sptr;
          int       f;

          src  = io+n;
          sptr = src->ptr;
          rid  = src->rid;
          while (rid < nanel)    //  While next profile is in range
            {
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
                    { fprintf(stderr,"\n%s: A Cannot open external file %s in %s\n",
                                     Prog_Name,fname,SORT_PATH);
                      Clean_Exit(1);
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
      if (nanel-panel < PAN_SIZE)
        { printf("Loaded %lld (%d)\n",nanel-panel,data->wch);
          for (n = 0; n < PAN_SIZE; n++)
            if (chord[n].len > 0)
              printf("  %d -- %lld\n",n,panel+n);
        }
      fflush(stdout);
#endif

      //  Compress and join the profile fragments in the current range and output

      { uint8 *o, *ptr;
        int64  p;
        uint8  lz;

        lz = 0;
        o = dbuf;
        for (n = 0, p = panel; p < nanel; n++, p++)
          { len = chord[n].len;
            ptr = chord[n].frag;

            if (len == 0)
              { if (iridx+nbase != p)           //  read < KMER bp long
                  { if (aptr >= atop)
                      { if (write(afile,abuf,BUFLEN_IBYTE) < 0)
                          { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",
                                           Prog_Name,data->aname);
                            Clean_Exit(1);
                          }
                        aptr = abuf;
                      }
                    *aptr++ = offset + (o-dbuf);
#ifdef SHOW_RUN
                    if (wlast)
                      printf("READ %d\n  %5d:: SHORT\n",nreads,n);
                    else
                      printf("  %5d:: SHORT\n",n);
#endif
                    nreads += 1;
                  }
                else
                  { if (wlast)
                      { *o++ = 0;
                        lz = 0;
#ifdef SHOW_RUN
                        printf("READ %d\n  %5d:: N=%d {0} [%02x]",nreads,n,ileng,o[-1]);
#endif
                      }
                    else
                      { d = -lcont;
#ifdef SHOW_RUN
                        printf("  %5d:: N=%d {0}",n,ileng);
#endif
                        if (d == 0)
                          { lz += 1;
                            if (lz >= 63)
                              { *o++ = 63;
                                lz = 0;
#ifdef SHOW_RUN
                                printf(" @");
#endif
                              }
                          }
                        else
                          { if (lz)
                              { *o++ = lz;
#ifdef SHOW_RUN
                                printf(" %02x",lz);
#endif
                                lz = 0;
                              }
                            if (d > 0xffe1 || d < 32)
                              { *o++ = 0x40 | (d & 0x3f);
#ifdef SHOW_RUN
                                printf(" %02x",o[-1]);
#endif
                              }
                            else
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
                              { *o++ = db[1] | 0x80;
                                *o++ = db[0];
#else
                              { *o++ = db[0] | 0x80;
                                *o++ = db[1];
#endif
#ifdef SHOW_RUN
                                printf(" %02x.%02x",o[-2],o[-1]);
#endif
                              }
                          }
                      } 

                    if (ileng > 1)
                      { ileng += lz;
                        for (ileng -= 1; ileng >= 63; ileng -= 63)
                          { *o++ = 63;
#ifdef SHOW_RUN
                            printf(" @");
#endif
                          }
                        lz = ileng;
#ifdef SHOW_RUN
                        printf(" {0} <%d>\n",lz);
#endif
                      }
#ifdef SHOW_RUN
                    else
                      printf(" {0} <%d>\n",lz);
#endif
                    lcont = 0;

                    if (ilast)
                      { if (lz)
                          { *o++ = lz;
#ifdef SHOW_RUN
                            printf("  TAILI:: %02x\n",lz);
#endif
                            lz = 0;
                          }
                        if (aptr >= atop)
                          { if (write(afile,abuf,BUFLEN_IBYTE) < 0)
                              { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",
                                               Prog_Name,data->aname);
                                Clean_Exit(1);
                              }
                            aptr = abuf;
                          }
                        *aptr++ = offset + (o-dbuf);
                        nreads += 1;
                        wlast = 1;
                      }
                    else
                      wlast = 0;

                    if (fread(&iridx,sizeof(int64),1,nfile) < 1)
                      iridx = 0x7fffffffffffffffll;
                    else
                      { if ((iridx & 0x8000000000000000ll) != 0)
                          { ilast = 1;
                            iridx &= 0x7fffffffffffffffll;
                          }
                        else
                          ilast = 0;
                        fread(&ileng,sizeof(int),1,nfile);
                      }
                  }
                continue;
              }

            //  Get 1st count of profile n

            d = *((uint16 *) ptr);
            ptr += 2;

            //  If start of read then encode as absolute

            if (wlast)
              { if (d < 128)
                  { *o++ = d;
#ifdef SHOW_RUN
                    printf("READ %d\n  %5d:: %3d: {%hu} [%02x]",nreads,n,len,d,o[-1]);
#endif
                  }
                else
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
                  { *o++ = db[1] | 0x80;
                    *o++ = db[0];
#else
                  { *o++ = db[0] | 0x80;
                    *o++ = db[1];
#endif
#ifdef SHOW_RUN
                    printf("READ %d\n  %5d:: %3d: {%hu} [%02x.%02x]",nreads,n,len,d,o[-2],o[-1]);
#endif
                  }
                lcont = d;
                lz = 0;
              }

            //  Otherwise encode 1st forward difference

            else
              {
#ifdef SHOW_RUN
                printf("  %5d:: %3d: {%hd}",n,len,d);
#endif
                d -= lcont;
                if (d == 0)
                  { lz += 1;
                    if (lz >= 63)
                      { *o++ = 63;
                        lz = 0;
#ifdef SHOW_RUN
                        printf(" @");
#endif
                      }
                  }
                else
                  { if (lz)
                      { *o++ = lz;
#ifdef SHOW_RUN
                        printf(" %02x",lz);
#endif
                        lz = 0;
                      }
                    if (d > 0xffe1 || d < 32)
                      { *o++ = 0x40 | (d & 0x3f);
#ifdef SHOW_RUN
                        printf(" %02x",o[-1]);
#endif
                      }
                    else
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
                      { *o++ = db[1] | 0x80;
                        *o++ = db[0];
#else
                      { *o++ = db[0] | 0x80;
                        *o++ = db[1];
#endif
#ifdef SHOW_RUN
                        printf(" %02x.%02x",o[-2],o[-1]);
#endif
                      }
                  } 
                lcont += d;
              } 

            //  Transfer remainder of super-mer profile
          
            if (len > 2)
              { len -= 4;
                lcont = *((uint16 *) (ptr+len));

                if (lz)
                  { while (len > 0 && *ptr < 64)
                      { lz += *ptr;
                        if (lz >= 63)
                          { *o++ = 63;
#ifdef SHOW_RUN
                            printf(" @");
#endif
                            lz -= 63;
                          }
                        ptr += 1;
                        len -= 1;
                      }
                  }

                if (len > 0)
                  { if (lz)
                      { *o++ = lz;
#ifdef SHOW_RUN
                        printf(" %02x",lz);
#endif
                        lz = 0;
                      }

                    { int u, one;

                      one = 1;
                      for (u = 0; u < len; u++)
                        if ((ptr[u] & 0x80) != 0)
                          { one = 0;
                            u += 1;
                          }
                        else
                          one = 1;

                      if (one)
                        { lz = ptr[len-1];
                          if (lz < 63)
                            len -= 1;
                          else
                            lz = 0;
                        }
                    }

                    memmove(o,ptr,len);
#ifdef SHOW_RUN
                    { int u;

                      for (u = 0; u < len; u++)
                        if ((ptr[u] & 0x80) != 0)
                          { printf(" %02x.%02x",ptr[u],ptr[u+1]);
                            u += 1;
                          }
                        else
                          printf(" %02x",ptr[u]);
                    }
#endif
                    o += len;
                  }
              }
#ifdef SHOW_RUN
            printf(" {%d}",lcont);
            if (lz)
              printf(" <%d>",lz);
            printf("\n");
#endif

            //  If last profile of a read, then add to A-file index

            wlast = chord[n].lst;
            if (wlast)
              { if (lz)
                  { *o++ = lz;
#ifdef SHOW_RUN
                    printf("  TAILZ:: %02x\n",lz);
#endif
                    lz = 0;
                  }
                if (aptr >= atop)
                  { if (write(afile,abuf,BUFLEN_IBYTE) < 0)
                      { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",
                                       Prog_Name,data->aname);
                        Clean_Exit(1);
                      }
                    aptr = abuf;
                  }
                *aptr++ = offset + (o-dbuf);
                nreads += 1;
              }
          }

        //  Write the merged, compressed profiles for the range

        if (lz)
          { *o++ = lz;
#ifdef SHOW_RUN
            printf("  TAILE:: %02x\n",lz);
#endif
          }
        if (write(dfile,dbuf,o-dbuf) < 0)
          { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",
                           Prog_Name,data->dname);
             Clean_Exit(1);
          }
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
    if (write(afile,abuf,(aptr-abuf)*sizeof(int64)) < 0)
      { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",
                       Prog_Name,data->aname);
         Clean_Exit(1);
      }

  data->nreads = nreads;

  fclose(nfile);

#ifndef DEVELOPER
  sprintf(fname,"%s.NS.T%d",data->root,data->wch);
  unlink(fname);
#endif

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
  Track_Arg  *parmk;
#ifndef DEBUG
  THREAD     *threads;
#endif

  if (VERBOSE)
    { fprintf(stderr,"\nPhase 4 (-p option): Merging Profile Fragments\n");
      fflush(stderr);
    }

  //  Chek that input files are at sort path and determine the number of panels,
  //     threads, and partitions

  fname = Malloc(3*(strlen(dpwd) + strlen(SORT_PATH) + strlen(dbrt)) + 100,"File name buffer");
  if (fname == NULL)
    Clean_Exit(1);

  //  Allocate all working data structures

  parmk = Malloc(sizeof(Track_Arg)*ITHREADS,"Allocating control data");
  if (parmk == NULL)
    Clean_Exit(1);

#ifndef DEBUG
  threads = Malloc(sizeof(THREAD)*ITHREADS,"Allocating control data");
  if (threads == NULL)
    Clean_Exit(1);
#endif

  { int   t;
    int64 o, n;

    o = 0;
    for (t = 0; t < ITHREADS; t++)
      { int         f, i;
        struct stat info;

        sprintf(fname,"%s/%s.0.P%d.0",SORT_PATH,dbrt,t);
        f = open(fname,O_RDONLY);
        if (f == -1)
          { fprintf(stderr,"%s: A Cannot find file %s.0.P%d.0 in directory %s\n",
                           Prog_Name,dbrt,t,SORT_PATH);
            Clean_Exit(1);
          }
  
#ifdef DEVELOPER
        if (t == 0)
          { read(f,&RUN_BYTES,sizeof(int));
            read(f,&PLEN_BYTES,sizeof(int));
            read(f,&MAX_SUPER,sizeof(int));
          }
        read(f,&n,sizeof(int64));
#else
        n = NUM_RID[t];
#endif
        parmk[t].nbase = o;
        o += n;

        if (VERBOSE && t == 0)
          { fstat(f,&info);
            totin = info.st_size;
  
            for (i = 1; i < NPANELS; i++)
              { sprintf(fname,"%s/%s.0.P0.%d",SORT_PATH,dbrt,i);
                stat(fname,&info);
                totin += info.st_size;
              }
          }
        close(f);
      }
  }

  BUFLEN_UINT8 = SORT_MEMORY/((NPARTS+1)*ITHREADS);
  if (BUFLEN_UINT8 > 0x7fffffffll)
    BUFLEN_UINT8 = 0x7ffffff8ll;
  if (BUFLEN_UINT8 < 2*MAX_SUPER)
    BUFLEN_UINT8 = 2*MAX_SUPER;

  BUFLEN_INT64 = BUFLEN_UINT8 / sizeof(int64);
  BUFLEN_IBYTE = BUFLEN_INT64 * sizeof(int64);
  PAN_SIZE     = 1024*NPARTS;

  io     = (IO_block *) Malloc(sizeof(IO_block)*(NPARTS+1)*ITHREADS,"Allocating IO buffers");
  blocks = (uint8 *) Malloc(BUFLEN_UINT8*(NPARTS+1)*ITHREADS,"Allocating IO buffers");
  chord  = (Entry *) Malloc(PAN_SIZE*sizeof(Entry)*ITHREADS,"Allocating IO buffers");
  _chord = (uint8 *) Malloc(PAN_SIZE*2*(MAX_SUPER+1)*ITHREADS,"Allocating IO buffers");
  if (io == NULL || blocks == NULL || chord == NULL || _chord == NULL)
    Clean_Exit(1);

  //  Remove previous profile result if any

  sprintf(fname,"rm -f %s/%s.prof %s/.%s.pidx.* %s/.%s.prof.*",dpwd,dbrt,dpwd,dbrt,dpwd,dbrt);
  system(fname);

  //  Create new stub file

  { int f;

    sprintf(fname,"%s/%s.prof",dpwd,dbrt);
    f = open(fname,O_WRONLY|O_CREAT|O_TRUNC,S_IRWXU|S_IRWXG|S_IRWXO);
    if (f == -1)
      { fprintf(stderr,"%s: Cannot open external file %s for writing\n",Prog_Name,fname);
        Clean_Exit(1);
      }
    write(f,&KMER,sizeof(int));
    if (write(f,&ITHREADS,sizeof(int)) < 0)
      { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",Prog_Name,fname);
         Clean_Exit(1);
      }
    close(f);
  }

  //  Open up A- and D-files, assign blocks for the inputs, and setup thread params 

  { char *root;
    char *aname, *dname;
    int   t, n, p;
    
    root = Malloc(strlen(SORT_PATH) + strlen(dbrt) + 10,"File name buffer");
    if (root == NULL)
      Clean_Exit(1);
    sprintf(root,"%s/%s",SORT_PATH,dbrt);

    p = 0;
    for (t = 0; t < ITHREADS; t++)
      { int   f, g;
        int64 zero = 0;
            
        sprintf(fname,"%s/.%s.pidx.%d",dpwd,dbrt,t+1);
        f = open(fname,O_WRONLY|O_CREAT|O_TRUNC,S_IRWXU|S_IRWXG|S_IRWXO);
        if (f == -1)
          { fprintf(stderr,"%s: Cannot open external file %s\n",Prog_Name,fname);
            Clean_Exit(1);
          }
        aname = Strdup(fname,"Allocating stream names");

        sprintf(fname,"%s/.%s.prof.%d",dpwd,dbrt,t+1);
        g = open(fname,O_WRONLY|O_CREAT|O_TRUNC,S_IRWXU|S_IRWXG|S_IRWXO);
        if (g == -1)
          { fprintf(stderr,"%s: Cannot open external file %s\n",Prog_Name,fname);
            Clean_Exit(1);
          }
        dname = Strdup(fname,"Allocating stream names");

        if (aname == NULL || dname == NULL)
          Clean_Exit(1);
        
        parmk[t].root  = root;
        parmk[t].wch   = t;
        parmk[t].afile = f;
        parmk[t].aname = aname;
        parmk[t].dfile = g;
        parmk[t].dname = dname;
        parmk[t].io    = io + p;
        parmk[t].chord = chord + PAN_SIZE*t;

        for (n = 0; n <= NPARTS; n++, p++)
          io[p].block = blocks + p*BUFLEN_UINT8;

        write(f,&KMER,sizeof(int));
        write(f,&zero,sizeof(int64));
        if (write(f,&zero,sizeof(int64)) < 0)
          { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",Prog_Name,aname);
             Clean_Exit(1);
          }
      }

    //  Setup the fragment buffers for each range chord

    for (n = 0; n < PAN_SIZE*ITHREADS; n++)
      chord[n].frag = _chord + 2*(MAX_SUPER+1)*n;

    //  In parallel process each of the ITHREADS partitions

#ifdef DEBUG
    for (t = 0; t < ITHREADS; t++)
      merge_profile_thread(parmk+t);
#else
    for (t = 1; t < ITHREADS; t++)
      pthread_create(threads+t,NULL,merge_profile_thread,parmk+t);
    merge_profile_thread(parmk);
    for (t = 1; t < ITHREADS; t++)
      pthread_join(threads[t],NULL);
#endif

    if (VERBOSE)
      { int64 psize;

        psize = 0;
        for (t = 0; t < ITHREADS; t++)
          psize += lseek(parmk[t].afile,0,SEEK_CUR) + lseek(parmk[t].dfile,0,SEEK_CUR);

        if (psize >= 5.e8)
          fprintf(stderr,"  The profiles occupy %.2f GB\n",psize/1.e9);
        else if (psize >= 5.e5)
          fprintf(stderr,"  The profiles occupy %.2f MB\n",psize/1.e6);
        else
          fprintf(stderr,"  The profiles occupy %.2f KB\n",psize/1.e3);
        fflush(stderr);
      }

    //  Rewind and set the header of each A-file

    { int64 nreads;
      int   f;

      nreads = 0;
      for (t = 0; t < ITHREADS; t++)
        { f = parmk[t].afile;

          lseek(f,sizeof(int),SEEK_SET);
          write(f,&nreads,sizeof(int64));
          if (write(f,&(parmk[t].nreads),sizeof(int64)) < 0)
            { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",
                             Prog_Name,parmk[t].aname);
               Clean_Exit(1);
            }
          nreads += parmk[t].nreads;

          close(parmk[t].dfile);
          close(f);
          free(parmk[t].dname);
          free(parmk[t].aname);
        }
    }

    //  Release working data

    free(root);
    free(_chord);
    free(chord);
    free(blocks);
    free(io);
#ifndef DEBUG
    free(threads);
#endif
    free(parmk);
    free(fname);
  }
}
