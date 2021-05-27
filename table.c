/*******************************************************************************************
 *
 *  Phase 3 of FastK: Given NPARTS sorted k-mer count table files for NTHREADS threads,
 *    merge in k-mer order the NPARTS files (in directory SORT_PATH) for each thread
 *    into a single k-mer tables (in directory "path").
 *
 *  Author:  Gene Myers
 *  Date  :  October 2020
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
#undef    DEBUG_MERGE

#define THREAD  pthread_t


  //  Utilities for k-mers

#ifdef DEBUG

static void print_kmer(uint8 *str, int len)
{ static char DNA[4] = { 'a', 'c', 'g', 't' };
  int j, k;

  for (k = 0; k < (len+3)/4; k++)
    for (j = 6; j >= 0; j-= 2)
      printf("%c",DNA[(str[k] >> j) & 0x3]);
  for (j = 6; j > 6-2*(len%4); j-= 2)
    printf("%c",DNA[(str[k] >> j) & 0x3]);
}

#endif

static inline int mycmp(uint8 *a, uint8 *b, int n)
{ while (n--)
    { if (*a++ != *b++)
        return (a[-1] < b[-1] ? 0 : 1);
    }
  return (0);
}

static inline void mycpy(uint8 *a, uint8 *b, int n)
{ while (n--)
    *a++ = *b++;
}


  //  Input block data structure and block fetcher

static int64 BUFLEN_UINT8;   //  Size of IO buffers in bytes

typedef struct
  { int     stream;
    uint8  *block;
    uint8  *ptr;
    uint8  *top;
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


  //  Heap of input buffer pointers ordering on k-mer at ->ptr

static void reheap(int s, IO_block **heap, int hsize)
{ int       c, l, r;
  int       bigger;
  IO_block *hs, *hr, *hl;
  uint8    *hsp;

  c   = s;
  hs  = heap[s];
  hsp = hs->ptr;
  while ((l = 2*c) <= hsize)
    { r  = l+1;
      hl = heap[l];
      if (r > hsize)
        bigger = 1;
      else
        { hr = heap[r];
          bigger = mycmp(hr->ptr,hl->ptr,KMER_BYTES);
        }
      if (bigger)
        { if (mycmp(hsp,hl->ptr,KMER_BYTES))
            { heap[c] = hl;
              c = l;
            }
          else
            break;
        }
      else
        { if (mycmp(hsp,hr->ptr,KMER_BYTES))
            { heap[c] = hr;
              c = r;
            }
          else
            break;
        }
    }
  if (c != s)
    heap[c] = hs;
}

#ifdef DEBUG

static void showheap(IO_block **heap, int hsize, IO_block *io)
{ int i;
  printf("\n");
  for (i = 1; i <= hsize; i++)
    { printf("   %3d: [%2ld] ",i,heap[i]-io);
      print_kmer(heap[i]->ptr,KMER);
      printf(" %d\n",*((uint16 *) (heap[i]->ptr+KMER_BYTES)));
    }
}

#endif


  //  Thread to merge NPARTS input tables into one large table

typedef struct
  { IO_block **heap;    //  heap of (read,post) items
    IO_block  *in;      //  input & output buffers
    IO_block  *out;
    char      *oname;
    int        id;
    int64      tsize;   //  output table size in bytes
  } Track_Arg;

static int64 totin;  //  Total kmer record for thread 0 (if VERBOSE)

static int    PMER_WORD;  //  TMER_WORD - IDX_BYTES
static int64 *pindex;     //  IDX_BYTES prefix index
static int    pidxlen;    //  length of index

static void *merge_table_thread(void *arg)
{ Track_Arg   *data  = (Track_Arg *) arg;
  IO_block   **heap  = data->heap;
  IO_block    *in    = data->in;
  IO_block    *out   = data->out;
  char        *oname = data->oname;

  int    afile;
  int64  anum;
  uint8 *abuf, *aptr, *atop;
  int    hsize;
  int    p;

  int64  pct1, nextin;
  int    CLOCK;

  int    idx = 0;
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
  uint8  *idp = ((uint8 *) &idx) - 1;
#else
  uint8  *idp = ((uint8 *) &idx) + (sizeof(int)-IDX_BYTES);
#endif

  nextin = pct1 = 0;
  if (VERBOSE && data->id == 0)
    { nextin = pct1 = totin/100;
      fprintf(stderr,"\n    0%%");
      fflush(stderr);
      CLOCK = 1;
    }
  else
    CLOCK = 0;

  //  Setup output IO-block specifically

  afile = out->stream;
  abuf  = out->block;
  aptr  = abuf;
  atop  = abuf + BUFLEN_UINT8;

  //  Load 1st block of each input file

  for (p = 0; p < NPARTS; p++)
    { uint8 *iblock;

      iblock       = in[p].block;
      in[p].ptr    = iblock;
      in[p].top    = iblock + read(in[p].stream,iblock,BUFLEN_UINT8);
    }

  //  Fill output file prolog (rewind & complete at end)

  anum = 0;

  write(afile,&KMER,sizeof(int));
  if (write(afile,&anum,sizeof(int64)) < 0)
    { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",Prog_Name,oname);
      Clean_Exit(1);
    }

  //  Initialize the heap

  hsize = 0;
  for (p = 0; p < NPARTS; p++)
    if (in[p].ptr < in[p].top)
      { hsize       += 1;
        heap[hsize]  = in + p;
      }
    
  if (hsize > 3)
    for (p = hsize/2; p > 1; p--)
      reheap(p,heap,hsize);
#ifdef DEBUG
  showheap(heap,hsize,in);
#endif

  //  While the heap is not empty ...

  while (hsize > 0)
    { IO_block *src;
      uint8    *sptr;
    
      //  Get next k-mer in order from heap

      reheap(1,heap,hsize);
      src  = heap[1];
      sptr = src->ptr;

#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
      for (p = IDX_BYTES; p > 0; p--)
        idp[p] = *sptr++;
#else
      for (p = 0; p < IDX_BYTES; p++)
        idp[p] = *sptr++
#endif
      pindex[idx] += 1;

      //  Flush output buffer if needed

      if (aptr + PMER_WORD > atop)
        { int c = aptr-abuf;

#ifdef DEBUG
          printf("Write %d bytes\n",c);
#endif
          anum += c;
          if (write(afile,abuf,c) < 0)
            { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",Prog_Name,oname);
              Clean_Exit(1);
            }
          aptr = abuf;
          if (CLOCK && anum > nextin)
            { fprintf(stderr,"\r  %3d%%",(int) ((100.*anum)/totin));
              fflush(stderr);
              nextin = anum+pct1;
            }
        }

      //  Append k-mer to output

      mycpy(aptr,sptr,PMER_WORD);
      aptr += PMER_WORD;
      sptr += PMER_WORD;

#ifdef DEBUG
      printf(" %3ld:  %0*x",src-in,IDX_BYTES*2,idx);
      print_kmer(aptr-PMER_WORD,KMER-IDX_BYTES*4);
      printf(" %d",*((uint16 *) (aptr-2)));
      printf("\n");
#endif

      if (sptr + TMER_WORD > src->top)
        { src->ptr = sptr;
          reload(src);
#ifdef DEBUG
          printf("Reload %d\n",p);
#endif
          sptr = src->ptr;
        }

      //  If input is now empty, reduce heap size

      if (sptr >= src->top)
        { heap[1] = heap[hsize];
          hsize  -= 1;
#ifdef DEBUG
          printf("EOF\n");
#endif
        }
      else
        src->ptr = sptr;
    }

  //  Flush output buffer

  if (aptr > abuf)
    { int c = aptr-abuf;

      anum += c;
      if (write(afile,abuf,c) < 0)
        { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",Prog_Name,oname);
          Clean_Exit(1);
        }
    }

  //  Set # of k-mers into output file prolog

  data->tsize = anum;

  anum /= PMER_WORD;
  lseek(afile,sizeof(int),SEEK_SET);
  if (write(afile,&anum,sizeof(int64)) < 0)
    { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",Prog_Name,oname);
      Clean_Exit(1);
    }

  if (CLOCK)
    fprintf(stderr,"\r         \r");

  return (NULL);
}

  //  Top-Level

void Merge_Tables(char *path, char *root)
{ char       *fname;
  IO_block  **heap;
  IO_block   *io;
  uint8      *blocks;

#ifndef DEBUG_MERGE
  THREAD    threads[NTHREADS];
#endif
  Track_Arg parmk[NTHREADS];
  int       p, f, t, n;

  if (VERBOSE)
    { fprintf(stderr,"\nPhase 3 (-t option): Merging K-mer Table Parts\n");
      fflush(stderr);
    }

  fname = Malloc(2*(strlen(SORT_PATH) + strlen(path) + strlen(root)) + 100,"File name buffer");
  if (fname == NULL)
    Clean_Exit(1);

  //  Allocate all working data structures
 
  BUFLEN_UINT8 = SORT_MEMORY/((NPARTS+1)*NTHREADS);
  if (BUFLEN_UINT8 > 0x7fffffffll)
    BUFLEN_UINT8 = 0x7ffffff8ll;

  heap   = (IO_block **) Malloc(sizeof(IO_block *)*(NPARTS+1)*NTHREADS,"Allocating heap");
  io     = (IO_block *) Malloc(sizeof(IO_block)*(NPARTS+1)*NTHREADS,"Allocating IO buffers");
  blocks = (uint8 *) Malloc(BUFLEN_UINT8*(NPARTS+1)*NTHREADS,"Allocating IO buffers");
  if (heap == NULL || io == NULL || blocks == NULL)
    Clean_Exit(1);

  //  Open all input files

  p = NTHREADS;
  for (t = 0; t < NTHREADS; t++)
    for (n = 0; n < NPARTS; n++)
      { sprintf(fname,"%s/%s.%d.L%d",SORT_PATH,root,n,t);
        f = open(fname,O_RDONLY);
        if (f == -1)
          { fprintf(stderr,"\n%s: Cannot open external file %s in %s\n",
                           Prog_Name,fname,SORT_PATH);
            Clean_Exit(1);
          }
        io[p].block  = blocks + p*BUFLEN_UINT8;
        io[p].stream = f;
        p += 1;
      }

  //  Allocate and initialize prefix index

#ifdef DEVELOPER
  read(io[NTHREADS+NPARTS-1].stream,&IDX_BYTES,sizeof(int));
#endif
  PMER_WORD = TMER_WORD - IDX_BYTES;
  pidxlen   = (1 << (8*IDX_BYTES));
  pindex    = (int64 *) Malloc(sizeof(int64)*pidxlen,"Allocating index table");
  if (pindex == NULL)
    Clean_Exit(1);
  bzero(pindex,sizeof(int64)*pidxlen);

  //  Get size in bytes of table file thread 0 will produce for the clock

  if (VERBOSE)
    { struct stat info;

      totin = 0;
      for (p = NPARTS+NTHREADS-1; p >= NTHREADS; p--)
        { fstat(io[p].stream,&info);
          totin += info.st_size;
        }
#ifdef DEVELOPER
      totin = ((totin-sizeof(int))/TMER_WORD)*PMER_WORD;
#else
      totin = (totin/TMER_WORD)*PMER_WORD;
#endif
    }

  //  Remove previous table result if any

  sprintf(fname,"rm -f %s/%s.ktab %s/.%s.ktab.*",path,root,path,root);
  system(fname);

  //  Setup thread params and open output file

  for (t = 0; t < NTHREADS; t++)
    { sprintf(fname,"%s/.%s.ktab.%d",path,root,t+1);
      f = open(fname,O_CREAT|O_TRUNC|O_WRONLY,S_IRWXU);
      if (f == -1)
        { fprintf(stderr,"\n%s: Cannot open external file %s for writing\n",Prog_Name,fname);
          Clean_Exit(1);
        }
      io[t].block  = blocks + t*BUFLEN_UINT8;
      io[t].stream = f;

      parmk[t].out   = io + t;
      parmk[t].heap  = heap + t*(NPARTS+1);
      parmk[t].in    = io + t*NPARTS + NTHREADS;
      parmk[t].id    = t;
      parmk[t].oname = Strdup(fname,"Allocating stream name");
      if (parmk[t].oname == NULL)
        Clean_Exit(1);
    }

  //  In parallel merge part-files for each thread

#ifdef DEBUG_MERGE
  for (t = 0; t < NTHREADS; t++)
    merge_table_thread(parmk+t);
#else
  for (t = 1; t < NTHREADS; t++)
    pthread_create(threads+t,NULL,merge_table_thread,parmk+t);
  merge_table_thread(parmk);
  for (t = 1; t < NTHREADS; t++)
    pthread_join(threads[t],NULL);
#endif

  //  Close input files and if user-mode then remove them

  p = 0;
  for (t = 0; t < NTHREADS; t++)
    { for (n = 0; n <= NPARTS; n++)
        close(io[p++].stream);
      free(parmk[t].oname);
    }

#ifndef DEVELOPER
  for (p = 0; p < NPARTS; p++)
    for (t = 0; t < NTHREADS; t++)
      { sprintf(fname,"%s/%s.%d.L%d",SORT_PATH,root,p,t);
        unlink(fname);
      }
#endif

  //  Turn index counts to index offsets and create stub file

  { int x;

    for (x = 1; x < pidxlen; x++)
      pindex[x] += pindex[x-1];

    sprintf(fname,"%s/%s.ktab",path,root);
    f = open(fname,O_CREAT|O_TRUNC|O_WRONLY,S_IRWXU);
    if (f == -1)
      { fprintf(stderr,"\n%s: Cannot open external file %s for writing\n",Prog_Name,fname);
        Clean_Exit(1);
      }
    write(f,&KMER,sizeof(int));
    write(f,&NTHREADS,sizeof(int));
    write(f,&DO_TABLE,sizeof(int));
    write(f,&IDX_BYTES,sizeof(int));
    if (write(f,pindex,sizeof(int64)*pidxlen) < 0)
      { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",Prog_Name,fname);
        Clean_Exit(1);
      }

    close(f);
  }

  if (VERBOSE)
    { int64 tsize;

      tsize = 0;
      for (t = 0; t < NTHREADS; t++)
        tsize += parmk[t].tsize;

      fprintf(stderr,"  There are ");
      Print_Number(tsize/PMER_WORD,0,stderr);
      fprintf(stderr," %d-mers that occur %d-or-more times\n",KMER,DO_TABLE);

      tsize += 4*sizeof(int) + pidxlen*sizeof(int64) + NTHREADS*(sizeof(int)+sizeof(int64));

      if (tsize >= 5.e8)
        fprintf(stderr,"\n  The table occupies %.2f GB\n",tsize/1.e9);
      else if (tsize >= 5.e5)
        fprintf(stderr,"\n  The table occupies %.2f MB\n",tsize/1.e6);
      else
        fprintf(stderr,"\n  The table occupies %.2f KB\n",tsize/1.e3);
      fflush(stderr);
    }

  free(pindex);
  free(blocks);
  free(io);
  free(heap);
  free(fname);
}
