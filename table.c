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

#include "gene_core.h"
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
    int        id;
  } Track_Arg;

static int64 totin;  //  Total kmer record for theard 0

static void *merge_table_thread(void *arg)
{ Track_Arg   *data  = (Track_Arg *) arg;
  IO_block   **heap  = data->heap;
  IO_block    *in    = data->in;
  IO_block    *out   = data->out;

  int    afile;
  int64  anum;
  uint8 *abuf, *aptr, *atop;
  int    hsize;
  int    p;

  int64  pct1, nextin;
  int    CLOCK;

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
  write(afile,&anum,sizeof(int64));

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

      //  Flush output buffer if needed

      if (aptr + TMER_WORD > atop)
        { int c = aptr-abuf;

#ifdef DEBUG
          printf("Write %d bytes\n",c);
#endif
          anum += c;
          write(afile,abuf,c);
          aptr = abuf;
          if (CLOCK && anum > nextin)
            { fprintf(stderr,"\r  %3d%%",(int) ((100.*anum)/totin));
              fflush(stderr);
              nextin = anum+pct1;
            }
        }

      //  Append k-mer to output

      mycpy(aptr,sptr,TMER_WORD);
      aptr += TMER_WORD;
      sptr += TMER_WORD;

#ifdef DEBUG
      printf(" %3d: ",p);
      print_kmer(aptr-TMER_WORD,KMER);
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
      write(afile,abuf,c);
    }

  //  Set # of k-mers into output file prolog

  anum /= TMER_WORD;
  lseek(afile,sizeof(int),SEEK_SET);
  write(afile,&anum,sizeof(int64));

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

  fname = Malloc(strlen(SORT_PATH) + strlen(path) + strlen(root) + 100,"File name buffer");
  if (fname == NULL)
    exit (1);

  //  Get rid of any previous results for this DB in this directory with this KMER

  sprintf(fname,"rm -f %s/%s.K%d.T*",path,root,KMER);
  system(fname);

  //  Allocate all working data structures
 
  BUFLEN_UINT8 = SORT_MEMORY/((NPARTS+1)*NTHREADS);
  if (BUFLEN_UINT8 > 0x7fffffffll)
    BUFLEN_UINT8 = 0x7ffffff8ll;

  heap   = (IO_block **) Malloc(sizeof(IO_block *)*(NPARTS+1)*NTHREADS,"Allocating heap");
  io     = (IO_block *) Malloc(sizeof(IO_block)*(NPARTS+1)*NTHREADS,"Allocating IO buffers");
  blocks = (uint8 *) Malloc(BUFLEN_UINT8*(NPARTS+1)*NTHREADS,"Allocating IO buffers");
  if (heap == NULL || io == NULL || blocks == NULL)
    exit (1);

  //  Open all input files

  p = NTHREADS;
  for (t = 0; t < NTHREADS; t++)
    for (n = 0; n < NPARTS; n++)
      { sprintf(fname,"%s/%s.%d.L%d",SORT_PATH,root,n,t);
        f = open(fname,O_RDONLY);
        if (f == -1)
          { fprintf(stderr,"\n%s: Cannot open external file %s in %s\n",
                           Prog_Name,fname,SORT_PATH);
            exit (1);
          }
        io[p].block  = blocks + p*BUFLEN_UINT8;
        io[p].stream = f;
        p += 1;
      }

  if (VERBOSE)
    { struct stat info;

      totin = 0;
      for (p = NPARTS+NTHREADS-1; p >= NTHREADS; p--)
        { fstat(io[p].stream,&info);
          totin += info.st_size;
        }
    }

  //  Setup thread params and open output file

  sprintf(fname,"%s/%s.ktab",path,root);
  f = open(fname,O_CREAT|O_TRUNC|O_WRONLY,S_IRWXU);
  if (f == -1)
    { fprintf(stderr,"\n%s: Cannot open external file %s for writing\n",Prog_Name,fname);
      exit (1);
    }
  write(f,&KMER,sizeof(int));
  write(f,&NTHREADS,sizeof(int));
  close(f);

  for (t = 0; t < NTHREADS; t++)
    { sprintf(fname,"%s/.%s.ktab.%d",path,root,t+1);
      f = open(fname,O_CREAT|O_TRUNC|O_WRONLY,S_IRWXU);
      if (f == -1)
        { fprintf(stderr,"\n%s: Cannot open external file %s for writing\n",Prog_Name,fname);
          exit (1);
        }
      io[t].block  = blocks + t*BUFLEN_UINT8;
      io[t].stream = f;

      parmk[t].out   = io + t;
      parmk[t].heap  = heap + t*(NPARTS+1);
      parmk[t].in    = io + t*NPARTS + NTHREADS;
      parmk[t].id    = t;
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

  p = 0;
  for (t = 0; t < NTHREADS; t++)
    for (n = 0; n <= NPARTS; n++)
      close(io[p++].stream);

  //  If user-mode the remove input files

#ifndef DEVELOPER
  for (p = 0; p < NPARTS; p++)
    for (t = 0; t < NTHREADS; t++)
      { sprintf(fname,"%s/%s.%d.L%d",SORT_PATH,root,p,t);
        unlink(fname);
      }
#endif

  free(blocks);
  free(io);
  free(heap);
  free(fname);
}
