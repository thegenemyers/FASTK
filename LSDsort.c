/*******************************************************************************************
 *
 *  Fast threaded lexical sort routine.  Can be compiled to accommodate any element size
 *     (set WORD_SIZE), and makes only n+1 passes to sort n radix bytes.  The radix order
 *     for the bytes of an element may be sorted in any order as listed in the array bytes
 *     (that is -1 terminated).
 *
 *  Author :  Gene Myers
 *  First  :  May 2018
 *
 ********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>

#include "libfastk.h"
#include "FastK.h"

typedef unsigned char uint8;
typedef long long     int64;

#undef  TEST_SORT
#undef  CHECK_SORT

static int    RSIZE;    //  Span between records

//  Global variables for every "lex_thread"

static int      LEX_byte;   //  Current byte to sort on
static int      LEX_next;   //  Next byte to sort on (if >= 0)
static int64    LEX_zdiv;   //  Size of thread segments (in bytes)
static uint8   *LEX_src;    //  Source data goes to ...
static uint8   *LEX_trg;    //  Target data

//  Thread control record

typedef struct
  { int64  beg;           //  Sort [beg,end) of LEX_src
    int64  end;
    int    check[256];    //  Not all of bucket will go to the same thread in the next cycle?
    int    next[256];     //  Thread assignment for next cycle (updated if check true)
    int64  thresh[256];   //  If check then multiple of LEX_zdiv to check for thread assignment
    int64  tptr[256];     //  Finger for each 8-bit value
    int64 *sptr;          //  Conceptually [256][NTHREADS].  At end of sorting pass
  } Lex_Arg;              //    sprtr[b][n] = # of occurences of value b in rangd of
                          //    thread n for the *next* pass

//  Threaded sorting pass

static void *lex_thread(void *arg)
{ Lex_Arg *data   = (Lex_Arg *) arg;
  int64   *sptr   = data->sptr;
  int64   *tptr   = data->tptr;
  uint8   *src    = LEX_src;
  uint8   *dig    = LEX_src + LEX_byte;
  uint8   *nig    = LEX_src + LEX_next;
  uint8   *trg    = LEX_trg;
  int64    zdiv   = LEX_zdiv;
  int     *check  = data->check;
  int     *next   = data->next;
  int64   *thresh = data->thresh;

  int64       i, n, x;
  uint8       d;

  n = data->end;
  if (LEX_next < 0)
    for (i = data->beg; i < n; i += RSIZE)
      { d = dig[i];
        x = tptr[d];
        tptr[d] += RSIZE;
        memcpy(trg+x,src+i,RSIZE);
      }
  else
    for (i = data->beg; i < n; i += RSIZE)
      { d = dig[i];
        x = tptr[d];
        tptr[d] += RSIZE;
        memcpy(trg+x,src+i,RSIZE);
        if (check[d])
          { if (x >= thresh[d])
              { next[d]   += 0x100;
                thresh[d] += zdiv;
              }
          }
        sptr[next[d] | nig[i]] += 1;
      }
  return  (NULL);
}

//  Threaded sort initiation pass: count bucket sizes

static void *lexbeg_thread(void *arg)
{ Lex_Arg    *data  = (Lex_Arg *) arg;
  int64      *tptr  = data->tptr;
  uint8      *dig   = LEX_src + LEX_byte;

  int64       i, n;

  n = data->end;
  for (i = data->beg; i < n; i += RSIZE)
    tptr[dig[i]] += 1;
  return (NULL);
}

//  Radix sort the indicated "bytes" of src, using array trg as the secondary array
//    The arrays contains len elements each of "size" bytes.
//    Return a pointer to the array containing the final result.

void *LSD_Sort(int64 nelem, void *src, void *trg, int rsize, int *bytes)
{ pthread_t *threads;
  Lex_Arg   *parmx;   //  Thread control record for sorting

  uint8   *xch;
  int64    x, y, asize;
  int      i, j, z, b;

  asize = nelem*rsize;
  RSIZE = rsize;

  LEX_zdiv = ((nelem-1)/NTHREADS + 1)*RSIZE;
  LEX_src  = (uint8 *) src;
  LEX_trg  = (uint8 *) trg;

  parmx   = Malloc(sizeof(Lex_Arg)*NTHREADS,"LSD sort vectors");
  threads = Malloc(sizeof(pthread_t)*NTHREADS,"LSD sort vectors");
  parmx[0].sptr = Malloc(sizeof(int64)*256*NTHREADS*NTHREADS,"LSD sort vectors");
  if (parmx == NULL || threads == NULL || parmx[0].sptr == NULL)
    exit (1);

  for (i = 1; i < NTHREADS; i++)
    parmx[i].sptr = parmx[i-1].sptr + NTHREADS*256;

  //  For each requested byte b in order, radix sort

  for (b = 0; bytes[b] >= 0; b++)
    { LEX_byte  = bytes[b];
      LEX_next  = bytes[b+1];

      //  Setup beg, end, and zero tptr counters

      x = 0;
      for (i = 0; i < NTHREADS; i++)
        { parmx[i].beg = x;
          x = LEX_zdiv*(i+1);
          if (x > asize)
            x = asize;
          parmx[i].end = x;
          for (j = 0; j < 256; j++)
            parmx[i].tptr[j] = 0;
        }
      parmx[NTHREADS-1].end = asize;

      //  If first pass, then explicitly sweep to get tptr counts
      //    otherwise accumulate from sptr counts of last sweep

      if (b == 0)
        { for (i = 1; i < NTHREADS; i++)
            pthread_create(threads+i,NULL,lexbeg_thread,parmx+i);
          lexbeg_thread(parmx);
          for (i = 1; i < NTHREADS; i++)
            pthread_join(threads[i],NULL);
        }
      else
        { int64 *pxt, *pxs;

          for (i = 0; i < NTHREADS; i++)
            { pxt = parmx[i].tptr;
              for (z = 0; z < NTHREADS; z++)
                { pxs = parmx[z].sptr + (i<<8);
                  for (j = 0; j < 256; j++)
                    pxt[j] += pxs[j];
                }
            }
        }

      //   Zero sptr array counters in preparation of pass

      for (i = 0; i < NTHREADS; i++)
        for (z = (NTHREADS<<8)-1; z >= 0; z--)
          parmx[i].sptr[z] = 0;

      //  Convert tptr from counts to fingers, and determine thead assignment arrays
      //    to avoid a division in the inner most loop

      { int64 thr;
        int   nxt;

        thr = LEX_zdiv;
        nxt = 0;
        x = 0;
        for (j = 0; j < 256; j++)
          for (i = 0; i < NTHREADS; i++)
            { y = parmx[i].tptr[j]*RSIZE;
              parmx[i].tptr[j] = x;
              x += y;
              parmx[i].next[j] = nxt;
              if (x < thr)
                parmx[i].check[j] = 0;
              else
                { parmx[i].check[j]  = 1;
                  parmx[i].thresh[j] = thr;
                  while (x >= thr)
                    { thr += LEX_zdiv;
                      nxt += 0x100;
                    }
                }
            }
      }

      //  Threaded pass

      for (i = 1; i < NTHREADS; i++)
        pthread_create(threads+i,NULL,lex_thread,parmx+i);
      lex_thread(parmx);
      for (i = 1; i < NTHREADS; i++)
        pthread_join(threads[i],NULL);

      xch     = LEX_src;
      LEX_src = LEX_trg;
      LEX_trg = xch;

#ifdef TEST_SORT
      { int64  c;
        uint8 *psort = LEX_src-RSIZE;

        printf("\nLSORT %d\n",LEX_byte);
        for (c = 0; c < 1000*RSIZE; c += RSIZE)
          { printf(" %4lld: ",c/RSIZE);
            for (j = 0; j < RSIZE; j++)
              printf(" %02x",LEX_src[c+j]);
            printf("\n");
          }
      }
#endif

#ifdef CHECK_SORT
      { int64  c;
        uint8 *psort = LEX_src-RSIZE;

        printf("\nLSORT %d\n",LEX_byte);
        for (c = RSIZE; c < asize; c += RSIZE)
          { for (j = b; j >= 0; j--)
              if (LEX_src[c+bytes[j]] > psort[c+bytes[j]])
                break;
              else if (LEX_src[c+bytes[j]] < psort[c+bytes[j]])
                { printf("  Order: %lld",c/RSIZE);
                  for (x = 2; x <= LEX_byte; x++)
                    printf(" %02x",psort[c+x]);
                  printf(" vs");
                  for (x = 2; x <= LEX_byte; x++)
                    printf(" %02x",LEX_src[c+x]);
                  printf("\n");
                  break;
                }
          }
      }
#endif
    }

  free(parmx[0].sptr);
  free(threads);
  free(parmx);

  return ((void *) LEX_src);
}
