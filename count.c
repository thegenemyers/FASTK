/*********************************************************************************************\
 *
 *  Phase 2 of FastK:  For each set of the NTHREADS files for a each partition bucket produced
 *     by the split.c module, do the following:
 *       * Sort the super-mers.
 *       * Sort the k-mers from each super-mer weighted by the # of times that super-mer occurs.
 *       * During the 2nd sort accumulate the histogram of k-mer frequencies.
 *            Output this to <source_root>.K<kmer>.
 *       * if requesteda (-t) produce a table of all the k-mers with counts >= -t in NTHREADSs
 *            pieces in files SORT_PATH/<root>.<bucket>.L<thread>
 *       * if requested (-p) invert the first two sorts to produce a profile for every super-mer
 *            in the order in the source in files SORT_PATH/<root>.<bucket>.P<thread>.[0-3]
 *       * if requested (-h) print the histogram of k-mer frequencies.
 *
 *  Author:  Gene Myers
 *  Date  :  October, 2020
 *
 *********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <pthread.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "libfastk.h"
#include "FastK.h"

#undef  DEBUG_COMPRESS
#undef  DEBUG_SLIST
#undef  DEBUG_KLIST
#undef  DEBUG_CANONICAL
#undef  DEBUG_TABOUT
#undef  DEBUG_CLIST
#undef  DEBUG_CMERGE
#undef  EQUAL_MERGE
#undef  DEBUG_PLIST
#undef    SHOW_RUN
#undef  DEBUG_PWRITE

#define THREAD  pthread_t


/*******************************************************************************************
 *
 * static void *supermer_list_thread(Slist_Arg *arg)
 *     Each thread reads a file containing a part of the bit encoded super-mers for a given
 *     bucket and unpacks them into an array for sorting.  If -p each entry also contains
 *     it super-mer ordinal number in the input.
 *
 ********************************************************************************************/

  //  Unstuff ints and reads from a bit packed buffer

uint8 Comp[256];

#if defined(DEBUG_CANONICAL) || defined(DEBUG_KLIST) || defined(DEBUG_SLIST) \
      || defined(DEBUG_TABOUT) || defined(DEBUG_CMERGE) || defined(EQUAL_MERGE)

static char DNA[4] = { 'a', 'c', 'g', 't' };

static char *fmer[256], _fmer[1280];

#endif

#if defined(DEBUG_SLIST) || defined(DEBUG_TABOUT) || defined(DEBUG_CMERGE) || defined(EQUAL_MERGE)

static void write_Ascii(uint8 *bytes, int len)
{ int i;

  for (i = 0; i < len; i += 4)
    printf("%s",fmer[*bytes++]);
}

#endif

static inline IO_UTYPE *Unstuff_Int(int64 *val, int nbits, IO_UTYPE mask, IO_UTYPE *buf, int *bitp)
{ int rem;

  rem = *bitp;
  if (nbits < rem)
    { rem -= nbits;
      *val = ((*buf) >> rem) & mask; 
#ifdef DEBUG_COMPRESS
      printf(" b = %2d/%016llx -%d-> v = %lld\n",rem,*buf,nbits,*val);
#endif
    }
  else if (nbits == rem)
    { *val = (*buf++) & mask;
      rem  = IO_UBITS;
#ifdef DEBUG_COMPRESS
      printf(" b = %2d/%016llx -%d-> v = %lld\n",rem,buf[-1],nbits,*val);
#endif
    }
  else
    { IO_UTYPE x = (*buf++) << (nbits-rem);
      rem += IO_UBITS - nbits;
      *val = (x | (*buf >> rem)) & mask;
#ifdef DEBUG_COMPRESS
      printf(" b = %2d/%016llx|%016llx -%d-> v = %lld\n",rem,buf[-1],*buf,nbits,*val);
#endif
    }
  *bitp = rem;
  return (buf);
}

static inline IO_UTYPE *Unstuff_Code(uint8 *s, int len, IO_UTYPE *buf, int *bitp)
{ int      i, rem;
  IO_UTYPE x;

  rem = *bitp;
  len <<= 1;
  for (i = 8; i < len; i += 8)
    if (rem > 8)
      { rem -= 8;
        *s++ = ((*buf) >> rem) & 0xff;
      }
    else if (rem == 8)
      { *s++ = *buf++ & 0xff;
        rem  = IO_UBITS;
      }
    else
      { x = (*buf++) << (8-rem);
        rem += IO_UBITS - 8;
        *s++ = (x | (*buf >> rem)) & 0xff;
      }
  i = len - (i-8); 
  if (i > 0)
    { if (rem > i)
        { rem -= i;
          x = ((*buf) >> rem);
        }
      else if (rem == i)
        { x = *buf++;
          rem  = IO_UBITS;
        }
      else
        { x = (*buf++) << (i-rem);
          rem += IO_UBITS - i;
          x |= ((*buf) >> rem);
        }
      *s++ = (x << (8-i)) & 0xff;
    }
  *bitp = rem;
  return (buf);
}

  // Read smers stuffed super-mers from tfile producing kmer list, loading at fill

static int  Fixed_Reload[IO_UBITS+1];
static int  Runer_Reload[IO_UBITS+1];
static int *Super_Reload[IO_UBITS+1];

typedef struct
  { int     tfile;      //  Bit compressed super-mer input streaam for thread
    int64   nmers;      //  # of super-mers in the input
    int64   nbase;      //  if DO_PROFILE then start of indices for this thread
    int64   nidxs;      //  total number of super-mers for this thread slice
    uint8  *fours[256]; //  finger for filling list sorted on first super-mer byte
  } Slist_Arg;

static void *supermer_list_thread(void *arg)
{ Slist_Arg   *data = (Slist_Arg *) arg;
  int64        nmers = data->nmers;
  int64        nbase = data->nbase;
  int          in    = data->tfile;
  uint8      **fours = data->fours;

  uint8 *fill, *prev;
  int64  f, *fb = &f;
  int64  rmask, rlim;
  int    N, i, rbits;

  int64 r = 0;
  int64 n = 0;
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
  uint8  *nb = ((uint8 *) &n) - 1;
  uint8  *rb = ((uint8 *) &r) - 1;
#else
  uint8  *nb = ((uint8 *) &n) + (sizeof(int64)-SLEN_BYTES);
  uint8  *rb = ((uint8 *) &r) + (sizeof(int64)-RUN_BYTES);
#endif

  IO_UTYPE iobuf[IO_BUF_LEN], *ioend, *ptr;
  int    bit, clen; 
  int64  k;

  read(in,iobuf,IO_UBYTES*IO_BUF_LEN);
  ioend = iobuf + IO_BUF_LEN;
  ptr   = iobuf;
  bit   = IO_UBITS;
  rbits = 17;
  rlim  = 0x10000ll;
  rmask = 0x1ffffll;

#ifdef DEBUG_SLIST
  printf("Index at %d(%llx)\n",rbits,rmask);
#endif

  prev = fours[0];
  for (k = 0; k < nmers; k++)
    { while (1)
        { if (ptr + Fixed_Reload[bit] >= ioend)
            { int res = 0;
              while (ptr < ioend)
                iobuf[res++] = *ptr++;
              read(in,iobuf+res,IO_UBYTES*(IO_BUF_LEN-res));
              ptr = iobuf;
            }
          ptr = Unstuff_Int(&n,SLEN_BITS,SLEN_BIT_MASK,ptr,&bit);

          if (n < MAX_SUPER)
            break;

          *prev |= 0x80;
#ifdef DEBUG_SLIST
          printf("+++\n");
#endif
        }

      if (ptr + Super_Reload[bit][n] >= ioend)
        { int res = 0;
          while (ptr < ioend)
            iobuf[res++] = *ptr++;
          read(in,iobuf+res,IO_UBYTES*(IO_BUF_LEN-res));
          ptr = iobuf;
        }

      N = n+KMER;

      ptr = Unstuff_Int(fb,8,0xffllu,ptr,&bit);
      
      fill = fours[f];
      fours[f] += SMER_WORD;
      *fill++ = 0;

      ptr = Unstuff_Code(fill,N-4,ptr,&bit);
      clen = ((N-1)>>2);
      fill += clen;

      for (i = clen+1; i < SMER_BYTES; i++)
        *fill++ = 0;
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
      for (i = SLEN_BYTES; i > 0; i--)
        *fill++ = nb[i];
#else
      for (i = 0; i < SLEN_BYTES; i++)
        *fill++ = nb[i];
#endif

      if (DO_PROFILE)
        { while (1)
            { if (ptr + Runer_Reload[bit] >= ioend)
                { int res = 0;
                  while (ptr < ioend)
                    iobuf[res++] = *ptr++;
                  read(in,iobuf+res,IO_UBYTES*(IO_BUF_LEN-res));
                  ptr = iobuf;
                }
              ptr = Unstuff_Int(&r,rbits,rmask,ptr,&bit);
              if (r < rlim)
                 break;
              rbits  += 1;
              rlim  <<= 1;
              rmask   = (rmask << 1) + 1;
#ifdef DEBUG_SLIST
              printf("Index bumped to %d(%llx) r = %llu\n",rbits,rmask,r);
#endif
            }
          r += nbase;
          prev = fill;
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
          for (i = RUN_BYTES; i > 0; i--)
            *fill++ = rb[i];
#else
          for (i = 0; i < RUN_BYTES; i++)
            *fill++ = rb[i];
#endif
        }

#ifdef DEBUG_SLIST
      if (DO_PROFILE)
        printf(" %5lld: ",r);
      printf("[%02llx] %s ",f,fmer[f]);
      write_Ascii(fill-(SMER_WORD-1),N-4);
      printf(" : %d\n",N);
      printf("            ");
      for (i = -SMER_WORD; i < 0; i++)
        printf(" %02x",fill[i]);
      printf("\n");
#endif
    }

  if (ptr + Fixed_Reload[bit] >= ioend)
    { int res = 0;
      while (ptr < ioend)
        iobuf[res++] = *ptr++;
      read(in,iobuf+res,IO_UBYTES*(IO_BUF_LEN-res));
      ptr = iobuf;
    }
  ptr = Unstuff_Int(&n,SLEN_BITS,SLEN_BIT_MASK,ptr,&bit);
  if (n >= MAX_SUPER)
    { *prev |= 0x80;
#ifdef DEBUG_SLIST
      printf("+++\n");
#endif
    }

  return (NULL);
}


/*******************************************************************************************
 *
 * static void *kmer_list_thread(Klist_Arg *arg)
 *     Each thread takes the now sorted super-mers, counts the # of each unique super-mer,
 *     and places the canonical k-mers of each with weights in an array for sorting.
 *     If the -p option is set then each k-mer entry also has the ordinal number of the
 *     k-mer in order of generation.
 *
 ********************************************************************************************/

typedef struct
  { uint8    *sort;
    int64    *parts;
    int       beg;
    int       end;
    int64     off;
    uint8    *fours[256];
    int64     kidx;
    int64     overflow;
  } Klist_Arg;

static int kclip[4] = { 0xff, 0xc0, 0xf0, 0xfc };

static void *kmer_list_thread(void *arg)
{ Klist_Arg   *data   = (Klist_Arg *) arg;
  int          beg    = data->beg;
  int          end    = data->end;
  int64       *part   = data->parts;
  uint8      **fours  = data->fours;
  int64        overflow;

  int       KMp3   = KMER+3;
  int       KRS    = ((KMER-1)&0x3)<<1;
  int       KMd2   = (KMER_BYTES+1)>>1;
  int       fptr[SMER_BYTES+1], rptr[SMER_BYTES+1];
  int      *FPT    = fptr-1;
  int      *RPT    = rptr-(KMER_BYTES-1);
  int       KCLIP  = kclip[KMER&0x3];

  uint8    *sptr, *send, *lptr;
  int       x, ct, sbytes;
  uint8    *asp, *fill;

#ifdef DEBUG_KLIST
  int64  ridx = 0;
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
  uint8  *rbx = ((uint8 *) &ridx) - 1;
#else
  uint8  *rbx = ((uint8 *) &ridx) + (sizeof(int64)-RUN_BYTES);
#endif
#endif

  int       i, o;
  int      *f, *r;
  int       fb, rb;
  int       fs, rs;
  int       kb, hb;
  int       kf, hf;

  int   sln = 0;
  int64 idx = 0;
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
  uint8  *sb = ((uint8 *) &sln) - 1;
  uint8  *ib = ((uint8 *) &idx) - 1;
#else
  uint8  *sb = ((uint8 *) &sln) + (sizeof(int)-SLEN_BYTES);
  uint8  *ib = ((uint8 *) &idx) + (sizeof(int64)-KMAX_BYTES);
#endif

  overflow = 0;

  idx  = data->kidx;
  sptr = data->sort + data->off;
  for (x = beg; x < end; x++)
    for (send = sptr + part[x]; sptr < send; sptr = lptr)

      { asp = sptr + SMER_BYTES;
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
        for (i = SLEN_BYTES; i > 0; i--)
          sb[i] = *asp++;
#else
        for (i = 0; i < SLEN_BYTES; i++)
          sb[i] = *asp++;
#endif
 
#ifdef DEBUG_KLIST
        printf("Run ids:");
        lptr = sptr;
        do
          { lptr += SMER_WORD;
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
            for (i = RUN_BYTES; i > 0; i--)
              rbx[i] = lptr[-i];
            rbx[RUN_BYTES] &= 0x7ff;
#else
            for (i = 0; i < RUN_BYTES; i++)
              rbx[i] = lptr[i-RUN_BYTES];
            rbx[0] &= 0x7ff;
#endif
            printf(" %c%lld",lptr[-RUN_BYTES]&0x8?'+':' ',ridx);
          }
        while (*lptr == 0);
        printf("\n");
#endif

        lptr = sptr + SMER_WORD;
        ct   = 1;
        while (*lptr == 0)
          { lptr += SMER_WORD;
            ct   += 1;
          }

        *sptr  = x;
        sbytes = (KMp3 + sln) >> 2;

#ifdef DEBUG_KLIST
        printf("%02x/%lld: l=%d c=%d\n  ",x,idx,sln+1,ct);
        for (i = 0; i < sbytes; i++)
          printf(" %s",fmer[sptr[i]]);
        printf("\n");
        fflush(stdout);
#endif

        for (o = sbytes, i = 0; i <= sbytes; i++, o--)
          { fptr[i] = (sptr[i] << 8) | sptr[i+1];
            rptr[i] = (Comp[sptr[o]] << 8) | Comp[sptr[o-1]];
          }

#ifdef DEBUG_CANONICAL
        printf("   F = ");
        for (i = 0; i < sbytes; i++)
          printf(" %s%s",fmer[fptr[i]>>8],fmer[fptr[i]&0xff]);
        printf("\n   R = ");
        for (i = 0; i < sbytes; i++)
          printf(" %s%s",fmer[rptr[i]>>8],fmer[rptr[i]&0xff]);
        printf("\n");
        fflush(stdout);
#endif

        if (ct >= 0x8000)
          { overflow += ((int64) (ct-0x7fff))*(sln+1);
            ct = 0x7fff;
          }

        f  = FPT;
        fs = 2;
        fb = 0;
        r  = RPT + sbytes;
        rs = KRS;
        rb = *r;
        for (o = 0; o <= sln; o++)
          { fs -= 2;
            rs += 2;
            if (fs == 0)
              { fb = *++f;
                fs = 8;
              }
            kb = kf = (fb >> fs) & 0xff;
            if (rs == 8)
              { rb = *--r;
                rs = 0;
              }
            hb = hf = (rb >> rs) & 0xff;
#ifdef DEBUG_CANONICAL
            printf("   + %d / %s%s / %s\n",fs, fmer[fb>>8], fmer[fb&0xff], fmer[kb]);
            printf("   - %d / %s%s / %s\n",rs, fmer[rb>>8], fmer[rb&0xff], fmer[hb]);
            fflush(stdout);
#endif
            if (kb == hb)
              { for (i = 1; i < KMd2; i++)
                  { kb = (f[i] >> fs) & 0xff;
                    hb = (r[i] >> rs) & 0xff;
#ifdef DEBUG_CANONICAL
                    printf("      @ %d: %s vs %s\n",i,fmer[kb],fmer[hb]);
                    fflush(stdout);
#endif
                    if (kb != hb)
                      break;
                  }
              }

            if (kb < hb)
              { fill = fours[kf];
                fours[kf] = fill+KMER_WORD;
                *fill++ = 0;
                for (i = 1; i < KMER_BYTES; i++)
                  *fill++ = (f[i] >> fs) & 0xff;
              }
            else
              { fill = fours[hf];
                fours[hf] = fill+KMER_WORD;
                *fill++ = 0;
                for (i = 1; i < KMER_BYTES; i++)
                  *fill++ = (r[i] >> rs) & 0xff;
              }
            fill[-1] &= KCLIP;
            *((uint16 *) fill) = ct;
            fill += 2;
            if (DO_PROFILE)
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
              { for (i = KMAX_BYTES; i > 0; i--)
                  *fill++ = ib[i];
#else
              { for (i = 0; i < KMAX_BYTES; i++)
                  *fill++ = ib[i];
#endif
                idx  += 1;
              }

#ifdef DEBUG_KLIST
            printf("    [%02x/%s]:",kf < hf ? kf : hf,fmer[kf < hf ? kf : hf]);
            for (i = 1-KMER_WORD; i < KMER_BYTES-KMER_WORD; i++)
              printf(" %s",fmer[fill[i]]);
            for (i = KMER_BYTES-KMER_WORD; i < 0; i++)
              printf(" %02x",fill[i]);
            printf("\n");
            fflush(stdout);
#endif
          }

        *sptr = 1;
      }

  data->overflow = overflow;

  return (NULL);
}


/*******************************************************************************************
 *
 * static void *table_write_thread(Twrite_Arg *arg)
 *     Each thread takes the now sorted weighted k-mers, adds together the wgt.s of each
 *     unique k-mer, and outputs to a "L" file each k-mer and count pair that has weight
 *     greater thaan DO_TABLE (the -t option value, called only if set).
 *
 ********************************************************************************************/
typedef struct
  { uint8    *sort;
    int64    *parts;
    int       beg;
    int       end;
    int64     off;
    int       kfile;
    char     *kname;
    int64     tmers;
  } Twrite_Arg;

static void *table_write_thread(void *arg)
{ Twrite_Arg  *data   = (Twrite_Arg *) arg;
  int          beg    = data->beg;
  int          end    = data->end;
  int64       *part   = data->parts;
  int          kfile  = data->kfile;

  uint8  bufr[0x10000];
  uint8 *fill, *bend;
  int    x, ct;
  uint8 *kptr, *lptr, *kend;
  int64  tmer;

  fill = bufr;
  bend = bufr + (0x10000 - TMER_WORD);
  tmer = 0;

  kptr = data->sort + data->off;
  for (x = beg; x < end; x++)
    for (kend = kptr + part[x]; kptr < kend; kptr = lptr)
      { ct = *((uint16 *) (kptr+KMER_BYTES));
        lptr = kptr+KMER_WORD;
        while (*lptr == 0)
          lptr += KMER_WORD;
        if (ct >= DO_TABLE)
          { if (fill >= bend)
              { if (write(kfile,bufr,fill-bufr) < 0)
                  { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",
                                   Prog_Name,data->kname);
                    Clean_Exit(1);
                  }
                fill = bufr;
              }

            *kptr = x;
            memcpy(fill,kptr,TMER_WORD);
#ifdef DEBUG_TABOUT
            write_Ascii(kptr,KMER);
            printf(" %hd\n",*((uint16 *) (kptr+KMER_BYTES)));
#endif
            fill += TMER_WORD;
            tmer += 1;
          }
      }
  if (fill > bufr)
    if (write(kfile,bufr,fill-bufr) < 0)
      { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",Prog_Name,data->kname);
        Clean_Exit(1);
      }

  data->tmers = tmer;
  return (NULL);
}


/*******************************************************************************************
 *
 * static void *cmer_list_thread(Clist_Arg *arg)
 *     Each thread takes the now sorted weighted k-mers entries and reduces each to
 *     an entry containg the count and k-mer ordinal index, in preparation for the
 *     first of two "inverting" sorts to realize the -p option.
 *
 ********************************************************************************************/

typedef struct
  { uint8    *sort;
    int64    *parts;
    int       beg;
    int       end;
    int64     off;
    uint8    *fours[256];

    Kmer_Stream *stm; //  used only by cmer_merge_thread
  } Clist_Arg;

static void *cmer_list_thread(void *arg)
{ Clist_Arg   *data   = (Clist_Arg *) arg;
  int          beg    = data->beg;
  int          end    = data->end;
  int64       *part   = data->parts;
  uint8      **fours  = data->fours;

  int KM1 = KMER_WORD-1;

  int    x, d, k;
  uint8 *kptr, *kend;
  uint8 *fill;

  kptr = data->sort + data->off;
  for (x = beg; x < end; x++)
    for (kend = kptr + part[x]; kptr < kend; kptr += KMER_WORD)
      { d = kptr[KM1];
        fill = fours[d];
        for (k = KMER_BYTES; k < KM1; k++)
          *fill++ = kptr[k];
        fours[d] = fill;
      } 

  return (NULL);
}


/*******************************************************************************************
 *
 * static void *cmer_merge_thread(Clist_Arg *arg)
 *     Each thread takes the now sorted weighted k-mers entries and reduces each to
 *     an entry containg the count and k-mer ordinal index, in preparation for the
 *     first of two "inverting" sorts to realize the -p option.
 *
 ********************************************************************************************/

static void *cmer_merge_thread(void *arg)
{ Clist_Arg   *data   = (Clist_Arg *) arg;
  int          beg    = data->beg;
  int          end    = data->end;
  int64       *part   = data->parts;
  uint8      **fours  = data->fours;
  Kmer_Stream *S      = data->stm;
  int          hbyte  = S->hbyte;

  int KM1 = KMER_WORD-1;

  int    x, d, k, c;
  uint8 *kptr, *lptr, *kend;
  uint8 *fill;
  int    ct;

#ifdef EQUAL_MERGE
  int    skip = 0;
#endif

  if (beg > 0)
    S = Clone_Kmer_Stream(S);

  kptr = data->sort + data->off;
  *kptr = beg;
  GoTo_Kmer_Entry(S,kptr);
#ifdef DEBUG_CMERGE
  printf("\n\nStart: x = %s  k = ",fmer[beg]);
  write_Ascii(kptr,KMER);
  printf(" e = %s",fmer[S->cpre]);
  write_Ascii(S->csuf,KMER-4);
  printf(" %lld\n",S->cidx);
#endif
  for (x = beg; x < end; x++)

    { while (S->cpre < x)
        Next_Kmer_Entry(S);

      for (kend = kptr + part[x]; kptr < kend; )
        { lptr = kptr+KMER_WORD;
          while (*lptr == 0)
            lptr += KMER_WORD;

          while (1)
            { if (S->cpre > x)
                { c = 5000;
                  break;
                }
              c = memcmp(S->csuf,kptr+1,hbyte);
              if (c >= 0)
                break;
              Next_Kmer_Entry(S);
#ifdef EQUAL_MERGE
              printf("Skip?\n"); fflush(stdout);
              skip = 1;
#endif
            }
#ifdef EQUAL_MERGE
          if (c != 0 || skip)
            { int u;

              printf("\nx = %02x  kptr = %ld start = %ld\n",
                     x,kptr-data->sort,(kend-part[x])-data->sort);
              write_Ascii(kptr+1,KMER-4);
              printf("\n");
              write_Ascii(S->csuf,KMER-4);
              printf("\n");
              GoTo_Kmer_Index(S,S->cidx-3);
              kptr -= 3*KMER_WORD;
              for (u = 0; u < 7; u++)
                { printf("c = %4d k = ",c);
                  write_Ascii(kptr,KMER);
                  printf(" %d  e = %s",*((uint16 *) (kptr+KMER_BYTES)),fmer[S->cpre]);
                  write_Ascii(S->csuf,KMER-4);
                  printf("  %d %lld\n",*((uint16 *) (S->csuf+hbyte)),S->cidx);
                  kptr += KMER_WORD;
                  Next_Kmer_Entry(S);
                }
              kptr -= 4*KMER_WORD;
              GoTo_Kmer_Index(S,S->cidx-4);
              printf("\n");
	    }
#endif
#ifdef DEBUG_CMERGE
          *kptr = x;
          printf("c = %4d k = ",c);
          write_Ascii(kptr,KMER);
          printf(" %d  e = %s",*((uint16 *) (kptr+KMER_BYTES)),fmer[S->cpre]);
          write_Ascii(S->csuf,KMER-4);
          printf("  %d %lld\n",*((uint16 *) (S->csuf+hbyte)),S->cidx);
#endif
          if (c != 0)
            ct = 0;
          else // c == 0
            { ct = *((uint16 *) (S->csuf+hbyte));
              Next_Kmer_Entry(S);
            }

          while (kptr < lptr)
            { d = kptr[KM1];
              fill = fours[d];
              *((uint16 *) fill) = ct;
              fill += 2;
              for (k = TMER_WORD; k < KM1; k++)
                *fill++ = kptr[k];
              fours[d] = fill;
              kptr += KMER_WORD;
            }
        } 
    } 

#ifndef DEVELOPER
  if (beg > 0)
    Free_Kmer_Stream(S);
#endif

  return (NULL);
}


/*******************************************************************************************
 *
 * static void *profile_list_thread(Plist_Arg *arg)
 *     Each thread takes the now "un"sorted weighted k-mers entries that are now in order
 *     of the k-mers along each super-mer.  Produce a profile for each unique super-mer
 *     in place and build an array of entries consisting of a pointer to the profile fragment
 *     and the ordinal super-mer id of each of the equall super-mers with that profile.
 *     This array will be "un"sorted on super-mer id.
 *
 ********************************************************************************************/

typedef struct
  { uint8    *sort;
    int64    *parts;
    int       beg;
    int       end;
    int64     off;
    uint8    *prol;
    int64     cnts;
    uint8    *fill;  //  profile list uses this (basaed on the counts)
  } Plist_Arg;

static void *profile_list_thread(void *arg)
{ Plist_Arg   *data   = (Plist_Arg *) arg;
  int          beg    = data->beg;
  int          end    = data->end;
  int64       *part   = data->parts;

  uint8    *cnts, *fill, *prof;
  uint8    *sptr, *send;
  int       i, x;
  uint8    *asp;
  uint8     b[2*MAX_SUPER+2];
  uint64    pidx;

  int       STOT = SMER_BYTES + SLEN_BYTES;

  int sln = 0;
  int len = 0;
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
  uint8  *sb = ((uint8 *) &sln) - 1;
  uint8  *lb = ((uint8 *) &len) - 1;
#else
  uint8  *sb = ((uint8 *) &sln) + (sizeof(int)-SLEN_BYTES);
  uint8  *lb = ((uint8 *) &len) + (sizeof(int)-PLEN_BYTES);
#endif

#ifdef SHOW_RUN
  int64  ridx = 0;
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
  uint8  *rbx = ((uint8 *) &ridx) - 1;
#else
  uint8  *rbx = ((uint8 *) &ridx) + (sizeof(int64)-RUN_BYTES);
#endif
#endif

  prof = data->prol;
  pidx = CMER_WORD*data->cnts;
  cnts = data->prol + pidx;
  sptr = data->sort + data->off;
  fill = data->fill;
  for (x = beg; x < end; x++)
    for (send = sptr + part[x]; sptr < send; )

      { asp = sptr + SMER_BYTES;
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
        for (i = SLEN_BYTES; i > 0; i--)
          sb[i] = *asp++;
#else
        for (i = 0; i < SLEN_BYTES; i++)
          sb[i] = *asp++;
#endif

        { int j, p, c, d;
          int run;
          uint8 *db = (uint8 *) &d;

          p = *((uint16 *) cnts);
          cnts += CMER_WORD;

          *((uint16 *) b) = p;
          len = 2;
#ifdef SHOW_RUN
          printf("  {%d}",p);
#endif

          run = 0;
          for (j = 1; j <= sln; j++)
            { c = *((uint16 *) cnts);
              cnts += CMER_WORD;

              if (c == p)
                if (run > 0)
                  { if (run >= 63)
                      { b[len++] = run;
#ifdef SHOW_RUN
                        printf(" [%d]",run);
#endif
                        run = 1;
                      }
                    else
                      run += 1;
                  }
                else
                  run = 1;
              else
                { if (run > 0)
                    { b[len++] = run;
#ifdef SHOW_RUN
                      printf(" [%d]",run);
#endif
                      run = 0;
                    }
                  d = c-p;
#ifdef SHOW_RUN
                  printf(" %d",d);
#endif
                  if (abs(d) < 32)
                    b[len++] = 0x40 | (d & 0x3f);
                  else
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
                    { b[len++] = db[1] | 0x80;
                      b[len++] = db[0];
#else
                    { b[len++] = db[0] | 0x80;
                      b[len++] = db[1];
#endif
#ifdef SHOW_RUN
                      printf("+");
#endif
                    }
                }

              p = c;
            }

          if (run > 0)
            { b[len++] = run;
#ifdef SHOW_RUN
              printf(" [%d]",run);
#endif
            }
          if (len > 2)
            { *((uint16 *) (b+len)) = p;
#ifdef SHOW_RUN
              printf(" {%d}",p);
#endif
            }

#ifdef SHOW_RUN
          printf(" += %d / %d\n",len,sln+1);
#endif
        }

#ifdef SHOW_RUN
        printf("Run ids:");
#endif

        do
          { if (sptr[STOT] & 0x80)
              { *((uint64 *) fill) = (pidx<<1) | 0x1;
                sptr[STOT] &= 0x7f;
#ifdef SHOW_RUN
                printf(" +");
#endif
              }
            else
              { *((uint64 *) fill) = (pidx<<1);
#ifdef SHOW_RUN
                printf("  ");
#endif
              }
            fill += sizeof(uint64);
            for (i = STOT; i < SMER_WORD; i++)
              *fill++ = sptr[i];

#ifdef SHOW_RUN
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
            for (i = RUN_BYTES; i > 0; i--)
              rbx[i] = sptr[SMER_WORD-i];
#else
            for (i = 0; i < PLEN_BYTES; i++)
              rbx[i] = sptr[STOT+i];
#endif
            printf("%lld",ridx);
#endif

            sptr += SMER_WORD;
          }
        while (*sptr == 0);

#ifdef SHOW_RUN
        printf(" <%ld>\n",(sptr-asp)/SMER_WORD+1);
#endif

#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
        for (i = PLEN_BYTES; i > 0; i--)
          prof[pidx++] = lb[i];
#else
        for (i = 0; i < PLEN_BYTES; i++)
          prof[pidx++] = lb[i];
#endif
        if (len > 2)
          len += 2;
        for (i = 0; i < len; i++)
          prof[pidx++] = b[i];
      }

  return (NULL);
}


/*******************************************************************************************
 *
 * static void *profile_write_thread(Pwrite_Arg *arg)
 *     Each thread takes the now "un"sorted list of pointers to profile fragments and
 *     outputs them to a "P" file in the SORT_PATH directory.
 *
 ********************************************************************************************/

typedef struct
  { uint8    *sort;
    char     *root;
    int       wch;
    int64     beg;
    int64     end;
    uint8    *prol;
    int64     nidxs;
  } Pwrite_Arg;

static void *profile_write_thread(void *arg)
{ Pwrite_Arg  *data   = (Pwrite_Arg *) arg;
  uint8       *prol   = data->prol;
  int64        beg    = data->beg;
  int64        rng    = data->end - beg;

  int    pfile;
  char  *fname;
  uint8 *psort, *next;
  uint8 *prof;
  uint8  bufr[0x10000];
  uint8 *fill, *bend;
  uint64 pidx;
  int    last;
  int    t, k;

  int len = 0;
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
  uint8  *lb = ((uint8 *) &len) - 1;
#else
  uint8  *lb = ((uint8 *) &len) + (sizeof(int)-PLEN_BYTES);
#endif

#ifdef DEBUG_PWRITE
  int64   r = 0;
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
  uint8  *rb = ((uint8 *) &r) - 1;
#else
  uint8  *rb = ((uint8 *) &r) + (sizeof(int64)-RUN_BYTES);
#endif
#endif

  fname = Malloc(strlen(data->root) + 100,"File Name");
  bend  = bufr + (0x10000 - (RUN_BYTES+PLEN_BYTES+2*MAX_SUPER+2));

  psort = data->sort + beg*PROF_BYTES;
  for (t = 1; t <= NPANELS; t++)
    { next  = data->sort + (beg + (rng*t)/NPANELS) * PROF_BYTES;

#ifdef DEBUG_PWRITE
      printf("Panel %d\n",t);
#endif
      sprintf(fname,"%s%d.%d",data->root,data->wch,t-1);
      pfile = open(fname,O_WRONLY|O_CREAT|O_TRUNC,S_IRWXU|S_IRWXG|S_IRWXO);
      if (pfile < 0)
        { fprintf(stderr,"\n%s: Could not open %s for writing\n",Prog_Name,fname);
          Clean_Exit(1);
        }

#ifdef DEVELOPER
      if (t == 1)
        { if (data->wch == 0)
            { write(pfile,&RUN_BYTES,sizeof(int));
              write(pfile,&PLEN_BYTES,sizeof(int));
              write(pfile,&MAX_SUPER,sizeof(int));
            }
          if (write(pfile,&data->nidxs,sizeof(int64)) < 0)
            { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",Prog_Name,fname);
              Clean_Exit(1);
            }
        }
#endif

      fill = bufr;
      while (psort < next)
        { if (fill >= bend)
            { if (write(pfile,bufr,fill-bufr) < 0)
                { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",Prog_Name,fname);
                  Clean_Exit(1);
                }
              fill = bufr;
            }

          pidx = *((uint64 *) psort);
          psort += sizeof(uint64);

          last = (pidx & 0x1);
          pidx >>= 1;

#ifdef DEBUG_PWRITE
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
          for (k = RUN_BYTES; k > 0; k--)
            rb[k] = psort[RUN_BYTES-k];
#else
          for (k = 0; k < RUN_BYTES; k++)
            rb[k] = psort[k];
#endif
#endif

          if (last)
            *fill++ = *psort++ | 0x80;
          else
            *fill++ = *psort++;
          for (k = 1; k < RUN_BYTES; k++)
            *fill++ = *psort++;

          prof = prol + pidx;

#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
          for (k = PLEN_BYTES; k > 0; k--)
            lb[k] = *fill++ = *prof++;
#else
          for (k = 0; k < PLEN_BYTES; k++)
            lb[k] = *fill++ = *prof++;
#endif

#ifdef DEBUG_PWRITE
          { uint16 x, d;

            printf("%8ld: %8lld%c (%9lld), %3d:",
                   (psort-data->sort)/PROF_BYTES,r,last?'+':' ',pidx,len);

            x = *((uint16 *) prof);
            k = 2;
            printf(" {%d}",x);
            while (k < len)
              { x = prof[k++];
                if ((x & 0x80) != 0)
                  { if ((x & 0x40) != 0)
                      d = (x << 8);
                    else
                      d = (x << 8) & 0x7fff;
                    x = prof[k++];
                    d |= x;
                    printf(" %hd+",d);
                  }
                else if ((x & 0x40) != 0)
                  { if ((x & 0x20) != 0)
                      printf(" -%d",32-(x&0x1fu));
                    else
                      printf(" +%d",x&0x1fu);
                  }
                else
                  printf(" [%hu]",x);
              }
            if (len > 2)
              printf(" {%d}",*((uint16 *) (prof+len)));
            printf("\n");
          }
#endif

          if (len > 2)
            len += 2;
          for (k = 0; k < len; k++)
            *fill++ = *prof++;
        }

      if (write(pfile,bufr,fill-bufr) < 0)
        { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",Prog_Name,fname);
          Clean_Exit(1);
        }
      close(pfile);
    }

  free(fname);
  return (NULL);
}


/*********************************************************************************************\
 *
 *  Sorting(char *path, char *root)
 *
 *     For each set of the NTHREADS files for a each partition bucket produced
 *     by the split.c module, do the following:
 *       * Sort the super-mers.
 *       * Sort the k-mers from each super-mer weighted by the # of times that super-mer occurs.
 *       * During the 2nd sort accumulate the histogram of k-mer frequencies.
 *            Output this to <source_root>.K<kmer>.
 *       * if requesteda (-t) produce a table of all the k-mers with counts >= -t in NTHREADSs
 *            pieces in files SORT_PATH/<root>.<bucket>.L<thread>
 *       * if requested (-p) invert the first two sorts to produce a profile for every super-mer
 *            in the order in the source in files SORT_PATH/<root>.<bucket>.P<thread>.[0-3]
 *       * if requested (-h) print the histogram of k-mer frequencies.
 *
 *********************************************************************************************/

void Sorting(char *path, char *root)
{ char  *fname;
  int64  counts[0x8000];
  int64  max_inst;
  int   *reload;

  if (VERBOSE)
    { fprintf(stderr,"\nPhase 2: Sorting & Counting K-mers in %d blocks\n\n",NPARTS);
      fflush(stderr);
    }

  fname = Malloc(2*(strlen(SORT_PATH) + strlen(path) + strlen(root)) + 100,"File name buffer");
  if (fname == NULL)
    Clean_Exit(1);

  //  First bundle: initialize all sizes & lookup tables

  { int  i, n;
    int  l0, l1, l2, l3;
    int *s;

    reload = Malloc(sizeof(int)*IO_UBITS*MAX_SUPER,"Allocating reload table");

    for (i = 0; i < 0x8000; i++)
      counts[i] = 0;
    max_inst = 0;

    s = reload;
    for (i = 0; i < IO_UBITS; i++)
      { Fixed_Reload[IO_UBITS-i] = (i + SLEN_BITS - 1) / IO_UBITS;
        Runer_Reload[IO_UBITS-i] = (i + RUN_BITS - 1) / IO_UBITS;
        Super_Reload[IO_UBITS-i] = s;
        for (n = 0; n < MAX_SUPER; n++)
          *s++ = (i + 2*(n+KMER) - 1) / IO_UBITS;
      }

    i = 0;
    for (l0 = 3; l0 >= 0; l0 -= 1)
     for (l1 = 12; l1 >= 0; l1 -= 4)
      for (l2 = 48; l2 >= 0; l2 -= 16)
       for (l3 = 192; l3 >= 0; l3 -= 64)
         Comp[i++] = (l3 | l2 | l1 | l0);

#if defined(DEBUG_CANONICAL) || defined(DEBUG_KLIST) || defined(DEBUG_SLIST) \
      || defined(DEBUG_TABOUT) || defined(DEBUG_CMERGE) || defined(EQUAL_MERGE)
    { char *t;

      i = 0;
      t = _fmer;
      for (l3 = 0; l3 < 4; l3++)
       for (l2 = 0; l2 < 4; l2++)
        for (l1 = 0; l1 < 4; l1++)
         for (l0 = 0; l0 < 4; l0++)
           { fmer[i] = t;
             *t++ = DNA[l3];
             *t++ = DNA[l2];
             *t++ = DNA[l1];
             *t++ = DNA[l0];
             *t++ = 0;
             i += 1;
           }
    }
#endif
  }

  // For each partition, make k-mer list, sort on k-mer and count, resort on position,
  //    and finally output super-counts

  { Slist_Arg  *parms = Malloc(sizeof(Slist_Arg)*ITHREADS,"Allocating sort controls");
    Klist_Arg  *parmk = Malloc(sizeof(Klist_Arg)*NTHREADS,"Allocating sort controls");
    Clist_Arg  *parmc = Malloc(sizeof(Clist_Arg)*NTHREADS,"Allocating sort controls");
    Twrite_Arg *parmt = Malloc(sizeof(Twrite_Arg)*NTHREADS,"Allocating sort controls");
    Plist_Arg  *parmp = Malloc(sizeof(Plist_Arg)*NTHREADS,"Allocating sort controls");
    Pwrite_Arg *parmw = Malloc(sizeof(Pwrite_Arg)*ITHREADS,"Allocating sort controls");

    int   *Table_Split = Malloc(sizeof(int)*NTHREADS,"Allocating sort controls");
    int64 *Sparts      = Malloc(sizeof(int64)*256,"Allocating sort controls");
    int64 *Kparts      = Malloc(sizeof(int64)*256,"Allocating sort controls");
    Range *Panels      = Malloc(sizeof(Range)*NTHREADS,"Allocating sort controls");
    int64 *Wkmers      = Malloc(sizeof(int64)*NPARTS,"Allocating sort controls");
    int64 *Ukmers      = Malloc(sizeof(int64)*NPARTS,"Allocating sort controls");

#if !defined(DEBUG) || !defined(SHOW_RUN)
    THREAD *threads = Malloc(sizeof(THREAD)*NTHREADS,"Allocating sort controls");
#endif

    int         ODD_PASS = 0;

    uint8      *s_sort;
    uint8      *k_sort;
    uint8      *i_sort;
    uint8      *p_sort;
    uint8      *a_sort;

    int64  kmers;
    int64  nmers;
    int64  skmers;
    int64  tmers;
    int    t, p;

    s_sort = NULL;
    i_sort = NULL;

    if (parms == NULL || parmk == NULL || parmc == NULL ||
        parmt == NULL || parmp == NULL || parmw == NULL)
      Clean_Exit(1);

    if (Table_Split == NULL || Sparts == NULL || Kparts == NULL ||
        Panels == NULL || Wkmers == NULL || Ukmers == NULL)
      Clean_Exit(1);

#if !defined(DEBUG) || !defined(SHOW_RUN)
    if (threads == NULL)
      Clean_Exit(1);
#endif

#ifndef DEVELOPER
    if (DO_PROFILE)
      { KMER_WORD += KMAX_BYTES;
        CMER_WORD  = KMAX_BYTES+1;
        ODD_PASS   = (KMAX_BYTES % 2 == 1);
        RUN_BYTES  = (RUN_BITS+7) >> 3;
        SMER_WORD += RUN_BYTES;
        PROF_BYTES = RUN_BYTES + sizeof(uint64);
      }
    s_sort = Malloc((NMAX+1)*SMER_WORD+1,"Allocating super-mer sort array");
    if (s_sort == NULL)
      Clean_Exit(1);
    *s_sort++ = 0;
#endif

    tmers = 0;
    for (p = 0; p < NPARTS; p++)
      {
        //  Open and reada headers of super-mer files for part p

        kmers = 0;
        nmers = 0;
        for (t = 0; t < ITHREADS; t++)
          { int64 k, n;
            int   f;

            sprintf(fname,"%s/%s.%d.T%d",SORT_PATH,root,p,t);
            f = open(fname,O_RDONLY);
            if (f < 0)
              { fprintf(stderr,"\n%s: File %s should exist but doesn't?\n",Prog_Name,fname); 
                Clean_Exit(1);
              }

            parms[t].tfile = f;
#ifdef DEVELOPER
            read(f,&KMAX,sizeof(int64));
            read(f,&NMAX,sizeof(int64));
            read(f,&KMAX_BYTES,sizeof(int));
            read(f,&RUN_BITS,sizeof(int));
            read(f,&n,sizeof(int64));
            parms[t].nidxs = n;
#else
            parms[t].nidxs = NUM_RID[t];
#endif
            read(f,&k,sizeof(int64));
            read(f,&n,sizeof(int64));
            parms[t].nmers = n;
            kmers += k;
            nmers += n;

            read(f,Panels[t].khist,sizeof(int64)*256);
          }

#ifdef DEVELOPER
        if (p == 0)
          { int i;

            if (DO_PROFILE)
              { KMER_WORD += KMAX_BYTES;
                CMER_WORD  = KMAX_BYTES+1;
                ODD_PASS   = (KMAX_BYTES % 2 == 1);
                RUN_BYTES  = (RUN_BITS+7) >> 3;
                SMER_WORD += RUN_BYTES;
                PROF_BYTES = RUN_BYTES + sizeof(uint64);
                for (i = 0; i < IO_UBITS; i++)
                  Runer_Reload[IO_UBITS-i] = (i + RUN_BITS - 1) / IO_UBITS;
              }
            s_sort = Malloc((NMAX+1)*SMER_WORD+1,"Allocating super-mer sort array");
            if (s_sort == NULL)
              Clean_Exit(1);
            *s_sort++ = 0;
          }
#endif

        //  Build super-mer list

        if (VERBOSE)
          { fprintf(stderr,"\r  Processing block %d: Sorting super-mers     ",p+1); 
            fflush(stderr);
          }

        { int64 o, x;
          int   j;

          o = 0;
          for (j = 0; j < 256; j++)
            for (t = 0; t < ITHREADS; t++)
              { parms[t].fours[j] = s_sort + o*SMER_WORD;
                o += Panels[t].khist[j];
              }

          o = 0;
          for (t = 0; t < ITHREADS; t++)
            { x = parms[t].nidxs;
              parms[t].nbase = o;
              o += x;
            }
        }

#ifdef DEBUG_SLIST
        for (t = 0; t < ITHREADS; t++)
          supermer_list_thread(parms+t);
#else
        for (t = 1; t < ITHREADS; t++)
          pthread_create(threads+t,NULL,supermer_list_thread,parms+t);
	supermer_list_thread(parms);
        for (t = 1; t < ITHREADS; t++)
          pthread_join(threads[t],NULL);
#endif

        for (t = 0; t < ITHREADS; t++)
          close(parms[t].tfile);

#ifndef DEVELOPER
        for (t = 0; t < ITHREADS; t++)
          { sprintf(fname,"%s/%s.%d.T%d",SORT_PATH,root,p,t);
            unlink(fname);
          }
#endif

        //  Sort super-mer list

        { uint8 *o, *x;
          int    j;

          o = s_sort;
          for (j = 0; j < 256; j++)
            { x = parms[ITHREADS-1].fours[j];
              Sparts[j] = x-o;
              o = x;
            }
        }

        Supermer_Sort(s_sort,nmers,SMER_WORD,SMER_BYTES+SLEN_BYTES,Sparts,NTHREADS,Panels);

        //  Allocate and fill in weighted k-mer list from sorted supermer list

        { int64 o, x;
          int   j;

          if (DO_PROFILE)
            { o = 0;
              for (t = 0; t < NTHREADS; t++)
                { parmk[t].kidx = o;
                  for (j = 0; j < 256; j++)
                    o += Panels[t].khist[j];
                }
            }

          o = 0;
          for (j = 0; j < 256; j++)
            for (t = 0; t < NTHREADS; t++)
              { x = Panels[t].khist[j];
                Panels[t].khist[j] = o;
                o += x;
              }
          skmers = o;
        }

        if (VERBOSE)
          { fprintf(stderr,"\r  Processing block %d: Sorting weighted k-mers",p+1); 
            Wkmers[p] = skmers;
            Ukmers[p] = kmers;
            fflush(stderr);
          }

        if (DO_PROFILE)
          if (ODD_PASS)
            { i_sort = Malloc(skmers*(CMER_WORD+KMER_WORD)+2,"Weighted k-mer & inverse list");
              k_sort = i_sort + skmers*CMER_WORD + 1;
            }
          else
            { k_sort = Malloc(skmers*(KMER_WORD+CMER_WORD)+2,"Weighted k-mer & inverse list");
              i_sort = k_sort + skmers*KMER_WORD + 1;
            }
        else
          k_sort = Malloc(skmers*KMER_WORD+1,"Weighted k-mer list");

        for (t = 0; t < NTHREADS; t++)
          { int j;

            parmk[t].sort   = s_sort;
            parmk[t].parts  = Sparts;
            parmk[t].beg    = Panels[t].beg;
            parmk[t].end    = Panels[t].end;
            parmk[t].off    = Panels[t].off;
            for (j = 0; j < 256; j++)
              parmk[t].fours[j] = k_sort + Panels[t].khist[j]*KMER_WORD;
          }

#if defined(DEBUG_KLIST) || defined(DEBUG_CANONICAL)
        for (t = 0; t < NTHREADS; t++)
          kmer_list_thread(parmk+t);
#else
        for (t = 1; t < NTHREADS; t++)
          pthread_create(threads+t,NULL,kmer_list_thread,parmk+t);
	kmer_list_thread(parmk);
        for (t = 1; t < NTHREADS; t++)
          pthread_join(threads[t],NULL);
#endif

        //  Sort weighted k-mer list

        { uint8 *o, *x;
          int    j;

          o = k_sort;
          for (j = 0; j < 256; j++)
            { x = parmk[NTHREADS-1].fours[j];
              Kparts[j] = x-o;
              o = x;
            }
        }

        Weighted_Kmer_Sort(k_sort,skmers,KMER_WORD,KMER_BYTES,Kparts,NTHREADS,Panels);

        //  Accumulate frequency histogram across threads

        if (PRO_TABLE == NULL)
          { int64 *ncnt;
            int    x;

            for (t = 0; t < NTHREADS; t++)
              { ncnt = Panels[t].count;
                for (x = 1; x < 0x8000; x++)
                  counts[x] += ncnt[x];
                max_inst += Panels[t].max_inst + parmk[t].overflow;
              }
          }

        if (DO_TABLE > 0)
          {
            //  Threaded write of sorted kmer+count table
            //    1st time determine partition bytes values.  Therafter, spit on said

            if (p == 0)
              { for (t = 0; t < NTHREADS; t++)
                  { Table_Split[t] = Panels[t].beg;
                    parmt[t].off = Panels[t].off;
                  }
              }
            else
              { int64 off;
                int   beg;

                off = 0;
                beg = 0;
                for (t = 0; t < NTHREADS; t++)
                  { while (beg < Table_Split[t])
                      { off += Kparts[beg];
                        beg += 1;
                      }
                    parmt[t].off = off;
                  }
              }

            for (t = 0; t < NTHREADS; t++)
              { parmt[t].sort  = k_sort;
                parmt[t].parts = Kparts;
                parmt[t].beg   = Table_Split[t];
                if (t < NTHREADS-1)
                  parmt[t].end  = Table_Split[t+1];
                else
                  parmt[t].end = 256;
                sprintf(fname,"%s/%s.%d.L%d",SORT_PATH,root,p,t);
                parmt[t].kfile = open(fname,O_WRONLY|O_CREAT|O_TRUNC,S_IRWXU|S_IRWXG|S_IRWXO);
                parmt[t].kname = Strdup(fname,"Allocating stream name");
                if (parmt[t].kname == NULL)
                  Clean_Exit(1);
              }
#ifdef DEVELOPER
            if (p == NPARTS-1)
              { int zero = 0;
                if (write(parmt[0].kfile,&zero,sizeof(int)) < 0)
                  { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",
                                   Prog_Name,parmt[0].kname);
                    Clean_Exit(1);
                  }
              }
#endif

#ifdef DEBUG_TABOUT
            for (t = 0; t < NTHREADS; t++)
              table_write_thread(parmt+t);
#else
            for (t = 1; t < NTHREADS; t++)
              pthread_create(threads+t,NULL,table_write_thread,parmt+t);
	    table_write_thread(parmt);
            for (t = 1; t < NTHREADS; t++)
              pthread_join(threads[t],NULL);
#endif

            for (t = 0; t < NTHREADS; t++)
              tmers += parmt[t].tmers;

            if (p == NPARTS-1)
              { if (tmers > 0x4000000ll && KMER >= 12)
                  IDX_BYTES = 3;
                else if (tmers >= 0x40000ll && KMER >= 8)
                  IDX_BYTES = 2;
                else // tmers < 0x40000ll, KMER always >= 7
                  IDX_BYTES = 1;
#ifdef DEVELOPER
                lseek(parmt[0].kfile,0,SEEK_SET);
                if (write(parmt[0].kfile,&IDX_BYTES,sizeof(int)) < 0)
                  { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",
                                   Prog_Name,parmt[0].kname);
                    Clean_Exit(1);
                  }
#endif
              }
 
            for (t = 0; t < NTHREADS; t++)
              { free(parmt[t].kname);
                close(parmt[t].kfile);
              }
          }

        if (! DO_PROFILE)
	  { free(k_sort);
            continue;
          }

        if (VERBOSE)
          { fprintf(stderr,"\r  Processing block %d: Inverse profile sorting",p+1); 
            fflush(stderr);
          }

        //  Fill in count/index list from sorted k-mer list, pre-sorted on
        //    LSD byte of index.

        { int64 o, x;
          int   j;

          o = 0;
          for (j = 0; j < 256; j++)
            for (t = 0; t < NTHREADS; t++)
              { x = Panels[t].khist[j];
                Panels[t].khist[j] = o;
                o += x;
              }
        }

        for (t = 0; t < NTHREADS; t++)
          { int j;
   
            parmc[t].sort   = k_sort;
            parmc[t].parts  = Kparts;
            parmc[t].beg    = Panels[t].beg;
            parmc[t].end    = Panels[t].end;
            parmc[t].off    = Panels[t].off;
            for (j = 0; j < 256; j++)
              parmc[t].fours[j] = i_sort + Panels[t].khist[j]*CMER_WORD;
          }

        if (PRO_TABLE != NULL)
          {
            // Relative profile: also set up table file for merges

            sprintf(fname,"%s/%s.U%d",SORT_PATH,root,p);
            parmc[0].stm = Open_Kmer_Stream(fname);
            if (parmc[0].stm == NULL)
              { fprintf(stderr,"\n%s: Table %s should exist but doesn't?\n",Prog_Name,fname); 
                Clean_Exit(1);
              }
            for (t = 1; t < NTHREADS; t++)
              parmc[t].stm = parmc[0].stm;

#if defined(DEBUG_CMERGE) || defined(EQUAL_MERGE)
            for (t = 0; t < NTHREADS; t++)
              cmer_merge_thread(parmc+t);
#else
            for (t = 1; t < NTHREADS; t++)
              pthread_create(threads+t,NULL,cmer_merge_thread,parmc+t);
            cmer_merge_thread(parmc);
            for (t = 1; t < NTHREADS; t++)
              pthread_join(threads[t],NULL);
#endif

            Free_Kmer_Stream(parmc[0].stm);

#ifndef DEVELOPER
            sprintf(fname,"rm -f %s/%s.U%d.ktab %s/.%s.U%d.ktab.*",
                          SORT_PATH,root,p,SORT_PATH,root,p);
            system(fname);
#endif
          }

        else // PRO_TABLE == 0
          {
#ifdef DEBUG_CLIST
            for (t = 0; t < NTHREADS; t++)
              cmer_list_thread(parmc+t);
#else
            for (t = 1; t < NTHREADS; t++)
              pthread_create(threads+t,NULL,cmer_list_thread,parmc+t);
            cmer_list_thread(parmc);
            for (t = 1; t < NTHREADS; t++)
              pthread_join(threads[t],NULL);
#endif
          }

        //  LSD sort count/index list on index and then tidy up memory

        { int i, x;
          int bytes[KMAX_BYTES+1];

          x = 0;
          for (i = KMAX_BYTES; i >= 2; i--)
            bytes[x++] = i;
          bytes[x] = -1;

          i_sort = LSD_Sort(skmers,i_sort,k_sort,CMER_WORD,bytes);

          if (ODD_PASS)
            i_sort = Realloc(i_sort,skmers*CMER_WORD,"Pruning count list");
          else
            i_sort = Realloc(k_sort,skmers*CMER_WORD,"Pruning count list");
        }

        //  Use i_sort & k_sort again to build list of compressed profile fragments
        //    in place in i_sort, and build reference list

        p_sort = Malloc(nmers*PROF_BYTES*2,"Allocating profile link array");

        for (t = 0; t < NTHREADS; t++)
          { parmp[t].sort   = s_sort;
            parmp[t].parts  = Sparts;
            parmp[t].beg    = parmk[t].beg;
            parmp[t].end    = parmk[t].end;
            parmp[t].off    = parmk[t].off;
            parmp[t].prol   = i_sort;
            parmp[t].cnts   = parmk[t].kidx;
            parmp[t].fill   = p_sort + (parmk[t].off / SMER_WORD) * PROF_BYTES;
          }

#if defined(DEBUG_PLIST) || defined(SHOW_RUN)
        for (t = 0; t < NTHREADS; t++)
          profile_list_thread(parmp+t);
#else
        for (t = 1; t < NTHREADS; t++)
          pthread_create(threads+t,NULL,profile_list_thread,parmp+t);
        profile_list_thread(parmp);
        for (t = 1; t < NTHREADS; t++)
          pthread_join(threads[t],NULL);
#endif

        //  LSD sort profile links on super-mer idx

        { int   i, x;
          int   bytes[PROF_BYTES+1];

          x = 0;
          for (i = PROF_BYTES-1; i >= (int) sizeof(uint64); i--)
            bytes[x++] = i;
          bytes[x] = -1;

          a_sort = LSD_Sort(nmers,p_sort,p_sort+nmers*PROF_BYTES,PROF_BYTES,bytes);
        }

        //  Output profile fragments in order of a_sort links

        { int64 o;

          sprintf(fname,"%s/%s.%d.P",SORT_PATH,root,p);
          o = 0;
          for (t = 0; t < ITHREADS; t++)
            { parmw[t].sort  = a_sort;
              parmw[t].beg   = o;
              o += parms[t].nmers;
              parmw[t].end   = o;
              parmw[t].nidxs = parms[t].nidxs;
#ifdef DEBUG_PWRITE
              printf("Partition %2d: %10lld [%lld]\n",t,o,nmers);
#endif
              parmw[t].prol  = i_sort;
              parmw[t].root  = fname;
              parmw[t].wch   = t;
            }
        }

#ifdef DEBUG_PWRITE
        for (t = 0; t < ITHREADS; t++)
          profile_write_thread(parmw+t);
#else
        for (t = 1; t < ITHREADS; t++)
          pthread_create(threads+t,NULL,profile_write_thread,parmw+t);
        profile_write_thread(parmw);
        for (t = 1; t < ITHREADS; t++)
          pthread_join(threads[t],NULL);
#endif

	free(p_sort);
	free(i_sort);
      }

    free(s_sort-1);

    if (VERBOSE)
      { int64  wtot, utot;
        double psav;
        int    wwide, awide;

        wtot = utot = 0;
        psav = 0.;
        for (p = 0; p < NPARTS; p++)
          { wtot += Wkmers[p];
            if (Wkmers[p] > 0)
              { if ((1.*Ukmers[p])/Wkmers[p] > psav)
                  psav = (1.*Ukmers[p])/Wkmers[p];
              }
            utot += Ukmers[p];
          }
        if (wtot > 0)
          { if ((1.*utot)/wtot > psav)
              psav = (1.*utot)/wtot;
          }

        wwide = Number_Digits(wtot);
        awide = Number_Digits((int64) psav); 
        wwide += (wwide-1)/3;
        awide += 2;
        if (12 > wwide)
          wwide = 12;
        if (7 > awide)
          awide = 7;

        fprintf(stderr,"\r                                               \r"); 
        fprintf(stderr,"      Part:%*swgt'd k-mers%*ssavings\n",wwide-10,"",awide-5,"");
        fflush(stderr);
        for (p = 0; p < NPARTS; p++)
          { fprintf(stderr,"     %5d:  ",p);
            Print_Number(Wkmers[p],wwide,stderr);
            if (Wkmers[p] > 0)
              fprintf(stderr,"  %*.1f\n",awide,(1.*Ukmers[p])/Wkmers[p]);
            else
              fprintf(stderr,"  %*sNA\n",awide-2,"");
          }
        fprintf(stderr,"       All:  ");
        Print_Number(wtot,wwide,stderr);
        if (wtot > 0)
          fprintf(stderr,"  %*.1f\n",awide,(1.*utot)/wtot);
        else
          fprintf(stderr,"  %*sNA\n",awide-2,"");
        fflush(stderr);
      }

#if !defined(DEBUG) || !defined(SHOW_RUN)
    free(threads);
#endif

    free(Ukmers);
    free(Wkmers);
    free(Panels);
    free(Kparts);
    free(Sparts);
    free(Table_Split);

    free(parmw);
    free(parmp);
    free(parmt);
    free(parmc);
    free(parmk);
    free(parms);
  }

  //  Output histogram

  if (PRO_TABLE == NULL)
    { int   i, f;

      sprintf(fname,"%s/%s.hist",path,root);
      f = open(fname,O_WRONLY|O_CREAT|O_TRUNC,S_IRWXU|S_IRWXG|S_IRWXO);
      write(f,&KMER,sizeof(int));
      i = 1;
      write(f,&i,sizeof(int));
      i = 0x7fff;
      write(f,&i,sizeof(int));
      write(f,counts+1,sizeof(int64));
      write(f,&max_inst,sizeof(int64));
      if (write(f,counts+1,0x7fff*sizeof(int64)) < 0)
        { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",Prog_Name,fname);
          Clean_Exit(1);
        }
      close(f);
    }

  free(reload);
  free(fname);
}
