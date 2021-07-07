/*********************************************************************************************\
 *
 *  Given a conanonical table, produce a "symmetric" table that includes non-canonical k-mers
 *
 *  Author:  Gene Myers
 *  Date  :  April, 2021
 *
 *********************************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>

#undef  DEBUG

#include "libfastk.h"
#include "FastK.h"

int   VERBOSE;
int   NTHREADS;
char *SORT_PATH;

static char *Usage = " [-v] [-T<int(4)>] [-P<dir(/tmp)] <source_root>[.ktab] <dest_root>[.ktab]";


/****************************************************************************************
 *
 *  Print & compare utilities
 *
 *****************************************************************************************/

#define  COUNT_PTR(p) ((uint16 *) (p+kbyte))

static char dna[4] = { 'a', 'c', 'g', 't' };

static char *fmer[256], _fmer[1280];

static void setup_fmer_table()
{ char *t;
  int   i, l3, l2, l1, l0;

  i = 0;
  t = _fmer;
  for (l3 = 0; l3 < 4; l3++)
   for (l2 = 0; l2 < 4; l2++)
    for (l1 = 0; l1 < 4; l1++)
     for (l0 = 0; l0 < 4; l0++)
       { fmer[i] = t;
         *t++ = dna[l3];
         *t++ = dna[l2];
         *t++ = dna[l1];
         *t++ = dna[l0];
         *t++ = 0;
         i += 1;
       }
}

static uint8 comp[256];

static void setup_comp_table()
{ int  i, l3, l2, l1, l0;

  i = 0;
  for (l3 = 3; l3 >= 0; l3 -= 1)
   for (l2 = 12; l2 >= 0; l2 -= 4)
    for (l1 = 48; l1 >= 0; l1 -= 16)
     for (l0 = 192; l0 >= 0; l0 -= 64)
       comp[i++] = l0 | l1 | l2 | l3;
}

static void print_seq(uint8 *seq, int len)
{ int i, b, k;

  b = (len >> 2);
  for (i = 0; i < b; i++)
    printf("%s",fmer[seq[i]]);
  k = 6;
  for (i = b << 2; i < len; i++)
    { printf("%c",dna[(seq[b] >> k) & 0x3]);
      k -= 2;
    }
}

/****************************************************************************************
 *
 *  Distribute duplicated entries for sorts
 *
 *****************************************************************************************/

#define BUFFER_SIZE 0x1000000;

typedef struct
  { int64  nels;
    int    fid;
    int    bptr;
    char  *buff;
  } O_Block;


  //  Writes over 2GB don't work on some systems, patch to overcome said

static inline int64 big_write(int f, uint8 *buffer, int64 bytes)
{ int64 v, x;

  v = 0;
  while (bytes > 0x70000000)
    { x = write(f,buffer,0x70000000);
      if (x < 0)
        return (-1);
      v += x;
      bytes  -= x;
      buffer += x;
    }
  x = write(f,buffer,bytes);
  if (x < 0)
    return (-1);
  return (v+x);
}

  //  Reads over 2GB don't work on some systems, patch to overcome said

static inline int64 big_read(int f, uint8 *buffer, int64 bytes)
{ int64 v, x;

  v = 0;
  while (bytes > 0x70000000)
    { x = read(f,buffer,0x70000000);
      if (x < 0)
        return (-1);
      v += x;
      bytes  -= x;
      buffer += x;
    }
  x = read(f,buffer,bytes);
  if (x < 0)
    return (-1);
  return (v+x);
}


static void Double_Up(Kmer_Stream *T, int nbits, int nblocks, char *output)
{ int kbyte  = T->kbyte;
  int tbyte  = T->tbyte;
  int pbyte  = T->pbyte;
  int ibyte  = T->ibyte;

  int     nbyte;
  uint32  nshift, tmask;
  int     kb1, kshift, lshift;
  uint8  *ent, *alt;

  O_Block *block;
  int64    blen;

  int64    max_el, sum_el;
  uint8   *array;
  uint8   *sarray, *tarray;
  int     *bytes;

  int      nid;
  int      ixlen;
  int64   *prefix;
  char    *root, *path;

  int      i;

  (void) print_seq;

  setup_comp_table();
  setup_fmer_table();

  path = PathTo(output);
  root = Root(output,".ktab");

  kshift = (8 - 2*(T->kmer & 0x3)) & 0x7;
  lshift = 8-kshift;
  tmask  = (0xff << kshift) & 0xff;
  kb1    = kbyte - 1;

#ifdef DEBUG
  printf("kmer = %d kbyte = %d kshift = %d lshift = %d\n",T->kmer,T->kbyte,kshift,lshift);
#endif

  nshift = 8 - (nbits & 0x7);  //  Mask block #
  nbyte  = (nbits-1)/8 + 1;    //  # of bytes to encode block number

  //  Allocatte nblock output buffers for distribution

  blen  = tbyte*BUFFER_SIZE;
  block = Malloc(sizeof(O_Block)*nblocks,"Allocating buckets");
  if (block == NULL)
    exit (1);
  block->buff = Malloc(blen*nblocks,"Allocating buckets");
  if (block->buff == NULL)
    exit (1);

  if (VERBOSE)
    fprintf(stderr,"Distributing %d-mers to %d bucket files at %s\n",
                  T->kmer,nblocks,SORT_PATH);

  for (i = 0; i < nblocks; i++)
    { block[i].nels = 0;
      block[i].fid  = open(Catenate(SORT_PATH,"/",root,Numbered_Suffix(".U.",i,"")),
                           O_CREAT|O_TRUNC|O_WRONLY,S_IRWXU);
      block[i].bptr = 0;
      block[i].buff = block->buff + i*blen;
    }

#ifdef DEBUG
  printf("nbits %d nblocks = %d (%x,%d)\n",nbits,nblocks,nmask,nbyte);
#endif

  //  Stream table T shipping each entry to appropriate file along with its complement
  //    (formed in 'alt') if not a palindrome.  Keep track of first bytes for each block
  //    and the # of elements for the ensuing sorts

  ent = Current_Entry(T,NULL); 
  alt = Current_Entry(T,NULL); 
  for (First_Kmer_Entry(T); T->csuf != NULL; Next_Kmer_Entry(T))
    { int      i, j, id;
      uint32   x, e0, e1, a0;
      O_Block *b;

      Current_Entry(T,ent); 

      x = 0;
      for (i = 0; i < nbyte; i++)
        x = (x << 8 | ent[i]);
      x >>= nshift;

#ifdef DEBUG
      print_seq(ent,T->kmer);
      printf(" %0*x\n",2*nbyte,x);
#endif

      b = block+x;
      if (b->bptr >= blen)
        { if (write(b->fid,b->buff,blen) < 0)
            { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",
                             Prog_Name,Catenate(SORT_PATH,"/",root,Numbered_Suffix(".U.",x,"")));
              exit (1);
            }
          b->bptr = 0;
        }
      memcpy(b->buff+b->bptr,ent,tbyte);
      b->bptr += tbyte;
      b->nels += 1;
    
      id = 1;
      e1 = 0;
      for (i = 0, j = kb1; i <= kb1; i++, j--)
        { e0 = ent[i];
          alt[j] = a0 = comp[(e0 >> kshift) | ((e1 << lshift) & 0xff)]; 
          if (a0 != ent[j])
            id = 0;
          e1 = e0;
        }
      alt[kb1] &= tmask;

#ifdef DEBUG
      print_seq(alt,T->kmer);
#endif

      if (id)
        {
#ifdef DEBUG
          printf(" ID\n");
#endif
          continue;
        }

      x = 0;
      for (i = 0; i < nbyte; i++)
        x = (x << 8 | alt[i]);
      x >>= nshift;

#ifdef DEBUG
      printf(" %0*x C\n",2*nbyte,x);
#endif

      *COUNT_PTR(alt) = *COUNT_PTR(ent);
    
      b = block+x;
      if (b->bptr >= blen)
        { if (write(b->fid,b->buff,blen) < 0)
            { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",
                             Prog_Name,Catenate(SORT_PATH,"/",root,Numbered_Suffix(".U.",x,"")));
              exit (1);
            }
          b->bptr = 0;
        }
      memcpy(b->buff+b->bptr,alt,tbyte);
      b->bptr += tbyte;
      b->nels += 1;
    }

  max_el = sum_el = 0;
  for (i = 0; i < nblocks; i++)
    { if (block[i].bptr > 0)
        { if (write(block[i].fid,block[i].buff,block[i].bptr) < 0)
            { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",
                             Prog_Name,Catenate(SORT_PATH,"/",root,Numbered_Suffix(".U.",i,"")));
              exit (1);
            }
        }
      close(block[i].fid);
      if (block[i].nels > max_el)
        max_el = block[i].nels;
      sum_el += block[i].nels;
    }
  free(block->buff);

  nid = open(Catenate(path,"/.",root,".ktab.1"),O_CREAT|O_TRUNC|O_WRONLY,S_IRWXU);
  if (nid == -1)
    { fprintf(stderr,"\n%s: Cannot open external file %s for writing\n",
                     Prog_Name,Catenate(path,"/.",root,".ktab.1"));
      exit (1);
    }
  ixlen = 1 << 8*ibyte;

  { int x;

    x = 0;
    if (write(nid,&(T->kmer),sizeof(int)) < 0)
      x = -1;
    if (write(nid,&sum_el,sizeof(int64)) < 0)
      x = -1;
    if (x < 0)
      { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",
                       Prog_Name,Catenate(path,"/.",root,".ktab.1"));
        exit (1);
      }
  }

  prefix = Malloc(sizeof(int64)*ixlen,"Allocating output table");
  array  = Malloc(2*tbyte*max_el,"ALlocating sort vectors");
  bytes  = Malloc(sizeof(int)*(kbyte+1),"Allocating sort vectors");
  if (prefix == NULL || array == NULL || bytes == NULL)
    exit (1);

  bzero(prefix,sizeof(int64)*ixlen);

  for (i = 0; i < kbyte; i++)
    bytes[i] = kbyte-(i+1);
  bytes[kbyte] = -1;

  for (i = 0; i < nblocks; i++)
    { int    fid, u;
      uint8 *p, *q, *t;
      uint32 x;

      sarray = array;
      tarray = array+max_el*tbyte;

      fid = open(Catenate(SORT_PATH,"/",root,Numbered_Suffix(".U.",i,"")),O_RDONLY),

      big_read(fid,sarray,tbyte*block[i].nels);

      unlink(Catenate(SORT_PATH,"/",root,Numbered_Suffix(".U.",i,"")));

      if (VERBOSE)
        fprintf(stderr,"Sorting %lld %d-mers in bucket %d\n",
                        block[i].nels,T->kmer,i+1);

      sarray = LSD_Sort(block[i].nels,sarray,tarray,tbyte,bytes);
      if (sarray != array)
        tarray = array;

      if (VERBOSE)
        fprintf(stderr,"Writing sorted %d-mers to output\n",T->kmer);

      t = tarray;
      q = sarray + tbyte*block[i].nels;
      for (p = sarray; p < q; p += pbyte)
        { x = 0;
          for (u = 0; u < ibyte; u++)
            x = (x << 8) | *p++;

          prefix[x] += 1;

          memcpy(t,p,pbyte);
          t += pbyte;
        }

      if (big_write(nid,tarray,pbyte*block[i].nels) < 0)
        { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",
                         Prog_Name,Catenate(path,"/.",root,".ktab.1"));
          exit (1);
        }
    }

  free(bytes);
  free(array);
  close(nid);

  //  Turn index counts to index offsets and create stub file

  { int x;
    int one = 1;

    for (x = 1; x < ixlen; x++)
      prefix[x] += prefix[x-1];

    nid = open(Catenate(path,"/",root,".ktab"),O_CREAT|O_TRUNC|O_WRONLY,S_IRWXU);
    if (nid == -1)
      { fprintf(stderr,"\n%s: Cannot open external file %s for writing\n",
                       Prog_Name,Catenate(path,"/",root,".ktab"));
        exit (1);
      }

    x = 0;
    if (write(nid,&(T->kmer),sizeof(int)) < 0)
      x = -1;
    if (write(nid,&one,sizeof(int)) < 0)
      x = -1;
    if (write(nid,&(T->minval),sizeof(int)) < 0)
      x = -1;
    if (write(nid,&ibyte,sizeof(int)) < 0)
      x = -1;
    if (write(nid,prefix,sizeof(int64)*ixlen) < 0)
      x = -1;
    if (x < 0)
      { fprintf(stderr,"%s: Cannot write to %s.  Enough disk space?\n",
                       Prog_Name,Catenate(path,"/",root,".ktab"));
        exit (1);
      }

    close(nid);
  }

  free(prefix);
  free(root);
  free(path);
}


/****************************************************************************************
 *
 *  Main
 *
 *****************************************************************************************/

int main(int argc, char *argv[])
{ Kmer_Stream *T;
  int          nbits, nblocks;

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    (void) flags;

    ARG_INIT("Symmex");

    SORT_PATH = "/tmp";
    NTHREADS  = 4;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
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

    if (argc != 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, output statistics as proceed.\n");
        fprintf(stderr,"      -T: Use -T threads.\n");
        fprintf(stderr,"      -P: Place all temporary files in directory -P.\n");
        exit (1);
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
                { fprintf(stderr,"\n%s: -P option: . not followed by /\n",Prog_Name);
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
        { fprintf(stderr,"\n%s: -P option: cannot open directory %s\n",Prog_Name,SORT_PATH);
          exit (1);
        }
      closedir(dirp);
    }
  }

  T = Open_Kmer_Stream(argv[1]);

  nblocks = T->nels / ((0x100000000 / T->tbyte));

#ifdef DEBUG
  printf("%d = %lld / (8B / %d)\n",nblocks,T->nels,T->tbyte);
#endif

  nbits = 1;
  while (nblocks > 1)
    { nblocks >>= 1;
      nbits  += 1;
    }

  nblocks = (0x1 << nbits);

  Double_Up(T,nbits,nblocks,argv[2]);

  Free_Kmer_Stream(T);

  exit (0);
}
