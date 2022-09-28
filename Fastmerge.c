/*******************************************************************************************
 *
 *  Merges tables, histograms, and profiles produced by independent HPC runs on sub-parts
 *     of a data set
 *
 *  Author:  Gene Myers
 *  Date  :  Sep. 20, 2021
 *  Date  :  Mar. 31, 2022  eq-priority queue, tailored prefix size, thread parts
 *  Date  :  June. 1, 2022  caching and slices for really big data sets
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <pthread.h>
#include <sys/resource.h>

#undef    DEBUG
#undef    DEBUG_THREADS
#undef    DEBUG_TRACE
#undef    DEBUG_PROF

#include "libfastk.h"

static char *Usage[] = { "[-ht] [-T<int(4)>] [#<int(1)>] [-P<dir(/tmp)>] [-S<N:int>of<D:int>]",
                         "<target> <source>[.hist|.ktab] ..."
                       };

static int   NTHREADS;
static int   NPARTS;
static char *OPATH;
static char *OROOT;
static int   PIVOT;
static char *SORT_PATH;

static int   TABLE_ERROR;

#define XFER_SIZE   0x8000000ll   // 128MB transfer buffer

  //  Special hooks into libfastk.c to support table caching

extern int          Open_Kmer_Cache(Kmer_Stream *, int64, int64, int, int, char *, uint8 *, int);
extern Kmer_Stream *Clone_Kmer_Cache(Kmer_Stream *, int64, int64, int);
extern void         Free_Kmer_Cache(Kmer_Stream *, int);


/****************************************************************************************
 *
 *  Streaming threaded merge of k-mer tables
 *
 *****************************************************************************************/

typedef struct
  { int           tid;     //  thread id
    Kmer_Stream **S;       //  array of open streams/caches (not clones)
    int           narg;    //  number of 
    int64        *prefx;   //  prefix index for this thread
    int           ibyte;   //  # of prefix bytes
    int64        *bidx;    //  merge [bidx,eidx)[c] for c in [0,narg)
    int64        *eidx;
    int          *bpre;    //  prefix for bidx[c] is bpre[c]
    int64        *fidx;    //  first index of slice is fidx[c]
    int64        *hist;    //  histogram count array
    int           dotab;   //  make table?
  } TP;

static inline int mycmp(uint8 *a, uint8 *b, int n)
{ while (n-- > 0)
    { if (*a++ != *b++)
        return (a[-1] < b[-1] ? -1 : 1);
    }
  return (0);
}    

  //  Heap of input buffer pointers ordering on table entry at ent[x]

static int KBYTE;

#define EQ_NONE  0x0   
#define EQ_RGHT  0x1   //  node value is equal to that of its right child
#define EQ_LEFT  0x2   //  node value is equal to that of its left child
#define EQ_BOTH  0x3

static void reheap(int s, int *heap, int *equb, int hsize, uint8 **ent)
{ int       c, l, r;
  int       hs, hr, hl;
  uint8    *es, *er, *el;
  int       s1, s2;
  
  c   = s;
  hs  = heap[s];
  es  = ent[hs];
  while ((l = (c<<1)) <= hsize)
    { r  = l+1;
      hl = heap[l];
      el = ent[hl];
      if (r > hsize)
        s1 = 1;
      else
        { hr = heap[r];
          er = ent[hr];
          s1 = mycmp(er,el,KBYTE);
        }
      if (s1 > 0)
        { s2 = mycmp(es,el,KBYTE);
          if (s2 > 0)
            { heap[c] = hl;
              equb[c] = (equb[l] ? EQ_LEFT : EQ_NONE);
              c = l;
            }
          else if (s2 == 0)
            { heap[c] = hs;
              equb[c] = EQ_LEFT;
              return;
            }
          else
            break;
        }
      else if (s1 == 0)
        { s2 = mycmp(es,el,KBYTE);
          if (s2 > 0)
            { heap[c] = hl;
              equb[c] = (equb[l] ? EQ_BOTH : EQ_RGHT);
              c = l;
            }
          else if (s2 == 0)
            { heap[c] = hs;
              equb[c] = EQ_BOTH;
              return;
            }
          else
            break;
        }
      else
        { s2 = mycmp(es,er,KBYTE);
          if (s2 > 0)
            { heap[c] = hr;
              equb[c] = (equb[r] ? EQ_RGHT : EQ_NONE);
              c = r;
            }
          else if (s2 == 0)
            { heap[c] = hs;
              equb[c] = EQ_RGHT;
              return;
            }
          else
            break;
        }
    }
  heap[c] = hs;
  equb[c] = EQ_NONE;
}

  //  make a list of all the nodes containing equal elements

static int next_group(int node, int *equb, int gtop, int *group)
{ if (equb[node] & EQ_RGHT)
    gtop = next_group(2*node+1,equb,gtop,group);
  if (equb[node] & EQ_LEFT)
    gtop = next_group(2*node,equb,gtop,group);
  group[++gtop] = node;
  return (gtop);
}

  //  table thread (see definition of TP)

static void *table_thread(void *args)
{ TP *parm = (TP *) args;
  int           tid   = parm->tid;
  int           ntabs = parm->narg;
  Kmer_Stream **S     = parm->S;
  int64        *fidx  = parm->fidx;
  int64        *bidx  = parm->bidx;
  int64        *eidx  = parm->eidx;
  int          *bpre  = parm->bpre;
  int64        *prefx = parm->prefx;
  int           ibyte = parm->ibyte;
  int           dotab = parm->dotab;

  int kbyte = S[0]->kbyte;
  int kmer  = S[0]->kmer;
  int hbyte = kbyte-ibyte;

  Kmer_Stream **T, *t, *p;
  char   *nbuf;
  FILE   *out;
  int64   nels, pend;
  int     npr;
  int64  *hist;
  int     hsize, *heap, *equb;
  uint8 **ent, *best;
  int    *grp, gtop;
  int64  *cnt, icnt;
  uint16  scnt;
  int     c, x, i, k;

#ifdef DEBUG_TRACE
  char *buffer;
#endif

  //  Allocate histogram, heap, and name buffer

  hist = Malloc(sizeof(int64)*0x8001,"Allocating histogram");
  bzero(hist,sizeof(int64)*0x8001);
  hist -= 1;

  cnt  = Malloc(sizeof(int64)*ntabs,"Allocating thread working memory");
  grp  = Malloc(sizeof(int64)*ntabs,"Allocating thread working memory");
  ent  = Malloc(sizeof(uint8 *)*ntabs,"Allocating thread working memory");
  heap = Malloc(sizeof(int)*(ntabs+1),"Allocating thread working memory");
  equb = Malloc(sizeof(int)*(ntabs+1),"Allocating thread working memory");
  nbuf = Malloc(strlen(OPATH)+strlen(OROOT)+20,"Allocating thread working memory");

#ifdef DEBUG_THREADS
  printf("\nDoing %d: pivot %d\n",tid,PIVOT);
  for (c = 0; c < ntabs; c++)
    printf("  %2d: [%lld-%lld]\n",c,bidx[c],eidx[c]);
#endif

  //  Setup tables if tid=0, otherwise make clones.  Set start position if not cached

  if (SORT_PATH == NULL)
    { if (tid != 0)
        { T = Malloc(sizeof(Kmer_Stream *)*ntabs,"Allocating thread working memory");
          if (T == NULL)
            { TABLE_ERROR = 1;
              return (NULL);
            }
          for (c = 0; c < ntabs; c++)
            T[c] = Clone_Kmer_Stream(S[c]);
        }
      else
        T = S;

      for (c = 0; c < ntabs; c++)
        GoTo_Kmer_Index(T[c],bidx[c]);
    }
  else
    { if (tid != 0)
        { T = Malloc(sizeof(Kmer_Stream *)*ntabs,"Allocating thread working memory");
          if (T == NULL)
            { TABLE_ERROR = 1;
              return (NULL);
            }
          for (c = 0; c < ntabs; c++)
            T[c] = Clone_Kmer_Cache(S[c],fidx[c],bidx[c],bpre[c]);
        }
      else
        T = S;
    }

  p = T[PIVOT];

#ifdef DEBUG_TRACE
  buffer = Current_Kmer(T[0],NULL);
#endif

  //  Start output of first table part, pend is the index in the pivot table to end this part

  if (dotab)
    { npr  = 1;
      pend = bidx[PIVOT] + ((eidx[PIVOT] - bidx[PIVOT]) * npr) / NPARTS;

      sprintf(nbuf,"%s/.%s.ktab.%d",OPATH,OROOT,tid*NPARTS+npr);
      out = fopen(nbuf,"w");
      if (out == NULL)
        { fprintf(stderr,"%s: Cannot open %s for writing\n",Prog_Name,nbuf);
          TABLE_ERROR = 1;
          return (NULL);
        }

#ifdef DEBUG_THREADS
      printf("  Making %s to %d:%lld\n",nbuf,npr,pend);
#endif

      nels = 0;
      if (fwrite(&kmer,sizeof(int),1,out) < 1)
        goto io_error;
      if (fwrite(&nels,sizeof(int64),1,out) < 1)
        goto io_error;
    }

  //  Init heap, element for an exhausted table has value 0xfff... with count 0

  hsize = ntabs;
  for (c = 0; c < ntabs; c++)
    { ent[c] = Current_Entry(T[c],NULL);
      if (T[c]->cidx >= eidx[c])
        { for (k = 0; k < kbyte; k++)
            ent[c][k] = 0xff;
          *((uint16 *) (ent[c] + kbyte)) = 0;
          hsize -= 1;
        }
    }

  KBYTE = kbyte;
  for (c = 1; c <= ntabs; c++)
    { heap[c] = c-1;
      equb[c] = EQ_NONE;
    }

  if (ntabs > 3)
    for (x = ntabs/2; x >= 1; x--)
      reheap(x,heap,equb,ntabs,ent);

  //  While the input tables are not exhausted, get the next =-group and process ...

  while (hsize > 0)
    { gtop = next_group(1,equb,-1,grp);
      best = ent[heap[1]];

      //  compute the sum of the counts of the equal elements in grp[0..gtop]

      icnt = (cnt[gtop] = *((uint16 *) (best + kbyte)));
      for (x = 0; x < gtop; x++)
        icnt += (cnt[x] = *((uint16 *) (ent[heap[grp[x]]] + kbyte)));

      //  be careful to handle overflow counts correctly

      if (icnt > 0x7fff)
        { scnt = 0x7fff;
          hist[0x7fff] += 1;
          for (x = 0; x <= gtop; x++)
            if (cnt[x] < 0x7fff)
              hist[0x8001] += cnt[x];
        }
      else
        { scnt = icnt;
          hist[icnt] += 1;
        }

#ifdef DEBUG_TRACE
      printf("%lld: ",icnt);
      printf("%s: %5d\n",Current_Kmer(T[heap[1]],buffer),scnt);
      for (x = 0; x <= gtop; x++)
        printf("    %2d: %2d: %2d: %10lld\n",x,grp[x],heap[grp[x]],cnt[x]);
      fflush(stdout);

/*
      for (i = 1; i <= ntabs; i++)      // Show heap
        { c = heap[i];
          if (T[c]->cidx < ends[c])
            printf(" %2d: %2d(%d) -> %s\n",i,c,equb[i],Current_Kmer(T[c],buffer));
          else
            printf(" %2d: %2d(%d) -> ttt...\n",i,heap[i],equb[i]);
        }
*/
#endif

      //  output table entry with count.  Stop and start a new table part if pend reached

      if (dotab)
        { if (fwrite(best+ibyte,hbyte,1,out) < 1)
            goto io_error;
          if (fwrite(&scnt,sizeof(uint16),1,out) < 1)
            goto io_error;

          x = best[0];
          for (i = 1; i < ibyte; i++)
            x = (x << 8) | best[i];
          prefx[x] += 1;
          nels += 1;

          if (p->cidx >= pend && npr < NPARTS)
            { rewind(out);
              fwrite(&kmer,sizeof(int),1,out);
              fwrite(&nels,sizeof(int64),1,out);
              fclose(out);

              npr += 1;
              pend = bidx[PIVOT] + ((eidx[PIVOT] - bidx[PIVOT]) * npr) / NPARTS;
              sprintf(nbuf,"%s/.%s.ktab.%d",OPATH,OROOT,tid*NPARTS+npr);
              out = fopen(nbuf,"w");
              if (out == NULL)
                { fprintf(stderr,"%s: Cannot open %s for writing\n",Prog_Name,nbuf);
                  TABLE_ERROR = 1;
                  return (NULL);
                }

#ifdef DEBUG_THREADS
              printf("  Making %s to %d:%lld\n",nbuf,npr,pend);
#endif

              nels = 0;
              if (fwrite(&kmer,sizeof(int),1,out) < 1)
                goto io_error;
              if (fwrite(&nels,sizeof(int64),1,out) < 1)
                goto io_error;
            }
        }

      //  refresh heap with next entries

      for (x = 0; x <= gtop; x++)
        { i = grp[x];
          c = heap[i];
          t = T[c];
          Next_Kmer_Entry(t);
          if (t->cidx < eidx[c])
            Current_Entry(t,ent[c]);
          else
            { for (k = 0; k < kbyte; k++)
                ent[c][k] = 0xff;
              *((uint16 *) (ent[c] + kbyte)) = 0;
              hsize -= 1;
            }
          reheap(i,heap,equb,ntabs,ent);
        }
    }

  //  finish current table part

  if (dotab)
    { rewind(out);
      fwrite(&kmer,sizeof(int),1,out);
      fwrite(&nels,sizeof(int64),1,out);
      fclose(out);
    }

  //  clean up and return histogram

  for (c = 0; c < ntabs; c++)
    free(ent[c]);

  if (tid != 0)
    { if (SORT_PATH == NULL)
        for (c = 0; c < ntabs; c++) 
          Free_Kmer_Stream(T[c]);
      else
        for (c = 0; c < ntabs; c++) 
          Free_Kmer_Cache(T[c],bpre[c]);
      free(T);
    }

  free(nbuf);
  free(equb);
  free(heap);
  free(ent);
  free(grp);
  free(cnt);

  parm->hist = hist;
  return (NULL);

io_error:
  fprintf(stderr,"%s: Could not write to file %s, out of disk space?\n",Prog_Name,nbuf);
  TABLE_ERROR = 1;
  return (NULL);
}


/****************************************************************************************
 *
 *  Main
 *
 *****************************************************************************************/

int main(int argc, char *argv[])
{ int           narg, kmer;
  int           FRAC_NUM, FRAC_DEN;
  int           DO_HIST;
  int           DO_TABLE;

  //  Process command line

  { int    i, j, k;
    int    flags[128];
    char  *eptr, *fptr;

    ARG_INIT("Fastmerge");

    NTHREADS  = 4;
    NPARTS    = 1;
    SORT_PATH = NULL;
    FRAC_NUM  = 0;
    FRAC_DEN  = 1;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("ht")
            break;
          case '#':
            ARG_POSITIVE(NPARTS,"Number of parts per thread")
            break;
          case 'P':
            SORT_PATH = argv[i]+2;
            break;
          case 'S':
            FRAC_NUM = strtol(argv[i]+2,&eptr,10);
            if (eptr > argv[i]+2 && strncmp(eptr,"of",2) == 0)
              { FRAC_DEN = strtol(eptr+2,&fptr,10);
                if (fptr > eptr+2 && *fptr == '\0')
                  { if (FRAC_DEN < 1)
                      { fprintf(stderr,"%s: Fraction denominator %d is not positive\n",
                                       Prog_Name,FRAC_DEN);
                        exit (1);
                      }
                    if (FRAC_NUM < 1 || FRAC_NUM > FRAC_DEN)
                      { fprintf(stderr,"%s: Fraction numerator %d is out of range\n",
                                       Prog_Name,FRAC_NUM);
                        exit (1);
                      }
                    FRAC_NUM -= 1;
                    break;
                  }
              }
            fprintf(stderr,"%s: Syntax of -S option invalid -S<int>of<int>\n",Prog_Name);
            exit (1);
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
        }   
      else
        argv[j++] = argv[i];
    argc = j;

    DO_HIST  = flags['h'];
    DO_TABLE = flags['t'];

    if (argc < 4)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -h: Produce a merged histogram.\n");
        fprintf(stderr,"      -t: Produce a merged k-mer table.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -T: Use -T threads.\n");
        fprintf(stderr,"      -#: Produce -# parts per thread.\n");
        fprintf(stderr,"      -P: Cache table inputs to this directory.\n");
        fprintf(stderr,"      -S: Divide into D slices and do slice N in [1,D].\n");
        exit (1);
      } 

    if (DO_HIST + DO_TABLE == 0)
      { fprintf(stderr,"%s: At least one of -h or -t must be set\n",Prog_Name);
        exit (1);
      }

    //  Get full path string for sorting subdirectory (in variable SORT_PATH)

    if (SORT_PATH != NULL)
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

    //  Remove FastK extensions from source arguments if any
    //  The target should not have one

    for (i = 1; i < argc; i++)
      { int dot = strlen(argv[i])-5;
        if (dot < 1)
          continue;
        if (strcmp(argv[i]+dot,".hist") == 0)
          argv[i][dot] = '\0';
        if (strcmp(argv[i]+dot,".ktab") == 0)
          argv[i][dot] = '\0';
        if (i == 1 && argv[1][dot] == '\0')
          { fprintf(stderr,"%s: Target name cannot have a .hist, .ktab, or .prof suffix\n",
                           Prog_Name);
            exit (1);
          }
      }

    //  The destination path and root

    OPATH = PathTo(argv[1]);
    OROOT = Root(argv[1],"");

    //  Make sure that sources have a full complement of FastK tables

    { FILE *f;
      int   has_table;
    
      narg = argc-2;
      argv += 2;

      has_table = 0;
      for (i = 0; i < narg; i++)
        { f = fopen(Catenate(argv[0],".ktab","",""),"r");
          if (f != NULL)
            { has_table += 1;
              fclose(f);
            } 
        }

      if (has_table != narg)
        { if (has_table == 0)
            fprintf(stderr,"%s: None of the sources have FastK table files?\n",Prog_Name);
          else
            fprintf(stderr,"%s: %d of the sources do not have FastK table files?\n",
                           Prog_Name,narg-has_table);
          exit (1);
        }
    }
  }   

  //  Make sure you can open max(4,(narg+1))*NTHREADS+tid files
  //    tid is typically 3 unless using valgrind or other instrumentation.

  { struct rlimit rlp;
    int           tid;
    uint64        nfiles;

    tid = open(".xxx",O_CREAT|O_TRUNC|O_WRONLY,S_IRWXU);
    close(tid);
    unlink(".xxx");

    if (narg >= 3)
      nfiles = (narg+1)*NTHREADS + tid;
    else
      nfiles = 4*NTHREADS + tid;
    
    getrlimit(RLIMIT_NOFILE,&rlp);
    if (nfiles > rlp.rlim_max)
      { fprintf(stderr,"\n%s: Cannot open %lld files simultaneously\n",Prog_Name,nfiles);
        exit(1);
      }
    rlp.rlim_cur = nfiles;
    setrlimit(RLIMIT_NOFILE,&rlp);
  }

  { Kmer_Stream **S;
    int64        *range[NTHREADS+1];
    int          *prefs[NTHREADS+1];
    TP            parm[NTHREADS];
#ifndef DEBUG_THREADS
    pthread_t     threads[NTHREADS];
#endif
    int64         tels;
    int           minval;
  
    //  Allocate table and partition vectors

    { int t;

      S   = Malloc(sizeof(Kmer_Stream *)*narg,"Allocating table pointers");
      range[0] = Malloc(sizeof(int64)*narg*(NTHREADS+1),"Allocating table partition");
      prefs[0] = Malloc(sizeof(int)*narg*(NTHREADS+1),"Allocating table partition");
      if (S == NULL || range[0] == NULL || prefs[0] == NULL)
        exit (1);
      for (t = 1; t <= NTHREADS; t++)
        { range[t] = range[t-1] + narg;
          prefs[t] = prefs[t-1] + narg;
        }
    }

    //  Read each source table header to determine pivot and output header values

    { int   c, f;
      int   smer, smin, ibyte;
      int64 nels, npiv;
      char *dir, *root, *full;
  
#ifdef DEBUG
      printf("\nSizes:\n");
#endif
      tels   = 0;
      npiv   = 0;
      kmer   = 0;
      minval = 0x10000;
      for (c = 0; c < narg; c++)
        { dir  = PathTo(argv[c]);
          root = Root(argv[c],".ktab");
          full = Malloc(strlen(dir)+strlen(root)+20,"Histogram name allocation");
          if (full == NULL)
            exit (1);
          sprintf(full,"%s/%s.ktab",dir,root);
          f = open(full,O_RDONLY);
          free(full);
          free(root);
          free(dir);
          if (f < 0)
            { fprintf(stderr,"%s: Cannot open table %s\n",Prog_Name,argv[c]);
              exit (1);
            }

          read(f,&smer,sizeof(int));
          read(f,&ibyte,sizeof(int));
          read(f,&smin,sizeof(int));
          read(f,&ibyte,sizeof(int));
          lseek(f,sizeof(int64)*((0x1<<(8*ibyte))-1),SEEK_CUR);
          read(f,&nels,sizeof(int64));

#ifdef DEBUG
          printf("  %10lld: %s\n",nels,argv[c]);
#endif

          if (smin < minval)
            minval = smin;
          if (nels > npiv)
            { npiv  = nels;
              PIVOT = c;
            }
          tels += nels;

          if (c == 0)
            kmer = smer;
          else
            { if (smer != kmer)
                { fprintf(stderr,"%s: K-mer tables do not involve the same K\n",Prog_Name);
                  exit (1);
                }
            }
          close(f);
        }
#ifdef DEBUG
      printf("%lld: %d %d\n\n",tels,kmer,PIVOT);
      fflush(stdout);
#endif
    }
  
    { int64       *prefx = NULL;
      int          ixlen = 0;
      int          ibyte = 0;

      //  Determine prefix length and allocate suitably large, zero'd vector

      if (DO_TABLE)
        { if (tels >= 0x8000000)
            { ixlen = 0x1000000;
              ibyte = 3;
            }
          else if (tels >= 0x80000)
            { ixlen = 0x10000;
              ibyte = 2;
            }
          else
            { ixlen = 0x100;
              ibyte = 1;
            }
          prefx = Malloc(sizeof(int64)*ixlen,"Allocating prefix table");
          bzero(prefx,sizeof(int64)*ixlen);
        }

      { uint8       *pivot[NTHREADS+1];
        uint8       *ent;
        int          t, i;
        int64        p;
        Kmer_Stream *ktab;
        uint8       *xferbuf;

        //  If caching will need a big transfer buffer

        if (SORT_PATH != NULL)
          { xferbuf = Malloc(XFER_SIZE,"Allocating cache buffer");
            if (xferbuf == NULL)
              exit (1);
          }

        //  Open the pivot table and determine partition for the requested slice
        //    The partition must be at a prefix boundary

        S[PIVOT] = ktab = Open_Kmer_Stream(argv[PIVOT]);
        for (t = 0; t <= NTHREADS; t++)
          { p = (ktab->nels*(FRAC_NUM*NTHREADS+t))/(NTHREADS*FRAC_DEN); 
            if (p >= ktab->nels)
              { range[t][PIVOT] = ktab->nels;
                prefs[t][PIVOT] = (1 << (8*ktab->ibyte));
                pivot[t] = NULL;
#ifdef DEBUG
                printf("\n%d: %0*x\n",t,2*ktab->ibyte,prefs[t][PIVOT]);
                printf("  %10lld/%0*x: EOT\n",ktab->nels,2*ktab->ibyte,prefs[t][PIVOT]);
                fflush(stdout);
#endif
              }
            else
              { GoTo_Kmer_Index(ktab,p);
#ifdef DEBUG
                { printf("\n%d: %0*x\n",t,2*ktab->ibyte,ktab->cpre);
                  printf(" %10lld:",p);
                  char *seq = Current_Kmer(ktab,NULL);
                  printf(" %s\n",seq);
                  free(seq);
                }
#endif
                ent = Current_Entry(ktab,NULL);
                for (i = ibyte; i < ktab->kbyte; i++)
                  ent[i] = 0;
                GoTo_Kmer_Entry(ktab,ent);
                pivot[t] = ent;
                range[t][PIVOT] = ktab->cidx;
                prefs[t][PIVOT] = ktab->cpre;
#ifdef DEBUG
                { printf("  %10lld/%0*x:",ktab->cidx,2*ktab->ibyte,ktab->cpre);
                  char *seq = Current_Kmer(ktab,NULL);
                  printf(" %s\n",seq);
                  free(seq);
                  fflush(stdout);
                }
#endif
              }
          }

        //  If caching then create the cache for the pivot slice now

        if (SORT_PATH != NULL)
          { if (Open_Kmer_Cache(ktab,range[0][PIVOT],range[NTHREADS][PIVOT],
                                prefs[0][PIVOT],prefs[NTHREADS][PIVOT],
                                SORT_PATH,xferbuf,XFER_SIZE) )
              { fprintf(stderr,"%s: Directory %s appears to be out of space\n",
                               Prog_Name,SORT_PATH);
                exit (1);
              }
          }

        //  Find corresponding partition points in all the other tables

        for (i = 0; i < narg; i++)
          if (i != PIVOT)
            { S[i] = ktab = Open_Kmer_Stream(argv[i]);
#ifdef DEBUG
              printf("\n");
#endif
              for (t = 0; t <= NTHREADS; t++)
                if (pivot[t] == NULL)
                  { range[t][i] = ktab->nels;
                    prefs[t][i] = (1 << (8*ktab->ibyte));
#ifdef DEBUG
                    printf("  %10lld/%0*x: EOT\n",ktab->nels,2*ktab->ibyte,prefs[t][i]);
                    fflush(stdout);
#endif
                  }
                else
                  { GoTo_Kmer_Entry(ktab,pivot[t]);
                    range[t][i] = ktab->cidx;
                    prefs[t][i] = ktab->cpre;
#ifdef DEBUG
                    printf("  %10lld/%0*x:",ktab->cidx,2*ktab->ibyte,ktab->cpre);
                    if (ktab->cidx < ktab->nels)
                      { char *seq = Current_Kmer(ktab,NULL);
                        printf(" %s\n",seq);
                        free(seq);
                      }
                    else
                      printf(" EOT\n");
                    fflush(stdout);
#endif
                  }

              //  If caching then create the cache for this table now

              if (SORT_PATH != NULL)
                { if (Open_Kmer_Cache(ktab,range[0][i],range[NTHREADS][i],
                                           prefs[0][i],prefs[NTHREADS][i],
                                           SORT_PATH,xferbuf,XFER_SIZE) )
                    { int j;

                      fprintf(stderr,"%s: Directory %s appears to be out of space\n",
                                     Prog_Name,SORT_PATH);
                      for (j = 0; j < i; j++)
                        Free_Kmer_Cache(S[j],prefs[0][j]);
                      if (PIVOT > i)
                        Free_Kmer_Cache(S[PIVOT],prefs[0][PIVOT]);
                      exit (1);
                    }
                }
            }

        for (t = 0; t < NTHREADS; t++)
          free(pivot[t]);
        if (SORT_PATH != NULL)
          free(xferbuf);
      }

      //  Call a thread to do the merge on each partition of the requested slice

      { int t;

        TABLE_ERROR = 0;
        for (t = 0; t < NTHREADS; t++)
          { parm[t].tid   = t;
            parm[t].S     = S;
            parm[t].narg  = narg;
            parm[t].fidx  = range[0];
            parm[t].bidx  = range[t];
            parm[t].eidx  = range[t+1];
            parm[t].bpre  = prefs[t];
            parm[t].prefx = prefx;
            parm[t].ibyte = ibyte;
            parm[t].dotab = DO_TABLE;
          }

#ifdef DEBUG_THREADS
        for (t = 0; t < NTHREADS; t++)
          table_thread(parm+t);
#else
        for (t = 1; t < NTHREADS; t++)
          pthread_create(threads+t,NULL,table_thread,parm+t);
        table_thread(parm);
        for (t = 1; t < NTHREADS; t++)
          pthread_join(threads[t],NULL);
#endif
      }

      { int c;
  
        if (SORT_PATH == NULL)
          for (c = 0; c < narg; c++)
            Free_Kmer_Stream(S[c]);
        else
          for (c = 0; c < narg; c++)
            Free_Kmer_Cache(S[c],prefs[0][c]);
        free(prefs[0]);
        free(range[0]);
        free(S);
      }

      //  Write the output table header stub (provided no error occurred)

      if (DO_TABLE & !TABLE_ERROR)
        { int nparts;
          int i, f;

          nparts = NTHREADS * NPARTS;
          for (i = 1; i < ixlen; i++)
            prefx[i] += prefx[i-1];

          f = open(Catenate(OPATH,"/",OROOT,".ktab"),O_CREAT|O_TRUNC|O_WRONLY,S_IRWXU);
          if (f < 0)
            { fprintf(stderr,"%s: Cannot open %s/%s.ktab\n",Prog_Name,OPATH,OROOT);
              TABLE_ERROR = 1;
            }
          else
            { TABLE_ERROR |= (write(f,&kmer,sizeof(int)) < 0);
              TABLE_ERROR |= (write(f,&nparts,sizeof(int)) < 0);
              TABLE_ERROR |= (write(f,&minval,sizeof(int)) < 0);
              TABLE_ERROR |= (write(f,&ibyte,sizeof(int)) < 0);
              TABLE_ERROR |= (write(f,prefx,sizeof(int64)*ixlen) < 0);
              close(f);
              if (TABLE_ERROR)
                fprintf(stderr,"%s: Cannot write to %s/%s.ktab\n",Prog_Name,OPATH,OROOT);
            }
  
          free(prefx);
        }

      //  If an error occured in any thread, remove any output table parts and quit

      if (TABLE_ERROR)
        { int   j;
          FILE *f;

          for (j = 0; j < NTHREADS*NPARTS; j++)
            { f = fopen(Catenate(OPATH,"/.",OROOT,Numbered_Suffix(".ktab",j,"")),"r");
              if (f != NULL)
                { fclose(f);
                  unlink(Catenate(OPATH,"/.",OROOT,Numbered_Suffix(".ktab",j,"")));
                }
            }
          exit (1);
        }

      //  Error free to this point, output histogram if requested

      if (DO_HIST)
        { int64 *hist, *gist;
          int    j, low, high;
          int    f, t;

          hist = parm[0].hist;
          for (t = 1; t < NTHREADS; t++)
            { gist = parm[t].hist;
              for (j = 1; j <= 0x8001; j++)
                hist[j] += gist[j];
              free(gist+1);
	    } 

          for (t = 0; t < narg; t++)
            { Histogram *H = Load_Histogram(argv[t]);
              if (H == NULL)
                { if (t == 0)
                    { fprintf(stderr,"%s: Warning: no input histograms => overflow count low\n",
                                     Prog_Name);
                      break;
                    }
                  fprintf(stderr,"%s: Cannot open histogram %s\n",Prog_Name,argv[t]);
                  exit (1);
                }
              hist[0x8001] += H->hist[0x8001];
              Free_Histogram(H);
            }
          hist[0x8000] = hist[1];

          low  = 1;
          high = 0x7fff;
          f  = open(Catenate(OPATH,"/",OROOT,".hist"),O_CREAT|O_TRUNC|O_WRONLY,S_IRWXU);
          if (f < 0)
            { fprintf(stderr,"%s: Cannot open %s/%s.ktab\n",Prog_Name,OPATH,OROOT);
              exit (1);
            }
          else
            { TABLE_ERROR |= (write(f,&kmer,sizeof(int)) < 0);
              TABLE_ERROR |= (write(f,&low,sizeof(int)) < 0);
              TABLE_ERROR |= (write(f,&high,sizeof(int)) < 0);
              TABLE_ERROR |= (write(f,hist+0x8000,sizeof(int64)) < 0);
              TABLE_ERROR |= (write(f,hist+0x8001,sizeof(int64)) < 0);
              TABLE_ERROR |= (write(f,hist+1,sizeof(int64)*0x7fff) < 0);
              close(f);
              if (TABLE_ERROR)
                { fprintf(stderr,"%s: Cannot write to %s/%s.ktab\n",Prog_Name,OPATH,OROOT);
                  exit (1);
                }
            }

          free(hist+1);
        }
    }
  }

  free(OROOT);
  free(OPATH);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
