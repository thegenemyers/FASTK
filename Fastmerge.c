/*******************************************************************************************
 *
 *  Merges tables, histograms, and profiles produced by independent HPC runs on sub-parts of a data set
 *
 *  Author:  Gene Myers
 *  Date  :  Aug. 20, 2021
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <pthread.h>

#undef  DEBUG
#undef  DEBUG_THREADS
#undef  DEBUG_TRACE

#include "libfastk.h"

static char *Usage = " [-T<int(4)>] <out_root> <sources_root>[.ktab|.prof|.hist] ...";

static int NTHREADS;


/****************************************************************************************
 *
 *  Streaming eval
 *
 *****************************************************************************************/

typedef struct
  { int           tid;
    Kmer_Stream **S;
    int           narg;
    FILE         *out;
    int64        *prefx;
    int64        *begs;
    int64        *ends;
  } TP;

static inline int mycmp(uint8 *a, uint8 *b, int n)
{ while (n-- > 0)
    { if (*a++ != *b++)
        return (a[-1] < b[-1] ? -1 : 1);
    }
  return (0);
}

static void *merge_thread(void *args)
{ TP *parm = (TP *) args;
  int           tid   = parm->tid;
  int           ntabs = parm->narg;
  Kmer_Stream **S     = parm->S;
  Kmer_Stream **T;
  int64        *begs  = parm->begs;
  int64        *ends  = parm->ends;
  int64        *prefx = parm->prefx;
  FILE         *out   = parm->out;

  int hbyte = S[0]->hbyte;
  int kbyte = S[0]->kbyte;
  int kmer  = S[0]->kmer;

  int64   nels;
  uint8 **ent, *bst;
  int     itop, *in, cnt;
  uint16  scnt;
  int     c, x;

#ifdef DEBUG_TRACE
  char *buffer;
#endif

  in  = Malloc(sizeof(int)*ntabs,"Allocating thread working memory");
  ent = Malloc(sizeof(uint8 *)*ntabs,"Allocating thread working memory");

#ifdef DEBUG_THREADS
  printf("Doing %d:",tid);
  for (c = 0; c < ntabs; c++)
    printf(" [%lld-%lld]",begs[c],ends[c]);
  printf("\n");
#endif

  if (tid != 0)
    { T = Malloc(sizeof(Kmer_Stream *)*ntabs,"Allocating thread working memory");
      for (c = 0; c < ntabs; c++)
        T[c] = Clone_Kmer_Stream(S[c]);
    }
  else
    T = S;

#ifdef DEBUG_TRACE
  buffer = Current_Kmer(T[0],NULL);
#endif

   nels = 0;
   fwrite(&kmer,sizeof(int),1,out);
   fwrite(&nels,sizeof(int64),1,out);

  for (c = 0; c < ntabs; c++)
    GoTo_Kmer_Index(T[c],begs[c]);

  for (c = 0; c < ntabs; c++)
    ent[c] = Current_Entry(T[c],NULL);

  while (1)
    { for (c = 0; c < ntabs; c++)
        if (T[c]->cidx < ends[c])
          break;
      if (c >= ntabs)
        break;
      itop  = 1;
      in[0] = c;
      bst = ent[c];
      for (c++; c < ntabs; c++)
        { if (T[c]->cidx >= ends[c])
            continue;
          x = mycmp(ent[c],bst,kbyte);
          if (x == 0)
            in[itop++] = c;
          else if (x < 0)
            { itop  = 1;
              in[0] = c;
              bst = ent[c];
            }
        }

      cnt = 0;
      for (c = 0; c < itop; c++)
        cnt += Current_Count(T[in[c]]);
      if (cnt > 0x7fff)
        scnt = 0x7fff;
      else
        scnt = cnt;

#ifdef DEBUG_TRACE
      for (c = 0; c < itop; c++)
        { x = in[c];
          printf(" %d: %s %5d",x,Current_Kmer(T[x],buffer),Current_Count(T[x]));
        }
      printf("\n");
#endif

      fwrite(bst+3,hbyte,1,out);
      fwrite(&scnt,sizeof(uint16),1,out);
      x = (bst[0] << 16) | (bst[1] << 8) | bst[2];
      prefx[x] += 1;
      nels += 1;

      for (c = 0; c < itop; c++)
        { Kmer_Stream *t;

          x = in[c];
          t = T[x];
          Next_Kmer_Entry(t);
          if (t->csuf != NULL)
            Current_Entry(t,ent[x]);
        }
    }

   rewind(out);
   fwrite(&kmer,sizeof(int),1,out);
   fwrite(&nels,sizeof(int64),1,out);
   fclose(out);

  for (c = 0; c < ntabs; c++)
    free(ent[c]);

  if (tid != 0)
    { for (c = 0; c < ntabs; c++) 
        Free_Kmer_Stream(T[c]);
      free(T);
    }

  free(ent);
  free(in);

  return (NULL);
}


/****************************************************************************************
 *
 *  Main
 *
 *****************************************************************************************/

int main(int argc, char *argv[])
{ int           narg, kmer;
  char         *Opath, *Oroot;
  Kmer_Stream **S;
  int           DO_HIST;
  int           DO_TABLE;
  int           DO_PROF;

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    (void) flags;

    ARG_INIT("Fastmerge");

    NTHREADS = 4;
    DO_HIST  = 0;
    DO_TABLE = 0;
    DO_PROF  = 0;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
        }   
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc < 4)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -T: Use -T threads.\n");
        exit (1);
      } 

    if (strcmp(argv[1]+(strlen(argv[1])-5),".ktab") == 0
       || strcmp(argv[1]+(strlen(argv[1])-5),".hist") == 0
       || strcmp(argv[1]+(strlen(argv[1])-5),".prof") == 0)
      { fprintf(stderr,"%s: Output name %s has a FastK suffix ending??\n",Prog_Name,argv[1]);
        exit (1);
      }
    Opath = PathTo(argv[1]);
    Oroot = Root(argv[1],".ktab");

    { FILE *f;
      char *r;

      r = Root(argv[2],".hist");
      f = fopen(Catenate(r,".hist","",""),"r");
      if (f != NULL)
        { DO_HIST = 1;
          fclose(f);
        } 
      f = fopen(Catenate(r,".ktab","",""),"r");
      if (f != NULL)
        { DO_TABLE = 1;
          fclose(f);
        } 
      f = fopen(Catenate(r,".prof","",""),"r");
      if (f != NULL)
        { DO_PROF = 1;
          fclose(f);
        } 

      if (DO_HIST + DO_TABLE + DO_PROF == 0)
        { fprintf(stderr,"%s: There are no FastK objects with root %s !?\n",Prog_Name,r);
          exit (1);
        }

      free(r);
    }
  }   

  { int c;

    narg = argc-2;

    S = Malloc(sizeof(Kmer_Stream *)*narg,"Allocating table pointers");
    if (S == NULL)
      exit (1);

    kmer = 0;
    for (c = 2; c <= narg; c++)
      { Kmer_Stream *s = Open_Kmer_Stream(argv[c]);
        if (s == NULL)
          { fprintf(stderr,"%s: Cannot open table %s\n",Prog_Name,argv[c]);
            exit (1);
          }
        if (c == 2)
          kmer = s->kmer;
        else
          { if (s->kmer != kmer)
              { fprintf(stderr,"%s: K-mer tables do not involve the same K\n",Prog_Name);
                exit (1);
              }
          }
        S[c-2] = s;
      }
  }

  { int64     range[NTHREADS+1][narg];
#ifndef DEBUG_THREADS
    pthread_t threads[NTHREADS];
#endif
    TP        parm[NTHREADS];
    FILE     *out[NTHREADS];
    int64    *prefx;
    int       ixlen = 0;
    char     *seq;
    uint8    *ent;
    int       t, a, i;
    int64     p;

    ixlen = 0x1000000;
    prefx = Malloc(sizeof(int64)*ixlen,"Allocating prefix table");
    bzero(prefx,sizeof(int64)*ixlen);

    for (t = 0; t < NTHREADS; t++)
      out[t] = fopen(Catenate(Opath,"/.",Oroot,
                              Numbered_Suffix(".ktab.",t+1,"")),"w");

    for (a = 0; a < narg; a++)
      { range[0][a] = 0;
        range[NTHREADS][a] = S[a]->nels;
      }

    seq = Current_Kmer(S[0],NULL);
    ent = Current_Entry(S[0],NULL);
    for (t = 1; t < NTHREADS; t++)
      { p = (S[0]->nels*t)/NTHREADS; 
        GoTo_Kmer_Index(S[0],p);
#ifdef DEBUG
        printf("\n%d: %0*x\n",t,2*S[0]->ibyte,S[0]->cpre);
        printf(" %lld: %s\n",p,Current_Kmer(S[0],seq));
#endif
        ent = Current_Entry(S[0],ent);                //  Break at prefix boundaries
        for (i = S[0]->ibyte; i < S[0]->kbyte; i++)
          ent[i] = 0;
        for (a = 0; a < narg; a++)
          { GoTo_Kmer_Entry(S[a],ent);
#ifdef DEBUG
            printf(" %lld: %s\n",S[a]->cidx,Current_Kmer(S[a],seq));
#endif
            range[t][a] = S[a]->cidx;
          }
      }
    free(seq);

    for (t = 0; t < NTHREADS; t++)
      { parm[t].tid   = t;
        parm[t].S     = S;
        parm[t].narg  = narg;
        parm[t].begs  = range[t];
        parm[t].ends  = range[t+1];
        parm[t].prefx = prefx;
        parm[t].out   = out[t];
      }

#ifdef DEBUG_THREADS
    for (t = 0; t < NTHREADS; t++)
      merge_thread(parm+t);
#else
    for (t = 1; t < NTHREADS; t++)
      pthread_create(threads+t,NULL,merge_thread,parm+t);
    merge_thread(parm);
    for (t = 1; t < NTHREADS; t++)
      pthread_join(threads[t],NULL);
#endif

    { int one = 1;
      int three = 3;

      FILE  *f   = fopen(Catenate(Opath,"/",Oroot,".ktab"),"w");
      int64 *prf = prefx;

      fwrite(&kmer,sizeof(int),1,f);
      fwrite(&NTHREADS,sizeof(int),1,f);
      fwrite(&one,sizeof(int),1,f);
      fwrite(&three,sizeof(int),1,f);

      for (i = 1; i < ixlen; i++)
        prf[i] += prf[i-1];

      fwrite(prf,sizeof(int64),ixlen,f);
      fclose(f);

      free(out);
      free(prefx);
    }
  }

  { int c;

    for (c = 0; c < narg; c++)
      Free_Kmer_Stream(S[c]);
    free(S);
  }

  free(Oroot);
  free(Opath);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
