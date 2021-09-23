/*******************************************************************************************
 *
 *  Merges tables, histograms, and profiles produced by independent HPC runs on sub-parts
 *     of a data set
 *
 *  Author:  Gene Myers
 *  Date  :  Sep. 20, 2021
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
 *  Streaming threaded merge of k-mer tables
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

  int hbyte = S[0]->kbyte-3;
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
 *  Streaming threaded merge of profile indices and data
 *
 *****************************************************************************************/

typedef struct
  { int           tid;
    char        **argv;
    int           kmer;
    int64         fidx;
    int64         nidx;
    FILE         *pout;
    FILE         *dout;
    int           fpart;
    int64         first;
    int           lpart;
    int64         last;
    int64         buffer[16384];
  } PP;

static void *prof_thread(void *args)
{ PP    *parm   = (PP *) args;
  char **argv   = parm->argv;
  FILE  *dout   = parm->dout;
  FILE  *pout   = parm->pout;
  int    fpart  = parm->fpart;
  int64  first  = parm->first;
  int    lpart  = parm->fpart;
  int64  last   = parm->last;
  int64 *buffer = parm->buffer;

  int64  q, cindex, lindex;
  int64  base, offs;
  int64  fdata, ldata;
  char  *path, *root;
  int    imer, nthreads;
  FILE  *f;
  int    c, t;

  fwrite(&(parm->kmer),sizeof(int),1,pout);
  fwrite(&(parm->fidx),sizeof(int64),1,pout);
  fwrite(&(parm->nidx),sizeof(int64),1,pout);

  base = -1;
  for (c = fpart; c <= lpart; c++) 
    { path = PathTo(argv[c]);
      root = Root(argv[c],".prof");

      f = fopen(Catenate(path,"/",root,".prof"),"r");
      if (f == NULL)
        { fprintf(stderr,"%s: Cannot open profile %s\n",Prog_Name,argv[c]);
          exit (1);
        }
      fread(&imer,sizeof(int),1,f);
      fread(&nthreads,sizeof(int),1,f);
      fclose(f);

      for (t = 1; t <= nthreads; t++)
        { f = fopen(Catenate(path,"/.",root,Numbered_Suffix(".pidx.",t,"")),"r");
          if (f == NULL)
            { fprintf(stderr,"%s: Cannot open profile index for %s\n",Prog_Name,argv[c]);
              exit (1);
            }
          fread(&imer,sizeof(int),1,f);
          fread(&cindex,sizeof(int64),1,f);
          fread(&lindex,sizeof(int64),1,f);
          lindex += cindex;
          if (c == fpart && base < 0)
            { if (lindex <= first) 
                { fclose(f);
                  break;
                }
              if (cindex < first)
                fseeko(f,sizeof(int64)*(first-cindex),SEEK_CUR);
              fread(&offs,sizeof(int64),1,f);
              base = offs;
              cindex = first;
            }
          else
            fread(&offs,sizeof(int64),1,f);
          fdata = offs;
          if (c == lpart && lindex >= last)
            { lindex = last;
              t = nthreads;
            }
          for (q = cindex; q < lindex; q++)
            { offs -= base;
              fwrite(&offs,sizeof(int64),1,pout);
              fread(&offs,sizeof(int64),1,f);
            }
          ldata = offs;
          base = -(offs-base);
          fclose(f);

          f = fopen(Catenate(path,"/.",root,Numbered_Suffix(".prof.",t,"")),"r");
          if (f == NULL)
            { fprintf(stderr,"%s: Cannot open profile data for %s\n",Prog_Name,argv[c]);
              exit (1);
            }
          if (fdata > 0)
            fseeko(f,sizeof(int64)*fdata,SEEK_SET);
          for (q = ldata-fdata; q >= 16384; q -= 16384)
            { fread(buffer,sizeof(int64)*16384,1,f);
              fwrite(buffer,sizeof(int64)*16384,1,dout);
            }
          if (q > 0)
            { fread(buffer,sizeof(int64)*q,1,f);
              fwrite(buffer,sizeof(int64)*q,1,dout);
            }
          fclose(f);
        }

      free(root);
      free(path);
    }
  fwrite(&offs,sizeof(int64),1,pout);

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
    
      narg = argc-2;
      argv += 2;

      if (strcmp(argv[0]+(strlen(argv[0])-5),".hist") == 0)
        DO_HIST = 1;
      else if (strcmp(argv[0]+(strlen(argv[0])-5),".ktab") == 0)
        DO_TABLE = 1;
      else if (strcmp(argv[0]+(strlen(argv[0])-5),".prof") == 0)
        DO_PROF = 1;
      else
        { f = fopen(Catenate(argv[0],".hist","",""),"r");
          if (f != NULL)
            { DO_HIST = 1;
              fclose(f);
            } 

          f = fopen(Catenate(argv[0],".ktab","",""),"r");
          if (f != NULL)
            { DO_TABLE = 1;
              fclose(f);
            } 

          f = fopen(Catenate(argv[0],".prof","",""),"r");
          if (f != NULL)
            { DO_PROF = 1;
              fclose(f);
            } 

          if (DO_HIST + DO_TABLE + DO_PROF == 0)
            { fprintf(stderr,"%s: There are no FastK objects with root %s !?\n",Prog_Name,argv[0]);
              exit (1);
            }
        }
    }
  }   

  if (DO_HIST)
    { int c, k;

      Histogram *H = Load_Histogram(argv[0]);

      for (c = 1; c < narg; c++) 
        { Histogram *G = Load_Histogram(argv[c]);

          if (G == NULL)
            { fprintf(stderr,"%s: Cannot open histogram with root %s\n",Prog_Name,argv[c]);
              exit (1);
            }

          if (G->low != H->low || G->high != G->high || G->kmer != H->kmer)
            { fprintf(stderr,"%s: K-mer histograms do not involve the same K or range\n",Prog_Name);
              exit (1);
            }

          for (k = G->low; k < G->high+3; k++)
            H->hist[k] += G->hist[k];

          Free_Histogram(G);
        }

      Write_Histogram(Catenate(Opath,"/",Oroot,".hist"), H);

      Free_Histogram(H);
    }

  if (DO_TABLE)
    {
      { int c;
    
        S = Malloc(sizeof(Kmer_Stream *)*narg,"Allocating table pointers");
        if (S == NULL)
          exit (1);
    
        kmer = 0;
        for (c = 0; c < narg; c++)
          { Kmer_Stream *s = Open_Kmer_Stream(argv[c]);
            if (s == NULL)
              { fprintf(stderr,"%s: Cannot open table %s\n",Prog_Name,argv[c]);
                exit (1);
              }
            if (c == 0)
              kmer = s->kmer;
            else
              { if (s->kmer != kmer)
                  { fprintf(stderr,"%s: K-mer tables do not involve the same K\n",Prog_Name);
                    exit (1);
                  }
              }
            S[c] = s;
          }
      }
    
      { int64     range[NTHREADS+1][narg];
#ifndef DEBUG_THREADS
        pthread_t threads[NTHREADS];
#endif
        TP        parm[NTHREADS];
        int64    *prefx;
        int       ixlen = 0;
        char     *seq;
        uint8    *ent;
        int       t, a, i;
        int64     p;
    
        ixlen = 0x1000000;
        prefx = Malloc(sizeof(int64)*ixlen,"Allocating prefix table");
        bzero(prefx,sizeof(int64)*ixlen);
    
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
            parm[t].out   = fopen(Catenate(Opath,"/.",Oroot,
                                  Numbered_Suffix(".ktab.",t+1,"")),"w");
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
    
        { int minval;
          int three = 3;
    
          FILE  *f   = fopen(Catenate(Opath,"/",Oroot,".ktab"),"w");
          int64 *prf = prefx;

          minval = S[i]->minval;
          for (i = 1; i < narg; i++)
            if (S[i]->minval < minval)
              minval = S[i]->minval;
    
          fwrite(&kmer,sizeof(int),1,f);
          fwrite(&NTHREADS,sizeof(int),1,f);
          fwrite(&minval,sizeof(int),1,f);
          fwrite(&three,sizeof(int),1,f);
    
          for (i = 1; i < ixlen; i++)
            prf[i] += prf[i-1];
    
          fwrite(prf,sizeof(int64),ixlen,f);
          fclose(f);
    
          free(prefx);
        }
      }
    
      { int c;
    
        for (c = 0; c < narg; c++)
          Free_Kmer_Stream(S[c]);
        free(S);
      }
    }

  if (DO_PROF)
    { PP        parm[NTHREADS];
#ifndef DEBUG_THREADS
      pthread_t threads[NTHREADS];
#endif
      int64 psize[narg], totreads, rpt;

      { char *path, *root;
        int   imer, nthreads;
        int64 nreads;
        FILE *f;
        int   c, t;
    
        totreads = 0;
        for (c = 0; c < narg; c++)
          { path = PathTo(argv[c]);
            root = Root(argv[c],".prof");
            f = fopen(Catenate(path,"/",root,".prof"),"r");
            if (f == NULL)
              { fprintf(stderr,"%s: Cannot open profile %s\n",Prog_Name,argv[c]);
                exit (1);
              }
            fread(&imer,sizeof(int),1,f);
            fread(&nthreads,sizeof(int),1,f);
            fclose(f);

            if (c == 0)
              kmer = imer;
            else
              { if (imer != kmer)
                  { fprintf(stderr,"%s: K-mer profiles do not involve the same K\n",Prog_Name);
                    exit (1);
                  }
              }

            for (t = 1; t <= nthreads; t++)
              { f = fopen(Catenate(path,"/.",root,Numbered_Suffix(".pidx.",t,"")),"r");
                if (f == NULL)
                  { fprintf(stderr,"%s: Cannot open profile index for %s\n",Prog_Name,argv[c]);
                    exit (1);
                  }
                fread(&imer,sizeof(int),1,f);
                fread(&nreads,sizeof(int64),1,f);
                fread(&nreads,sizeof(int64),1,f);
                totreads += nreads;
                fclose(f);
              }

            psize[c] = totreads;
            free(root);
            free(path);
          }

        rpt = totreads/NTHREADS;
      }

      { int64 u;
        int c, t;

        c = 0;
        u = 0;
        for (t = 0; t < NTHREADS; t++)
          { parm[t].tid   = t;
            parm[t].argv  = argv;
            parm[t].kmer  = kmer;
            parm[t].fidx  = u;

            u = ((t+1)*totreads)/NTHREADS;
            while (psize[c] <= u)
              c += 1;

            parm[t].nidx  = u-parm[t].fidx;
            parm[t].lpart = c;
            parm[t].last  = u - psize[c];

            parm[t].pout = fopen(Catenate(Opath,"/.",Oroot,
                                  Numbered_Suffix(".pidx.",t+1,"")),"w");
            parm[t].dout = fopen(Catenate(Opath,"/.",Oroot,
                                 Numbered_Suffix(".prof.",t+1,"")),"w");
          }
        parm[0].fpart = 0;
        parm[0].first = 0;
        for (t = 1; t < NTHREADS; t++)
          { parm[t].fpart = parm[t-1].lpart;
            parm[t].first = parm[t-1].last;
          }

#ifdef DEBUG_THREADS
        printf("\n");
        for (t = 0; t < NTHREADS; t++)
          printf("%d/%lld - %d/%lld\n",parm[t].fpart,parm[t].first,parm[t].lpart,parm[t].last);
#endif
    
#ifdef DEBUG_THREADS
        for (t = 0; t < NTHREADS; t++)
          prof_thread(parm+t);
#else
        for (t = 1; t < NTHREADS; t++)
          pthread_create(threads+t,NULL,prof_thread,parm+t);
        prof_thread(parm);
        for (t = 1; t < NTHREADS; t++)
          pthread_join(threads[t],NULL);
#endif
      }
    }

  free(Oroot);
  free(Opath);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
