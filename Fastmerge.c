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

#undef    DEBUG
#undef    DEBUG_THREADS
#undef    DEBUG_TRACE
#undef    DEBUG_PROF

#include "libfastk.h"

static char *Usage = " [-htp] [-T<int(4)>] <target> <sources>[.hist|.ktab|.prof] ...";

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
    int64        *hist;
    int           dotab;
  } TP;

static inline int mycmp(uint8 *a, uint8 *b, int n)
{ while (n-- > 0)
    { if (*a++ != *b++)
        return (a[-1] < b[-1] ? -1 : 1);
    }
  return (0);
}

static void *table_thread(void *args)
{ TP *parm = (TP *) args;
  int           tid   = parm->tid;
  int           ntabs = parm->narg;
  Kmer_Stream **S     = parm->S;
  Kmer_Stream **T;
  int64        *begs  = parm->begs;
  int64        *ends  = parm->ends;
  int64        *prefx = parm->prefx;
  FILE         *out   = parm->out;
  int           dotab = parm->dotab;

  int hbyte = S[0]->kbyte-3;
  int kbyte = S[0]->kbyte;
  int kmer  = S[0]->kmer;

  int64  *hist;
  int64   nels;
  uint8 **ent, *bst;
  int     itop, *in, cnt;
  uint16  scnt;
  int     c, x;

#ifdef DEBUG_TRACE
  char *buffer;
#endif

  hist = Malloc(sizeof(int64)*0x8001,"Allocating histogram");
  bzero(hist,sizeof(int64)*0x8001);
  hist -= 1;

  in  = Malloc(sizeof(int)*ntabs,"Allocating thread working memory");
  ent = Malloc(sizeof(uint8 *)*ntabs,"Allocating thread working memory");

#ifdef DEBUG_THREADS
  printf("Doing %d:\n",tid);
  for (c = 0; c < ntabs; c++)
    printf("  %2d: [%lld-%lld]",c,begs[c],ends[c]);
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
   if (dotab)
     { fwrite(&kmer,sizeof(int),1,out);
       fwrite(&nels,sizeof(int64),1,out);
     }

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
      if (cnt >= 0x7fff)
        { scnt = 0x7fff;
          hist[0x7fff] += 1;
          for (c = 0; c < itop; c++)
            { cnt = Current_Count(T[in[c]]);
              if (cnt < 0x7fff)
                hist[0x8001] += cnt;
            }
        }
      else
        { scnt = cnt;
          hist[cnt] += 1;
        }

#ifdef DEBUG_TRACE
      for (c = 0; c < itop; c++)
        { x = in[c];
          printf(" %d: %s %5d",x,Current_Kmer(T[x],buffer),Current_Count(T[x]));
        }
      printf("\n");
#endif

      if (dotab)
        { fwrite(bst+3,hbyte,1,out);
          fwrite(&scnt,sizeof(uint16),1,out);
          x = (bst[0] << 16) | (bst[1] << 8) | bst[2];
          prefx[x] += 1;
          nels += 1;
        }

      for (c = 0; c < itop; c++)
        { Kmer_Stream *t;

          x = in[c];
          t = T[x];
          Next_Kmer_Entry(t);
          if (t->csuf != NULL)
            Current_Entry(t,ent[x]);
        }
    }

   if (dotab)
     { rewind(out);
       fwrite(&kmer,sizeof(int),1,out);
       fwrite(&nels,sizeof(int64),1,out);
       fclose(out);
     }

  for (c = 0; c < ntabs; c++)
    free(ent[c]);

  if (tid != 0)
    { for (c = 0; c < ntabs; c++) 
        Free_Kmer_Stream(T[c]);
      free(T);
    }

  free(ent);
  free(in);

  parm->hist = hist;
  return (NULL);
}


/****************************************************************************************
 *
 *  Streaming threaded merge of profile indices and data
 *
 *****************************************************************************************/

typedef struct
  { char        **argv;   //  input profiles
#ifdef DEBUG_PROF
    int           tid;
#endif
    int           kmer;   
    int64         fidx;   //  output will contain profiles [fidx,fidx+nidx) 
    int64         nidx;
    FILE         *pout;   //  output files for profile .pidx and .prof part
    FILE         *dout;
    int           fpart;  //  output is fpart:first to lpart:last in terms of input blocks
    int64         first;
    int           lpart;
    int64         last;
    int64         buffer[16384];  //  buffer for data transfer
  } PP;

static void *prof_thread(void *args)
{ PP    *parm   = (PP *) args;
  char **argv   = parm->argv;
  FILE  *dout   = parm->dout;
  FILE  *pout   = parm->pout;
  int    fpart  = parm->fpart;
  int64  first  = parm->first;
  int    lpart  = parm->lpart;
  int64  last   = parm->last;
  int64 *buffer = parm->buffer;

  int64  q, cindex, lindex;
  int64  base, offs;
  int64  fdata, ldata;
  char  *path, *root, *name;
  int    imer, nthreads;
  FILE  *f, *g;
  int    c, t;


  fwrite(&(parm->kmer),sizeof(int),1,pout);
  fwrite(&(parm->fidx),sizeof(int64),1,pout);
  fwrite(&(parm->nidx),sizeof(int64),1,pout);

  base = 0;
  for (c = fpart; c <= lpart; c++) 
    { path = PathTo(argv[c]);
      root = Root(argv[c],".prof");

      if (c == lpart && last == 0)
        break;

      name = Malloc(strlen(path) + strlen(root) + 100,"Allocating path name");

      sprintf(name,"%s/%s.prof",path,root);
      f = fopen(name,"r");
      if (f == NULL)
        { fprintf(stderr,"%s: Cannot open profile %s %d %lld\n",Prog_Name,argv[c],c,last);
          exit (1);
        }
      fread(&imer,sizeof(int),1,f);
      fread(&nthreads,sizeof(int),1,f);
      fclose(f);
#ifdef DEBUG_PROF
      printf("Opening %d: %d: n=%d\n",parm->tid,c,nthreads);
#endif

      for (t = 1; t <= nthreads; t++)
        { sprintf(name,"%s/.%s.pidx.%d",path,root,t);
          f = fopen(name,"r");
          if (f == NULL)
            { fprintf(stderr,"%s: Cannot open profile index for %s\n",Prog_Name,argv[c]);
              exit (1);
            }
          sprintf(name,"%s/.%s.prof.%d",path,root,t);
          g = fopen(name,"r");
          if (g == NULL)
            { fprintf(stderr,"%s: Cannot open profile data for %s\n",Prog_Name,argv[c]);
              exit (1);
            }

          fread(&imer,sizeof(int),1,f);
          fread(&cindex,sizeof(int64),1,f);
          fread(&lindex,sizeof(int64),1,f);
          lindex += cindex;
#ifdef DEBUG_PROF
          printf("Pidx %d: %d: %d: %lld - %lld base = %lld\n",parm->tid,c,t,cindex,lindex,base);
#endif
          fdata = 0;
          if (c == fpart && cindex <= first)
            { if (lindex <= first) 
                { if (lindex == first)
                    { fseeko(f,sizeof(int64)*((first-cindex)-1),SEEK_CUR);
                      fread(&fdata,sizeof(int64),1,f);
                    }
                  fclose(f);
                  continue;
                }
              if (cindex+1 < first)
                fseeko(f,sizeof(int64)*((first-cindex)-1),SEEK_CUR);
              if (cindex < first)
                { fread(&offs,sizeof(int64),1,f);
                  fdata = offs;
                }
              base = -fdata;
              cindex = first;
            }
          if (c == lpart && lindex >= last)
            { lindex = last;
              t = nthreads;
            }
#ifdef DEBUG_PROF
          printf("Pidx %d: %d: %d: %lld - %lld base = %lld\n",parm->tid,c,t,cindex,lindex,base);
#endif
          for (q = cindex; q < lindex; q++)
            { fread(&ldata,sizeof(int64),1,f);
              offs = ldata + base;
              fwrite(&offs,sizeof(int64),1,pout);
            }
          base += ldata;
          fclose(f);

#ifdef DEBUG_PROF
          printf("Prof %d: %d: %d: range = %10lld-%10lld\n",parm->tid,c,t,fdata,ldata);
#endif
          if (fdata > 0)
            fseeko(g,fdata,SEEK_SET);
          for (q = ldata-fdata; q >= 16384; q -= 16384)
            { fread(buffer,16384,1,g);
              fwrite(buffer,16384,1,dout);
            }
          if (q > 0)
            { fread(buffer,q,1,g);
              fwrite(buffer,q,1,dout);
            }
          fclose(g);
        }

      free(name);
      free(root);
      free(path);
    }
  fwrite(&offs,sizeof(int64),1,pout);

  fclose(pout);
  fclose(dout);

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

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("htp")
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
        }   
      else
        argv[j++] = argv[i];
    argc = j;

    DO_HIST  = flags['h'];
    DO_TABLE = flags['t'];
    DO_PROF  = flags['p'];

    if (argc < 4)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -h: Produce a merged histogram.\n");
        fprintf(stderr,"      -t: Produce a merged k-mer table.\n");
        fprintf(stderr,"      -p: Produce a merged profile.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -T: Use -T threads.\n");
        exit (1);
      } 

    for (i = 1; i < argc; i++)
      { int dot = strlen(argv[i])-5;
        if (strcmp(argv[i]+dot,".hist") == 0)
          { argv[i][dot] = '\0';
            continue;
          }
        if (strcmp(argv[i]+dot,".ktab") == 0)
          { argv[i][dot] = '\0';
            continue;
          }
        if (strcmp(argv[i]+dot,".prof") == 0)
          { argv[i][dot] = '\0';
            continue;
          }
      }

    Opath = PathTo(argv[1]);
    Oroot = Root(argv[1],".ktab");

    { FILE *f;
      int   has_hist, has_table, has_prof;
    
      narg = argc-2;
      argv += 2;

      has_hist = has_table = has_prof = 0;
      for (i = 0; i < narg; i++)
        { f = fopen(Catenate(argv[0],".hist","",""),"r");
          if (f != NULL)
            { has_hist += 1;
              fclose(f);
            } 
          f = fopen(Catenate(argv[0],".ktab","",""),"r");
          if (f != NULL)
            { has_table += 1;
              fclose(f);
            } 
          f = fopen(Catenate(argv[0],".prof","",""),"r");
          if (f != NULL)
            { has_prof += 1;
              fclose(f);
            } 
        }

      if (has_hist != 0 && has_hist != narg)
        { fprintf(stderr,"%s: Some sources have .hist files, others do not?\n",Prog_Name);
          exit (1);
	}
      if (has_table != 0 && has_table != narg)
        { fprintf(stderr,"%s: Some sources have .ktab files, others do not?\n",Prog_Name);
          exit (1);
	}
      if (has_prof != 0 && has_prof != narg)
        { fprintf(stderr,"%s: Some sources have .prof files, others do not?\n",Prog_Name);
          exit (1);
	}

      if (DO_HIST + DO_TABLE + DO_PROF == 0)
        { DO_HIST  = (has_hist > 0);
          DO_TABLE = (has_table > 0);
          DO_PROF  = (has_prof > 0);
        }
      else
        { if (DO_HIST && has_hist == 0)
            { fprintf(stderr,"%s: Requesting histogram merge but no source .hist's\n",Prog_Name);
              exit (1);
            }
          if (DO_TABLE && has_table == 0)
            { fprintf(stderr,"%s: Requesting table merge but no source .ktab's\n",Prog_Name);
              exit (1);
            }
          if (DO_PROF && has_prof == 0)
            { fprintf(stderr,"%s: Requesting profile merge but no source .ktab's\n",Prog_Name);
              exit (1);
            }
        }

      if (DO_HIST && has_table == 0)
        { fprintf(stderr,"%s: Need source tables to compute merged histograms\n",Prog_Name);
          exit (1);
        }
    }
  }   

  if (DO_HIST || DO_TABLE)
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
        int       pivot;
        int       t, a, i;
        int64     p;

        if (DO_TABLE)
          { ixlen = 0x1000000;
            prefx = Malloc(sizeof(int64)*ixlen,"Allocating prefix table");
            bzero(prefx,sizeof(int64)*ixlen);
          }
        else
          prefx = NULL;
    
        for (a = 0; a < narg; a++)
          { range[0][a] = 0;
            range[NTHREADS][a] = S[a]->nels;
          }

        pivot = 0;
        for (a = 1; a < narg; a++)
          if (S[a]->nels > S[pivot]->nels)
            pivot = a;

        seq = Current_Kmer(S[0],NULL);
        ent = Current_Entry(S[0],NULL);
        for (t = 1; t < NTHREADS; t++)
          { p = (S[pivot]->nels*t)/NTHREADS; 
            GoTo_Kmer_Index(S[pivot],p);
#ifdef DEBUG
            printf("\n%d: %0*x\n",t,2*S[pivot]->ibyte,S[pivot]->cpre);
            printf(" %lld:",p);
            if (p < S[pivot]->nels)
              printf(" %s\n",Current_Kmer(S[pivot],seq));
            else
              printf(" EOT\n");
#endif
            if (p >= S[pivot]->nels)
              for (a = 0; a < narg; a++)
                range[t][a] = S[a]->nels;
            else
              { ent = Current_Entry(S[pivot],ent);                //  Break at prefix boundaries
                for (i = S[0]->ibyte; i < S[0]->kbyte; i++)
                  ent[i] = 0;
                for (a = 0; a < narg; a++)
                  { GoTo_Kmer_Entry(S[a],ent);
#ifdef DEBUG
                    printf("  %lld:",S[a]->cidx);
                    if (S[a]->cidx < S[a]->nels)
                      printf(" %s\n",Current_Kmer(S[a],seq));
                    else
                      printf(" EOT\n");
#endif
                    range[t][a] = S[a]->cidx;
                  }
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
            parm[t].dotab = DO_TABLE;
            if (DO_TABLE)
              parm[t].out = fopen(Catenate(Opath,"/.",Oroot,
                                   Numbered_Suffix(".ktab.",t+1,"")),"w");
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

        if (DO_HIST)
          { int64 *hist, *gist;
            int    j, low, high;
            int    f;

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
                  { fprintf(stderr,"%s: Cannot open histogram %s\n",Prog_Name,argv[t]);
                    exit (1);
                  }
                hist[0x8001] += H->hist[0x8001];
                Free_Histogram(H);
              }
            hist[0x8000] = hist[1];

            low  = 1;
            high = 0x7fff;
            f  = open(Catenate(Opath,"/",Oroot,".hist"),O_CREAT|O_TRUNC|O_WRONLY,S_IRWXU);
            write(f,&kmer,sizeof(int));
            write(f,&low,sizeof(int));
            write(f,&high,sizeof(int));
            write(f,hist+0x8000,sizeof(int64));
            write(f,hist+0x8001,sizeof(int64));
            write(f,hist+1,sizeof(int64)*0x7fff);
            close(f);

            free(hist+1);
          }

        if (DO_TABLE)
          { int minval;
            int three = 3;
    
            FILE  *f   = fopen(Catenate(Opath,"/",Oroot,".ktab"),"w");
            int64 *prf = prefx;

            minval = S[0]->minval;
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
      int64 psize[narg+1], totreads;

      //  Determine total # of profiles (totreads), and # in each arg (psize[c])

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
#ifdef DEBUG
                if (c == 0)
                  printf(" %d: %d: %10lld %10lld %10lld\n",
                         c+1,t,nreads,totreads,totreads);
                else
                  printf(" %d: %d: %10lld %10lld %10lld\n",
                         c+1,t,nreads,totreads-psize[c-1],totreads);
#endif
                fclose(f);
              }

            psize[c] = totreads;
#ifdef DEBUG
            printf("   psize[%d] = %10lld\n",c,psize[c]);
#endif
            free(root);
            free(path);
          }
        psize[c] = totreads+1;
      }

      //  Setup division into threads in 'parm'

      { int64 u;
        int c, t;

        c = 0;
        u = 0;
        for (t = 0; t < NTHREADS; t++)
          { parm[t].argv  = argv;
#ifdef DEBUG_PROF
            parm[t].tid   = t;
#endif
            parm[t].kmer  = kmer;
            parm[t].fidx  = u;

            u = ((t+1)*totreads)/NTHREADS;
            while (psize[c] <= u)
              c += 1;

            parm[t].nidx  = u-parm[t].fidx;
            parm[t].lpart = c;
            if (c > 0)
              parm[t].last  = u - psize[c-1];
            else
              parm[t].last  = u;

            parm[t].pout = fopen(Catenate(Opath,"/.",Oroot,
                                  Numbered_Suffix(".pidx.",t+1,"")),"w");
            if (parm[t].pout == NULL)
              { fprintf(stderr,"%s: Cannot create part .pidx file for ouput %s\n",Prog_Name,Oroot);
                exit (1);
              }
            parm[t].dout = fopen(Catenate(Opath,"/.",Oroot,
                                 Numbered_Suffix(".prof.",t+1,"")),"w");
            if (parm[t].dout == NULL)
              { fprintf(stderr,"%s: Cannot create part .prof file for ouput %s\n",Prog_Name,Oroot);
                exit (1);
              }
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
        printf("\n");
        for (t = 0; t < NTHREADS; t++)
          printf("  %lld + %lld = %10lld\n",parm[t].fidx,parm[t].nidx,parm[t].fidx+parm[t].nidx);
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

      //  Create stub file

      { FILE *f;

        f = fopen(Catenate(Opath,"/",Oroot,".prof"),"w");
        if (f == NULL)
          { fprintf(stderr,"%s: Cannot create stub file for ouput %s\n",Prog_Name,Oroot);
            exit (1);
          }
        fwrite(&kmer,sizeof(int),1,f);
        fwrite(&NTHREADS,sizeof(int),1,f);
        fclose(f);
      }
    }

  free(Oroot);
  free(Opath);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
