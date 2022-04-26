/*******************************************************************************************
 *
 *  Merges tables, histograms, and profiles produced by independent HPC runs on sub-parts
 *     of a data set
 *
 *  Author:  Gene Myers
 *  Date  :  Sep. 20, 2021
 *  Date  :  Mar. 31, 2022  priority queue, tailored prefix size, thread parts
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

static char *Usage = " [-htp] [-T<int(4)>] [-P<int(1)>} <target> <sources>[.hist|.ktab|.prof] ...";

static int   NTHREADS;
static int   NPARTS;
static char *OPATH;
static char *OROOT;
static int   PIVOT;

/****************************************************************************************
 *
 *  Streaming threaded merge of k-mer tables
 *
 *****************************************************************************************/

typedef struct
  { int           tid;
    Kmer_Stream **S;
    int           narg;
    int64        *prefx;
    int           ibyte;
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

  //  Heap of input buffer pointers ordering on k-mer at ->ptr

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


static int next_group(int node, int *equb, int gtop, int *group)
{ if (equb[node] & EQ_RGHT)
    gtop = next_group(2*node+1,equb,gtop,group);
  if (equb[node] & EQ_LEFT)
    gtop = next_group(2*node,equb,gtop,group);
  group[++gtop] = node;
  return (gtop);
}

static void *table_thread(void *args)
{ TP *parm = (TP *) args;
  int           tid   = parm->tid;
  int           ntabs = parm->narg;
  Kmer_Stream **S     = parm->S;
  int64        *begs  = parm->begs;
  int64        *ends  = parm->ends;
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
  printf("Doing %d: pivot %d\n",tid,PIVOT);
  for (c = 0; c < ntabs; c++)
    printf("  %2d: [%lld-%lld]\n",c,begs[c],ends[c]);
  printf("\n");
#endif

  if (tid != 0)
    { T = Malloc(sizeof(Kmer_Stream *)*ntabs,"Allocating thread working memory");
      for (c = 0; c < ntabs; c++)
        T[c] = Clone_Kmer_Stream(S[c]);
    }
  else
    T = S;
  p = T[PIVOT];

#ifdef DEBUG_TRACE
  buffer = Current_Kmer(T[0],NULL);
#endif

   if (dotab)
     { npr  = 1;
       pend = begs[PIVOT] + ((ends[PIVOT] - begs[PIVOT]) * npr) / NPARTS;

       sprintf(nbuf,"%s/.%s.ktab.%d",OPATH,OROOT,tid*NPARTS+npr);
       out = fopen(nbuf,"w");

#ifdef DEBUG_THREADS
       printf("  Making %s to %d:%lld\n",nbuf,npr,pend);
#endif

       nels = 0;
       fwrite(&kmer,sizeof(int),1,out);
       fwrite(&nels,sizeof(int64),1,out);
     }

  for (c = 0; c < ntabs; c++)
    GoTo_Kmer_Index(T[c],begs[c]);

  hsize = ntabs;
  for (c = 0; c < ntabs; c++)
    { ent[c] = Current_Entry(T[c],NULL);
      if (T[c]->cidx >= ends[c])
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

  while (hsize > 0)
    { gtop = next_group(1,equb,-1,grp);
      best = ent[heap[1]];

      icnt = (cnt[gtop] = *((uint16 *) (best + kbyte)));
      for (x = 0; x < gtop; x++)
        icnt += (cnt[x] = *((uint16 *) (ent[heap[grp[x]]] + kbyte)));

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
      printf("%lld: ",count);
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
  }
#endif

      if (dotab)
        { fwrite(best+ibyte,hbyte,1,out);
          fwrite(&scnt,sizeof(uint16),1,out);

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
              pend = begs[PIVOT] + ((ends[PIVOT] - begs[PIVOT]) * npr) / NPARTS;
              sprintf(nbuf,"%s/.%s.ktab.%d",OPATH,OROOT,tid*NPARTS+npr);
              out = fopen(nbuf,"w");

#ifdef DEBUG_THREADS
              printf("  Making %s to %d:%lld\n",nbuf,npr,pend);
#endif

              nels = 0;
              fwrite(&kmer,sizeof(int),1,out);
              fwrite(&nels,sizeof(int64),1,out);
            }
        }

      for (x = 0; x <= gtop; x++)
        { i = grp[x];
          c = heap[i];
          t = T[c];
          Next_Kmer_Entry(t);
          if (t->cidx < ends[c])
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

  free(nbuf);
  free(equb);
  free(heap);
  free(ent);
  free(grp);
  free(cnt);

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
    NPARTS   = 1;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("htp")
            break;
          case 'P':
            ARG_POSITIVE(NPARTS,"Number of parts per thread")
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
        fprintf(stderr,"      -P: Produce -P parts per thread.\n");
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

    OPATH = PathTo(argv[1]);
    OROOT = Root(argv[1],".ktab");

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
        int       ibyte = 3;
        char     *seq;
        uint8    *ent;
        int       t, a, i;
        int64     p;

        if (DO_TABLE)
          { int64 nels = 0;
            for (a = 1; a < narg; a++)
              nels = S[a]->nels;
            if (nels >= 0x8000000)
              { ixlen = 0x1000000;
                ibyte = 3;
              }
            else if (nels >= 0x80000)
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
        else
          prefx = NULL;
    
        for (a = 0; a < narg; a++)
          { range[0][a] = 0;
            range[NTHREADS][a] = S[a]->nels;
          }

        PIVOT = 0;
        for (a = 1; a < narg; a++)
          if (S[a]->nels > S[PIVOT]->nels)
            PIVOT = a;

        seq = Current_Kmer(S[0],NULL);
        ent = Current_Entry(S[0],NULL);
        for (t = 1; t < NTHREADS; t++)
          { p = (S[PIVOT]->nels*t)/NTHREADS; 
            GoTo_Kmer_Index(S[PIVOT],p);
#ifdef DEBUG
            printf("\n%d: %0*x\n",t,2*S[PIVOT]->ibyte,S[PIVOT]->cpre);
            printf(" %lld:",p);
            if (p < S[PIVOT]->nels)
              printf(" %s\n",Current_Kmer(S[PIVOT],seq));
            else
              printf(" EOT\n");
#endif
            if (p >= S[PIVOT]->nels)
              for (a = 0; a < narg; a++)
                range[t][a] = S[a]->nels;
            else
              { ent = Current_Entry(S[PIVOT],ent);                //  Break at prefix boundaries
                for (i = ibyte; i < S[0]->kbyte; i++)
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
            f  = open(Catenate(OPATH,"/",OROOT,".hist"),O_CREAT|O_TRUNC|O_WRONLY,S_IRWXU);
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
          { int minval, nparts;
    
            FILE  *f   = fopen(Catenate(OPATH,"/",OROOT,".ktab"),"w");
            int64 *prf = prefx;

            minval = S[0]->minval;
            for (i = 1; i < narg; i++)
              if (S[i]->minval < minval)
                minval = S[i]->minval;
            nparts = NTHREADS * NPARTS;
    
            fwrite(&kmer,sizeof(int),1,f);
            fwrite(&nparts,sizeof(int),1,f);
            fwrite(&minval,sizeof(int),1,f);
            fwrite(&ibyte,sizeof(int),1,f);
    
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

      { int64 u, plast;
        int c, t, p, plpart;

        c = 0;
        u = 0;
        plpart = 0;
        plast  = 0;

        for (p = 0; p < NPARTS; p++)
          { for (t = 0; t < NTHREADS; t++)
              { parm[t].argv  = argv;
#ifdef DEBUG_PROF
                parm[t].tid   = t;
#endif
                parm[t].kmer  = kmer;
                parm[t].fidx  = u;
    
                u = ((p*NTHREADS+t+1)*totreads)/(NTHREADS*NPARTS);
                while (psize[c] <= u)
                  c += 1;
    
                parm[t].nidx  = u-parm[t].fidx;
                parm[t].lpart = c;
                if (c > 0)
                  parm[t].last  = u - psize[c-1];
                else
                  parm[t].last  = u;
    
                parm[t].pout = fopen(Catenate(OPATH,"/.",OROOT,
                                      Numbered_Suffix(".pidx.",p*NTHREADS+t+1,"")),"w");
                if (parm[t].pout == NULL)
                  { fprintf(stderr,"%s: Cannot create part .pidx file for ouput %s\n",
                                   Prog_Name,OROOT);
                    exit (1);
                  }
                parm[t].dout = fopen(Catenate(OPATH,"/.",OROOT,
                                     Numbered_Suffix(".prof.",p*NTHREADS+t+1,"")),"w");
                if (parm[t].dout == NULL)
                  { fprintf(stderr,"%s: Cannot create part .prof file for ouput %s\n",
                                   Prog_Name,OROOT);
                    exit (1);
                  }
              }
            parm[0].fpart = plpart;
            parm[0].first = plast;
            for (t = 1; t < NTHREADS; t++)
              { parm[t].fpart = parm[t-1].lpart;
                parm[t].first = parm[t-1].last;
              }
            plpart = parm[NTHREADS-1].lpart;
            plast  = parm[NTHREADS-1].last;
    
#ifdef DEBUG_THREADS
            printf("\n");
            for (t = 0; t < NTHREADS; t++)
              printf("%d/%lld - %d/%lld\n",parm[t].fpart,parm[t].first,parm[t].lpart,parm[t].last);
            printf("\n");
            for (t = 0; t < NTHREADS; t++)
              printf("  %lld + %lld = %10lld\n",
                     parm[t].fidx,parm[t].nidx,parm[t].fidx+parm[t].nidx);
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

      //  Create stub file

      { FILE *f;
        int   nparts;

        f = fopen(Catenate(OPATH,"/",OROOT,".prof"),"w");
        if (f == NULL)
          { fprintf(stderr,"%s: Cannot create stub file for ouput %s\n",Prog_Name,OROOT);
            exit (1);
          }
        nparts = NTHREADS*NPARTS;
        fwrite(&kmer,sizeof(int),1,f);
        fwrite(&nparts,sizeof(int),1,f);
        fclose(f);
      }
    }

  free(OROOT);
  free(OPATH);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
