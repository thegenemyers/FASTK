/*********************************************************************************************\
 *
 *  Example code for opening and fetching compressed profiles produced by FastK
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
#include <math.h>

#undef SHOW_RUN

#include "gene_core.h"

static char *Usage = "<source_root>.K<k> <read:int> ...";


/****************************************************************************************
 *
 *  The interface you may want to lift for your own use
 *
 *****************************************************************************************/

typedef struct
  { int    kmer;     //  Kmer length
    int    nparts;   //  # of threads/parts for the profiles
    int    nreads;   //  # of threads/parts for the profiles
    int64 *nbase;    //  nbase[i] = id of last read in part i
    FILE **nfile;    //  nfile[i] = stream for "P" file of part i
    int64 *index;    //  Kmer+count entry in bytes
  } Profile_Index;

Profile_Index *Open_Profiles(char *name);

void Free_Profiles(Profile_Index *P);

int Fetch_Profile(Profile_Index *P, int64 id, int plen, uint16 *profile);


/****************************************************************************************
 *
 *  Open a profile as a Profile_Index.  Index to compressed profiles is in memory,
 *    but compressed profiles are left on disk and reaad only when requested.
 *
 *****************************************************************************************/

Profile_Index *Open_Profiles(char *name)
{ Profile_Index *P;
  int            kmer, nparts;
  int64          nreads, *nbase, *index;
  FILE         **nfile;

  FILE  *f;
  int64  n;

  //  Find all parts and accumulate total size

  nparts = 0;
  nreads = 0;
  for (nparts = 0; 1; nparts++)
    { f = fopen(Catenate(name,Numbered_Suffix(".A",nparts+1,""),"",""),"r");
      if (f == NULL)
        break;
      fread(&kmer,sizeof(int),1,f);
      fread(&n,sizeof(int64),1,f);
      fread(&n,sizeof(int64),1,f);
      nreads += n;
      fclose(f);
    }
  if (nparts == 0)
    { fprintf(stderr,"%s: Cannot find profile files for %s\n",Prog_Name,name);
      exit (1);
    }

  //  Allocate in-memory table

  index = Malloc((nreads+1)*sizeof(int64),"Allocating profile index");
  nbase = Malloc(nparts*sizeof(int64),"Allocating profile index");
  nfile = Malloc(nparts*sizeof(FILE *),"Allocating profile index");
  if (index == NULL || nbase == NULL || nfile == NULL)
    exit (1);

  nparts = 0;
  nreads = 0;
  index[0] = 0;
  for (nparts = 0; 1; nparts++)
    { f = fopen(Catenate(name,Numbered_Suffix(".A",nparts+1,""),"",""),"r");
      if (f == NULL)
        break;
      fread(&kmer,sizeof(int),1,f);
      fread(&n,sizeof(int64),1,f);
      fread(&n,sizeof(int64),1,f);
      fread(index+(nreads+1),sizeof(int64),n,f);
      nreads += n;
      nbase[nparts] = nreads;
      fclose(f);

      f = fopen(Catenate(name,Numbered_Suffix(".P",nparts+1,""),"",""),"r");
      if (f == NULL)
        { fprintf(stderr,"%s: Cannot find profile file %s.P%d\n",Prog_Name,name,nparts+1);
          exit (1);
        }
      nfile[nparts] = f;
    }

  P = Malloc(sizeof(Profile_Index),"Allocating profile record");
  if (P == NULL)
    exit (1);

  P->kmer   = kmer;
  P->nparts = nparts;
  P->nreads = nreads;
  P->index  = index;
  P->nbase  = nbase;
  P->nfile  = nfile;

  return (P);
}


/****************************************************************************************
 *
 *  Free a Profile_Index and fetch a profile
 *
 *****************************************************************************************/

void Free_Profiles(Profile_Index *P)
{ int i;

  free(P->index);
  free(P->nbase);
  for (i = 0; i < P->nparts; i++)
    fclose(P->nfile[i]);
  free(P->nfile);
  free(P);
}

  //  Places uncompressed profile for read id (0-based) in profile of length plen.
  //    Returns the length of the uncompressed profile.  If the plen is less than
  //    this then only the first plen counts are uncompressed into profile

int Fetch_Profile(Profile_Index *P, int64 id, int plen, uint16 *profile)
{ uint8 count[1000], *cend = count+999;
  FILE  *f;
  int    w, len;
  uint8 *p, *q;
  uint16 x, d, i;
  int    n;

  for (w = 0; w < P->nparts; w++)
    if (id < P->nbase[w])
      break;
  if (w >= P->nparts)
    { fprintf(stderr,"%s: Id %lld is out of range [1,%lld]\n",Prog_Name,id,P->nbase[P->nparts-1]);
      exit (1);
    }
  f = P->nfile[w];

  if (id == 0 || (w > 0 && id == P->nbase[w-1]))
    { fseek(f,0,SEEK_SET);
      len = P->index[id+1];
    }
  else
    { int64 off = P->index[id];
      fseek(f,off,SEEK_SET);
      len = P->index[id+1] - off;
    }

  fread(count,1000,1,f);

  p = count;
  q = count + len;

  x = *p++;
  if ((x & 0x80) != 0)
    d = ((x & 0x7f) << 8) | *p++;
  else
    d = x;
  n = 1;

  if (plen > 0)
    { profile[0] = d;
#ifdef SHOW_RUN
  printf(" %d\n",d);
#endif

      while (p < q)
        { if (p >= cend)
            { if (p == cend)
                { *count = *p; 
                  fread(count+1,999,1,f);
                  p  = count;
                  q -= 999;
                }
              else
                { fread(count,1000,1,f);
                  p = count;
                  q -= 1000;
                }
            }
          x = *p++;
          if ((x & 0xc0) == 0)
            { if (n+x > plen)
                { n += x;
                  break;
                }
              for (i = 0; i < x; i++)
                profile[n++] = d;
#ifdef SHOW_RUN
              printf(" [%hu]\n",x);
#endif
            }
          else
            { if ((x & 0x80) != 0)
                { if ((x & 0x40) != 0)
                    x <<= 8;
                  else
                    x = (x << 8) & 0x7fff;
                  x |= *p++;
                  d += x;
#ifdef SHOW_RUN
                  printf(" %hd+(%d)\n",x,d);
#endif
                }
              else
                { if ((x & 0x20) != 0)
                    d += (x & 0x1fu) | 0xffe0u;
                  else
                    d += (x & 0x1fu);
#ifdef SHOW_RUN
                  if ((x & 0x20) != 0)
                    printf(" -%d(%d)\n",32-(x&0x1fu),d);
                  else
                    printf(" +%d(%d)\n",x&0x1fu,d);
#endif
                }
              if (n >= plen)
                { n += 1;
                  break;
                }
              profile[n++] = d;
            }
        }
    }

  while (p < q)
    { if (p >= cend)
        { if (p == cend)
            { *count = *p; 
              fread(count+1,999,1,f);
              p  = count;
              q -= 999;
            }
          else
            { fread(count,1000,1,f);
              p = count;
              q -= 1000;
            }
        } 
      x = *p++;
      if ((x & 0xc0) == 0)
        n += x;
      else
        { if ((x & 0x80) != 0)
            p += 1;
          n += 1;
        }
    }

  return (n);
}


/****************************************************************************************
 *
 *  Test Stub
 *
 *****************************************************************************************/

int main(int argc, char *argv[])
{ Profile_Index *P;

  { int    i, j, k;
    int    flags[128];

    (void) flags;

    ARG_INIT("Profex");

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc < 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  P = Open_Profiles(argv[1]);

  { int     c, id;
    char   *eptr;
    uint16 *profile;
    int     plen, tlen;

    plen    = 20000;
    profile = Malloc(plen*sizeof(uint16),"Profile array");

    for (c = 2; c < argc; c++)
      { id = strtol(argv[c],&eptr,10);
        if (*eptr != '\0' || argv[c][0] == '\0')
          { fprintf(stderr,"%s: argument '%s' is not an integer\n",Prog_Name,argv[c]);
             exit (1);
          }
        if (id <= 0 || id > P->nbase[P->nparts-1])
          { fprintf(stderr,"%s: Id %d is out of range\n",Prog_Name,id);
            exit (1);
          }
        tlen = Fetch_Profile(P,(int64) id-1,plen,profile);
        if (tlen > plen)
          { plen    = 1.2*tlen + 1000;
            profile = Realloc(profile,plen*sizeof(uint16),"Profile array");
            Fetch_Profile(P,(int64) id-1,plen,profile);
          }
        printf("\nRead %d:\n",id);
        for (int i = 0; i < tlen; i++)
          printf(" %5d: %5d\n",i,profile[i]);
      }
    free(profile);
  }

  Free_Profiles(P);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
