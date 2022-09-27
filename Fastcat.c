/*******************************************************************************************
 *
 *  Catenates the tables, histograms, and profiles produced by Fastmerge with the -S option
 *     of a data set
 *
 *  Author:  Gene Myers
 *  Date  :  June 1, 2022
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <pthread.h>

#include "libfastk.h"

static char *Usage = "[-vk] [-htp] <target> <source>[.hist|.ktab|.prof] ...";

/****************************************************************************************
 *
 *  Main
 *
 *****************************************************************************************/

int main(int argc, char *argv[])
{ int   narg;
  int   KEEP;
  int   VERBOSE;
  char *OPATH;
  char *OROOT;
  int   DO_HIST, DO_TABLE, DO_PROF;

  //  Process the command line

  { int    i, j, k;
    int    flags[128];

    ARG_INIT("Fastcat");

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vkhtp")
            break;
        }   
      else
        argv[j++] = argv[i];
    argc = j;

    KEEP     = flags['k'];
    VERBOSE  = flags['v'];
    DO_HIST  = flags['h'];
    DO_TABLE = flags['t'];
    DO_PROF  = flags['p'];

    if (argc < 4)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Print progress as you go.\n");
        fprintf(stderr,"      -k: Keep source files (requires copying all parts).\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -h: Produce a merged histogram.\n");
        fprintf(stderr,"      -t: Produce a merged k-mer table.\n");
        fprintf(stderr,"      -p: Produce a merged profile.\n");

        exit (1);
      } 

    if (DO_HIST + DO_TABLE + DO_PROF == 0)
      { fprintf(stderr,"%s: At least one of the flags -h, -t, or -p must be set\n",Prog_Name);
        exit (1);
      }

    //  Remove FastK extensions from source arguments if any
    //  The target should not have one

    for (i = 1; i < argc; i++)
      { int dot = strlen(argv[i])-5;
        if (strcmp(argv[i]+dot,".hist") == 0)
          argv[i][dot] = '\0';
        if (strcmp(argv[i]+dot,".ktab") == 0)
          argv[i][dot] = '\0';
        if (strcmp(argv[i]+dot,".prof") == 0)
          argv[i][dot] = '\0';
        if (i == 1 && argv[1][dot] == '\0')
          { fprintf(stderr,"\n%s: Target name cannot have a .hist, .ktab, or .prof suffix\n",
                           Prog_Name);
            exit (1);
          }
      }

    //  The destination path and root

    OPATH = PathTo(argv[1]);
    OROOT = Root(argv[1],"");

    //  Make sure that sources have a full complement of each FastK file type and at least one

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

      if (DO_HIST && has_hist != narg)
        { if (has_hist == 0)
            fprintf(stderr,"\n%s: None of the sources have .hist files?\n",Prog_Name);
          else
            fprintf(stderr,"\n%s: Some (%d) of the sources do not have .hist files?\n",
                           Prog_Name,narg-has_hist);
          exit (1);
	}
      if (DO_TABLE && has_table != narg)
        { if (has_table == 0)
            fprintf(stderr,"\n%s: None of the sources have .ktab files?\n",Prog_Name);
          else
            fprintf(stderr,"\n%s: Some (%d) of the sources do not have .ktab files?\n",
                           Prog_Name,narg-has_table);
          exit (1);
	}
      if (DO_PROF && has_prof != narg)
        { if (has_prof == 0)
            fprintf(stderr,"\n%s: None of the sources have .prof files?\n",Prog_Name);
          else
            fprintf(stderr,"\n%s: Some (%d) of the sources do not have .prof files?\n",
                           Prog_Name,narg-has_prof);
          exit (1);
	}
    }
  }   

  if (DO_TABLE)
    { char *path, *root, *command;
      int64 ixlen, idx, ldx, *prefx;
      int   smer, kmer;
      int   sbyte, ibyte;
      int   spart, npart;
      int   scnt, mcnt;
      int   nmax, error;
      FILE *f;
      int   i, j, c;

      // In a first pass read each source table stub:
      //   check that the kmer length and entry prefix lengths are the same for all sources
      //   compute the minimum kmer count threshold
      //   build a combined prefix table
      //   make sure each table part exists and you can open it
      //   determine the length of the longest source name (path+root)

      if (VERBOSE)
        { printf("\n  Cat'ing %d tables\n",narg);
          fflush(stdout);
        }

      error = 0;
      kmer  = 0;
      npart = 0;
      mcnt  = 0x10000;
      nmax  = 0;
      for (c = 0; c < narg; c++)
        { path = PathTo(argv[c]);
          root = Root(argv[c],"");
          f = fopen(Catenate(path,"/",root,".ktab"),"r");
          if (f == NULL)
            { fprintf(stderr,"\n%s: Cannot find and open table %s/%s.ktab\n",Prog_Name,path,root);
              exit (1);
            }

          if ((int) (strlen(path)+strlen(root)) > nmax)
            nmax = strlen(path)+strlen(root);

          error |= (fread(&smer,sizeof(int),1,f) < 1);
          error |= (fread(&spart,sizeof(int),1,f) < 1);
          error |= (fread(&scnt,sizeof(int),1,f) < 1);
          error |= (fread(&sbyte,sizeof(int),1,f) < 1);
          if (error)
            { fprintf(stderr,"\n%s: Cannot read table %s/%s.ktab\n",Prog_Name,path,root);
              exit (1);
            }

          if (VERBOSE)
            { printf("  Processing table header of %s\n",root);
              fflush(stdout);
            }

          if (c == 0)
            { kmer  = smer;
              mcnt  = scnt;
              ibyte = sbyte;
              ixlen = (0x1 << (8*sbyte));
              prefx = Malloc(sizeof(int64)*ixlen,"Allocating prefix table");
              if (prefx == NULL)
                exit (1);
              bzero(prefx,sizeof(int64)*ixlen);
            }
          else
            { if (smer != kmer)
                { fprintf(stderr,"\n%s: K-mer tables do not involve the same K\n",Prog_Name);
                  exit (1);
                }
              if (sbyte != ibyte)
                { fprintf(stderr,"\n%s: K-mer tables do not involve the same K\n",Prog_Name);
                  exit (1);
                }
              if (scnt < mcnt)
                mcnt = scnt;
            }

          ldx = 0;
          for (i = 0; i < ixlen; i++)
            { if (fread(&idx,sizeof(int64),1,f) < 1)
                break;
              prefx[i] += idx-ldx;
              ldx = idx;
            }
          if (i < ixlen)
            { fprintf(stderr,"\n%s: Cannot read table %s/%s.ktab\n",Prog_Name,path,root);
              exit (1);
            }

          fclose(f);

          for (i = 1; i <= spart; i++)
            { f = fopen(Catenate(path,"/.",root,Numbered_Suffix(".ktab.",i,"")),"r");
              if (f == NULL)
                { fprintf(stderr,"\n%s: Cannot find and open table part %s/.%s.ktab.%d\n",
                                 Prog_Name,path,root,i);
                  exit (1);
                }
              fclose(f);
            }
          npart += spart;

          free(root);
          free(path);
        }

      for (i = 1; i < ixlen; i++)
        prefx[i] += prefx[i-1];

      command = Malloc(2*nmax+50,"Command buffer");
      if (command == NULL)
        exit (1);

      //  Create the target stub, being sure to back out if any error occurs

      f = fopen(Catenate(OPATH,"/",OROOT,".ktab"),"w");
      if (f == NULL)
        { fprintf(stderr,"\n%s: Cannot create and write table %s/%s.ktab\n",Prog_Name,OPATH,OROOT);
          exit (1);
        }

      error |= (fwrite(&kmer,sizeof(int),1,f) < 1);
      error |= (fwrite(&npart,sizeof(int),1,f) < 1);
      error |= (fwrite(&mcnt,sizeof(int),1,f) < 1);
      error |= (fwrite(&ibyte,sizeof(int),1,f) < 1);
      error |= (fwrite(prefx,sizeof(int64)*ixlen,1,f) < 1);
      fclose(f);

      if (error)
        { fprintf(stderr,"\n%s: Cannot write table %s/%s.ktab\n",Prog_Name,OPATH,OROOT);
          unlink(Catenate(OPATH,"/",OROOT,".ktab"));
          exit (1);
        }

      free(prefx);

      //  For each source part either move it or copy it to a target file with the
      //    catenated part #.  If any error occurs, remove the stub and any part files
      //    created to the point of the error.  For the destructive case, assume that if
      //    the first move worked then all subsequent moves will work.
    
      npart = 0;
      for (c = 0; c < narg; c++)
        { path  = PathTo(argv[c]);
          root = Root(argv[c],"");
          f = fopen(Catenate(path,"/",root,".ktab"),"r");

          fread(&smer,sizeof(int),1,f);
          fread(&spart,sizeof(int),1,f);

          fclose(f);

          if (VERBOSE)
            { if (KEEP)
                printf("  Copying %d parts of table %s\n",spart,root);
              else
                printf("  Moving %d parts of table %s\n",spart,root);
              fflush(stdout);
            }

          if (KEEP)
            for (i = 1; i <= spart; i++)
              { sprintf(command,"cp -f %s/.%s.ktab.%d %s/.%s.ktab.%d",
                                path,root,i,OPATH,OROOT,npart+i);
                if (system(command))
                  { fprintf(stderr,"\n%s: Could not copy %s/.%s.ktab.%d to %s/.%s.ktab.%d\n",
                                   Prog_Name,path,root,i,OPATH,OROOT,npart+i);
                    for (j = 1; j < npart+i; j++)
                      unlink(Catenate(OPATH,"/.",OROOT,Numbered_Suffix(".ktab.",j,"")));
                    unlink(Catenate(OPATH,"/",OROOT,".ktab"));
                    exit (1);
                  }
              }
          else
            { if (npart == 0)
                { sprintf(command,"mv -f %s/.%s.ktab.1 %s/.%s.ktab.1",path,root,OPATH,OROOT);
                  if (system(command))
                    { fprintf(stderr,"\n%s: Could not move %s/.%s.ktab.1 to %s/.%s.ktab.1\n",
                                     Prog_Name,path,root,OPATH,OROOT);
                      unlink(Catenate(OPATH,"/",OROOT,".ktab"));
                      exit (1);
                    }
                }
              for (i = 1 + (npart == 0); i <= spart; i++)
                { sprintf(command,"mv -f %s/.%s.ktab.%d %s/.%s.ktab.%d",
                                  path,root,i,OPATH,OROOT,npart+i);
                  system(command);
                }
            }
          npart += spart;

          free(root);
          free(path);
        }

      //  If a "destructive" cat, get rid of all the source stub files.

      if (!KEEP)
        { for (c = 0; c < narg; c++)
            { path = PathTo(argv[c]);
              root = Root(argv[c],"");
              unlink(Catenate(path,"/",root,".ktab"));
              free(root);
              free(path);
            }
        }
    }
    
  if (DO_HIST)
    { char      *path, *root;
      Histogram *H, *G;
      int        c, i;

      if (VERBOSE)
        { printf("\n  Cat'ing %d histograms\n",narg);
          fflush(stdout);
        }

      //  Open the first histogram in H

      H = Load_Histogram(argv[0]);
      if (H == NULL)
        { fprintf(stderr,"\n%s: Cannot open histogram %s\n",Prog_Name,argv[0]);
          exit (1);
        }

      //  Accumulate each subsequent histogram to H
 
      for (c = 1; c < narg; c++)
        { G = Load_Histogram(argv[c]);
          if (G == NULL)
            { fprintf(stderr,"\n%s: Cannot open histogram %s\n",Prog_Name,argv[c]);
              exit (1);
            }
          for (i = 0; i <= 0x8001; i++)
            H->hist[i] += G->hist[i];
          Free_Histogram(G);
        }

      //  Write the final histogram

      Write_Histogram(Catenate(OPATH,"/",OROOT,".hist"),H);
      Free_Histogram(H);

      //  If not -k then remove all the source histograms

      if (!KEEP)
        { for (c = 0; c < narg; c++)
            { path = PathTo(argv[c]);
              root = Root(argv[c],"");
              unlink(Catenate(path,"/",root,".hist"));
              free(root);
              free(path);
            }
        }
    }

  if (DO_PROF)
    { char *path, *root, *command;
      int   smer, kmer;
      int   spart, npart;
      int   sread, nread;
      int   nmax, error;
      FILE *f;
      int   i, j, c;

      // In a first pass read each source table stub:
      //   check that the kmer length are the same for all sources
      //   make sure each table part exists and you can open it
      //   determine the length of the longest source name (path+root)

      if (VERBOSE)
        { printf("\n  Cat'ing %d profiles\n",narg);
          fflush(stdout);
        }

      error = 0;
      kmer  = 0;
      npart = 0;
      nmax  = 0;
      for (c = 0; c < narg; c++)
        { path = PathTo(argv[c]);
          root = Root(argv[c],"");
          f = fopen(Catenate(path,"/",root,".prof"),"r");
          if (f == NULL)
            { fprintf(stderr,"\n%s: Cannot find and open index %s/%s.prof\n",Prog_Name,path,root);
              exit (1);
            }

          if ((int) (strlen(path)+strlen(root)) > nmax)
            nmax = strlen(path)+strlen(root);

          error |= (fread(&smer,sizeof(int),1,f) < 1);
          error |= (fread(&spart,sizeof(int),1,f) < 1);
          if (error)
            { fprintf(stderr,"\n%s: Cannot read profile %s/%s.prof\n",Prog_Name,path,root);
              exit (1);
            }

          if (VERBOSE)
            { printf("  Processing profile header of %s\n",root);
              fflush(stdout);
            }

          if (c == 0)
            kmer = smer;
          else
            { if (smer != kmer)
                { fprintf(stderr,"\n%s: K-mer tables do not involve the same K\n",Prog_Name);
                  exit (1);
                }
            }

          fclose(f);

          for (i = 1; i <= spart; i++)
            { f = fopen(Catenate(path,"/.",root,Numbered_Suffix(".pidx.",i,"")),"r");
              if (f == NULL)
                { fprintf(stderr,"\n%s: Cannot find and open profile part %s/.%s.pidx.%d\n",
                                 Prog_Name,path,root,i);
                  exit (1);
                }
              fclose(f);
              f = fopen(Catenate(path,"/.",root,Numbered_Suffix(".prof.",i,"")),"r");
              if (f == NULL)
                { fprintf(stderr,"\n%s: Cannot find and open profile part %s/.%s.prof.%d\n",
                                 Prog_Name,path,root,i);
                  exit (1);
                }
              fclose(f);
            }
          npart += spart;

          free(root);
          free(path);
        }

      command = Malloc(2*nmax+50,"Command buffer");
      if (command == NULL)
        exit (1);

      //  Create the target stub, being sure to back out if any error occurs

      f = fopen(Catenate(OPATH,"/",OROOT,".prof"),"w");
      if (f == NULL)
        { fprintf(stderr,"\n%s: Cannot create and write table %s/%s.prof\n",Prog_Name,OPATH,OROOT);
          exit (1);
        }

      error |= (fwrite(&kmer,sizeof(int),1,f) < 1);
      error |= (fwrite(&npart,sizeof(int),1,f) < 1);
      fclose(f);

      if (error)
        { fprintf(stderr,"\n%s: Cannot write profile %s/%s.prof\n",Prog_Name,OPATH,OROOT);
          unlink(Catenate(OPATH,"/",OROOT,".prof"));
          exit (1);
        }

      //  For each source part either move it or copy it to a target file with the
      //    catenated part #.  If any error occurs, remove the stub and any part files
      //    created to the point of the error.  For the destructive case, assume that if
      //    the first part moves worked then all subsequent moves will work.
    
      npart = 0;
      for (c = 0; c < narg; c++)
        { path  = PathTo(argv[c]);
          root = Root(argv[c],"");
          f = fopen(Catenate(path,"/",root,".ktab"),"r");

          fread(&smer,sizeof(int),1,f);
          fread(&spart,sizeof(int),1,f);

          fclose(f);

          if (VERBOSE)
            { if (KEEP)
                printf("  Copying %d parts of profile %s\n",spart,root);
              else
                printf("  Moving %d parts of profile %s\n",spart,root);
              fflush(stdout);
            }

          if (KEEP)
            for (i = 1; i <= spart; i++)
              { sprintf(command,"cp -f %s/.%s.pidx.%d %s/.%s.pidx.%d",
                                path,root,i,OPATH,OROOT,npart+i);
                if (system(command))
                  { fprintf(stderr,"\n%s: Could not copy %s/.%s.pidx.%d to %s/.%s.pidx.%d\n",
                                   Prog_Name,path,root,i,OPATH,OROOT,npart+i);
                    for (j = 1; j < npart+i; j++)
                      { unlink(Catenate(OPATH,"/.",OROOT,Numbered_Suffix(".pidx.",j,"")));
                        unlink(Catenate(OPATH,"/.",OROOT,Numbered_Suffix(".prof.",j,"")));
                      }
                    unlink(Catenate(OPATH,"/",OROOT,".prof"));
                    exit (1);
                  }
                sprintf(command,"cp -f %s/.%s.prof.%d %s/.%s.prof.%d",
                                path,root,i,OPATH,OROOT,npart+i);
                if (system(command))
                  { fprintf(stderr,"\n%s: Could not copy %s/.%s.prof.%d to %s/.%s.prof.%d\n",
                                   Prog_Name,path,root,i,OPATH,OROOT,npart+i);
                    for (j = 1; j < npart+i; j++)
                      { unlink(Catenate(OPATH,"/.",OROOT,Numbered_Suffix(".pidx.",j,"")));
                        unlink(Catenate(OPATH,"/.",OROOT,Numbered_Suffix(".prof.",j,"")));
                      }
                    unlink(Catenate(OPATH,"/.",OROOT,Numbered_Suffix(".pidx.",npart+i,"")));
                    unlink(Catenate(OPATH,"/",OROOT,".prof"));
                    exit (1);
                  }
              }
          else
            { if (npart == 0)
                { sprintf(command,"mv -f %s/.%s.pidx.1 %s/.%s.pidx.1",path,root,OPATH,OROOT);
                  if (system(command))
                    { fprintf(stderr,"\n%s: Could not move %s/.%s.pidx.1 to %s/.%s.pidx.1\n",
                                     Prog_Name,path,root,OPATH,OROOT);
                      unlink(Catenate(OPATH,"/",OROOT,".prof"));
                      exit (1);
                    }
                  sprintf(command,"mv -f %s/.%s.prof.1 %s/.%s.prof.1",path,root,OPATH,OROOT);
                  if (system(command))
                    { fprintf(stderr,"\n%s: Could not move %s/.%s.prof.1 to %s/.%s.prof.1\n",
                                     Prog_Name,path,root,OPATH,OROOT);
                      unlink(Catenate(OPATH,"/.",OROOT,".prof.1"));
                      unlink(Catenate(OPATH,"/",OROOT,".prof"));
                      exit (1);
                    }
                }
              for (i = 1 + (npart == 0); i <= spart; i++)
                { sprintf(command,"mv -f %s/.%s.pidx.%d %s/.%s.pidx.%d",
                                  path,root,i,OPATH,OROOT,npart+i);
                  system(command);
                  sprintf(command,"mv -f %s/.%s.prof.%d %s/.%s.prof.%d",
                                  path,root,i,OPATH,OROOT,npart+i);
                  system(command);
                }
            }
          npart += spart;

          free(root);
          free(path);
        }

      nread = 0;
      for (i = 1; i <= npart; i++)
        { int g = open(Catenate(OPATH,"/.",OROOT,Numbered_Suffix(".pidx.",i,"")),O_RDWR);
          read(g,&smer,sizeof(int));
          read(g,&sread,sizeof(int64));
          read(g,&sread,sizeof(int64));
          lseek(g,sizeof(int),SEEK_SET);
          write(g,&nread,sizeof(int64));
          close(g);
          nread += sread;
        }
        
      //  If a "destructive" cat, get rid of all the source stub files.

      if (!KEEP)
        { for (c = 0; c < narg; c++)
            { path = PathTo(argv[c]);
              root = Root(argv[c],"");
              unlink(Catenate(path,"/",root,".prof"));
              free(root);
              free(path);
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
