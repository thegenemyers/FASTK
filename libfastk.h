/*******************************************************************************************
 *
 *  C library routines to access and operate upon FastK histogram, k-mer tables, and profiles
 *
 *  Author:  Gene Myers
 *  Date  :  November 2020
 *
 *******************************************************************************************/

#ifndef _LIBFASTK
#define _LIBFASTK

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <unistd.h>
#include <dirent.h>

#include "gene_core.h"

  //  HISTOGRAM

typedef struct
  { int    kmer;  //  Histogram is for k-mers of this length
    int    low;   //  Histogram is for range [low,hgh]
    int    high;
    int64 *hist;  //  hist[i] for i in [low,high] = # of k-mers occuring i times
  } Histogram;

Histogram *Load_Histogram(char *name);
void       Subrange_Histogram(Histogram *H, int low, int high);
void       Free_Histogram(Histogram *H);


  //  K-MER TABLE

typedef struct
  { int     kmer;    //  Kmer length
    int     kbyte;   //  Kmer encoding in bytes
    int     tbyte;   //  Kmer+count entry in bytes
    int64   nels;    //  # of unique, sorted k-mers in the table
    uint8  *table;   //  The (huge) table in memory
    int64  *index;   //  Search accelerator if needed
  } Kmer_Table;

Kmer_Table *Load_Kmer_Table(char *name);
void        Cut_Kmer_Table(Kmer_Table *t, int cut_freq);
void        Free_Kmer_Table(Kmer_Table *T);

char       *Fetch_Kmer(Kmer_Table *T, int i);
int         Fetch_Count(Kmer_Table *T, int i);

int         Find_Kmer(Kmer_Table *T, char *kseq);

void        List_Kmer_Table(Kmer_Table *T, FILE *out);
int         Check_Kmer_Table(Kmer_Table *T);


   //  PROFILES

typedef struct
  { int    kmer;     //  Kmer length
    int    nparts;   //  # of threads/parts for the profiles
    int    nreads;   //  total # of reads in data set
    int64 *nbase;    //  nbase[i] for i in [0,nparts) = id of last read in part i + 1
    FILE **nfile;    //  nfile[i] for i in [0,nparts) = stream for "P" file of part i
    int64 *index;    //  index[i] for i in [0,nreads) = offset in relevant part of
                     //    compressed profile for read i
  } Profile_Index;

Profile_Index *Open_Profiles(char *name);

void Free_Profiles(Profile_Index *P);

int Fetch_Profile(Profile_Index *P, int64 id, int plen, uint16 *profile);

#endif // _LIBFASTK
