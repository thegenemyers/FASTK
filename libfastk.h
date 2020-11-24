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
#include <ctype.h>
#include <stdint.h>
#include <math.h>
#include <dirent.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <errno.h>

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
  { int     kmer;         //  Kmer length
    int     minval;       //  The minimum count of a k-mer in the table
    int     kbyte;        //  Kmer encoding in bytes
    int     tbyte;        //  Kmer+count entry in bytes
    int64   nels;         //  # of unique, sorted k-mers in the table
    uint8  *table;        //  The (huge) table in memory
    void   *private[1];   //  Private fields
  } Kmer_Table;

Kmer_Table *Load_Kmer_Table(char *name, int cut_off);
void        Free_Kmer_Table(Kmer_Table *T);

char       *Fetch_Kmer(Kmer_Table *T, int i);
int         Fetch_Count(Kmer_Table *T, int i);

int         Find_Kmer(Kmer_Table *T, char *kseq);

void        List_Kmer_Table(Kmer_Table *T, FILE *out);
int         Check_Kmer_Table(Kmer_Table *T);


  //  K-MER STREAM

typedef struct
  { int    kmer;       //  Kmer length
    int    minval;     //  The minimum count of a k-mer in the stream
    int    kbyte;      //  Kmer encoding in bytes
    int    tbyte;      //  Kmer+count entry in bytes
    int64  nels;       //  # of elements in entire table
    void  *private[8]; //  Private fields
  } Kmer_Stream;

Kmer_Stream *Open_Kmer_Stream(char *name, int cut_freq);
void         Free_Kmer_Stream(Kmer_Stream *S);

uint8       *First_Kmer_Entry(Kmer_Stream *S);
uint8       *Next_Kmer_Entry(Kmer_Stream *S);

char        *Current_Kmer(Kmer_Stream *entry);
int          Current_Count(Kmer_Stream *entry);


  //  PROFILES

typedef struct
  { int    kmer;     //  Kmer length
    int    nparts;   //  # of threads/parts for the profiles
    int    nreads;   //  total # of reads in data set
    int64 *nbase;    //  nbase[i] for i in [0,nparts) = id of last read in part i + 1
    int   *nfile;    //  nfile[i] for i in [0,nparts) = stream for ".prof" file of part i
    int64 *index;    //  index[i] for i in [0,nreads) = offset in relevant part of
                     //    compressed profile for read i
  } Profile_Index;

Profile_Index *Open_Profiles(char *name);

void Free_Profiles(Profile_Index *P);

int Fetch_Profile(Profile_Index *P, int64 id, int plen, uint16 *profile);

#endif // _LIBFASTK
