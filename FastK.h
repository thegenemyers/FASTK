 /********************************************************************************************
 *
 *  FastK: a rapid disk-based k-mer counter for high-fidelity shotgun data sets.
 *     Uses a novel minimizer-based distribution scheme that permits problems of
 *     arbitrary size, and a two-staged "super-mer then weighted k-mer" sort to acheive
 *     greater speed when error rates are low (1% or less).  Directly produces read
 *     profiles.
 *
 *  Author:  Gene Myers
 *  Date  :  October, 2020
 *
 *********************************************************************************************/

#ifndef _KMERS
#define _KMERS

#undef DEVELOPER

#define IO_BUF_LEN   4096       // number of uint's in bit stuffed IO buffer for each part+thread
#define IO_UBITS       64
#define IO_UBYTES       8
#define IO_UTYPE   uint64

#define NPANELS 4

#ifdef DEVELOPER

extern int DO_STAGE;  // Stage to run if code development (not for users)

#endif

  //  Option Settings

extern int    VERBOSE;     //  show progress
extern int64  SORT_MEMORY; // Memory available for each k-mer sort 
extern int    KMER;        //  desired K-mer length
extern int    NTHREADS;    //  # of threads to run with
extern int      ITHREADS;    //  # of threads possible for input
extern char  *SORT_PATH;   //  where to put external files

extern int    DO_TABLE;    // Zero or table cutoff
extern int    DO_PROFILE;  // Do or not
extern Kmer_Stream *PRO_TABLE;   //  Kmer stream of profile option (only if relative profile)
extern char        *PRO_NAME;    //  Name of profile table
extern int    BC_PREFIX;   // Ignore prefix of each read of this length
extern int    COMPRESS;    // Homopolymer compress the input


  //  Sizes and numbers of items (k-mers, super-mers, reads, positions)

extern int    NPARTS;      //  number of k-mer buckets
extern int    SMER;        //  max size of a super-mer (= MAX_SUPER + KMER - 1)ZZ
extern int64  KMAX;        //  max k-mers in any part
extern int64  NMAX;        //  max super-mers in any part

extern int    MOD_LEN;     //  length of minimizer buffer (power of 2)
extern int    MOD_MSK;     //  mask for minimzer buffer

extern int    MAX_SUPER;   //  = (KMER - PAD_LEN) + 1

extern int    SLEN_BITS;      //  # of bits needed to encode the length of a super-mer (less KMER-1)
extern uint64 SLEN_BIT_MASK;  //  Bit-mask for super-mer lengths

extern int RUN_BITS;     //  # of bits encoding largest run index
extern int RUN_BYTES;    //  # of bytes encoding largest run index

extern int SLEN_BYTES;   //  # of bytes encoding Length of a super-mer
extern int PLEN_BYTES;   //  # of bytes encoding length of a compressed profile segment
extern int PROF_BYTES;   //  # of bytes encoding RUN+PLEN_BYTES
extern int IDX_BYTES;    //  # of bytes for table prefix index

extern int KMER_BYTES;   //  # of bytes encoding a KMER 
extern int SMER_BYTES;   //  # of bytes encoding a super-mer 
extern int KMAX_BYTES;   //  # of bytes to encode an int in [0,KMAX)

extern int KMER_WORD;    //  # of bytes in a k-mer entry
extern int SMER_WORD;    //  # of bytes in a super-mer entry
extern int TMER_WORD;    //  bytes to hold a k-mer/count table entry
extern int CMER_WORD;    //  bytes to hold a count/position entry

extern int64 *NUM_RID;   //  [i] for i in [0,ITHREADS) = # of super-mers per vertical stripe

extern uint8 Comp[256];  //  complement of 4bp byte code

  //  IO Module Interface

typedef struct
  { int64       totlen;   //  total # of bases in data set
    int         nreads;   //  # of reads in data set
    double      ratio;    //  ratio of file size to portion read (first block only)
    int64       maxbps;   //  size of bases array
    int         maxrds;   //  size of boff array
    char       *bases;    //  concatenation of 0-terminated read strings
    int        *boff;     //  read i is at bases+boff[i], boff[n] = total bytes
    int         rem;      //  Length of remainder of current sequence to process
    char       *next;     //  Remainder of current input sequence (overlaps by KMER-1) with last
                          //     sequence in this block buffer now.
  } DATA_BLOCK;

typedef void *Input_Partition;

  Input_Partition *Partition_Input(int argc, char *argv[]);

  DATA_BLOCK *Get_First_Block(Input_Partition *part, int64 numbps);

  void Free_First_Block(DATA_BLOCK *block);

  char *First_Root(Input_Partition *part);
  char *First_Pwd (Input_Partition *part);

  void Scan_All_Input(Input_Partition *part);

  void Free_Input_Partition(Input_Partition *part);

  //  Stages

void Clean_Exit(int status);

int Determine_Scheme(DATA_BLOCK *block);

void Split_Kmers(Input_Partition *io, char *root);

  void Distribute_Block(DATA_BLOCK *block, int tid);

void Split_Table(char *root);

void Sorting(char *path, char *root);

void Merge_Tables(char *path, char *root);

void Merge_Profiles(char *path, char *root);

  //  Sorts

typedef struct
  { int    beg;     //  [beg,end) first byte with sorted data at off
    int    end;
    int64  off;
    int64  khist[256];
    int64  count[0x8000];
    int64  max_inst;
    int    byte1;   //  used internallly by sort
  } Range;

void Supermer_Sort(uint8 *array, int64 nelem, int rsize, int ksize,
                   int64 *part, int nthreads, Range *panels);

void Weighted_Kmer_Sort(uint8 *array, int64 nelem, int rsize, int ksize,
                        int64 *part, int nthreads, Range *panels);

void MSD_Sort(uint8 *array, int64 nelem, int rsize, int ksize,
              int64 *part, int nthreads, Range *panels);

void *LSD_Sort(int64 nelem, void *src, void *trg, int rsize, int *bytes);

#endif // _KMERS
