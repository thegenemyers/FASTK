# FastK: A K-mer counter (for HQ assembly data sets)
  
## _Author:  Gene Myers_
## _First:   July 22, 2020_
## _Current: October 20, 2020_

FastK is a k-mer counter that is optimized for processing high quality DNA assembly data
sets such as those produced with an Illumina instrument or a PacBio run in HiFi mode.
For example it is about 2 times faster than KMC3 when counting 40-mers in a 50X HiFi data
set.  Its relative speedup decreases with increasing error rate or increasing values of k,
but regardless is a general program that works for any DNA sequence data set and choice of k.
It is further designed to handle data sets of arbitrarily large size, e.g. a 100X data
set of a 32GB Axolotl genome (3.2Tbp) can be performed on a machine with just 12GB of memory
provided it has ~6.5TB of disk space available.

FastK can produce the following outputs:

1. a histogram of the frequency with which each k-mer in the data set occurs.
2. a table of k-mer/count pairs sorted lexicographically on the k-mer where a < c < g < t.
3. a k-mer count profile of every sequence in the data set.  A **profile** is the sequence of counts of the n-(k-1) consecutive k-mers of a sequence of length n. 

Note carefully, that in order to accommodate the unknown orientation of a sequencing read,
a k-mer and its Watson Crick complement are considered to be the same k-mer by FastK, where the
lexicograpahically smaller of the two alternatives is termed **canonical**.
The histogram is always produced whereas the production of a k-mer table (2.) and profiles (3.)
are controlled by command
line options.  The table (2.) is over just the *canonical* k-mers present in the data set.  Producing profiles (3.) as part of the underlying sort is much more efficient than producing them after the fact using a table or hash of all k-mers such as is necessitated when using other k-mer counter programs.  The profiles are recorded in a space-efficient compressed form, e.g.
about 4.7-bits per base for a recent 50X HiFi asssembly data set.

```
1. FastK [-k<int(40)] [-h[<int(1)>:]<int>] [-t<int(3)>] [-cp] [-bc<int>]
          [-v] [-N<path_name>] [-P<dir(/tmp)>] [-M<int(12)>] [-T<int(4)>]
            <source>[.cram|.[bs]am|.db|.dam|.f[ast][aq][.gz]] ...
```

FastK counts the number of k-mers in a corpus of DNA sequences over the alphabet {a,c,g,t} for a specified k-mer size, 40 by default.
The input data can be in one or more CRAM, BAM, SAM, fasta, or fastq files, where the later two can be gzip'd, preferrably
by [VGPzip](https://github.com/VGP/vgp-tools/tree/master/VGP). The data can also be in [Dazzler databases](https://github.com/thegenemyers/DAZZ_DB).  The type of the
file is determined by its extension (and not its
contents).  The extension need not be given if the root name suffices
to uniquely identify a file.  If more than one source file is given
they must all be of the same type in the current implementation.

FastK produces a number of outputs depending on the setting of its options.  By default, the
outputs will be placed in the same directory as that of the first input and begin with the
prefix \<source> which is the first path name absent any suffix extensions.  For example, if
the input is <code>../BLUE/foo.fastq</code> then \<source> is <code>../BLUE/foo</code>, the
outputs will be placed in directory <code>../BLUE</code>, and all result file names will begin with <code>foo</code>.  If the -N option is specified then the path name specified is used
as \<source>.

One can select any value of k &ge; 5 with the -k option.
FastK always outputs a file <code>\<source>.hist</code> that contains a histogram of the k-mer frequency
distribution where the highest possible count is 2<sup>15</sup>-1 = 32,767 -- FastK clips all higher values to this upper limit.  Its exact format is described in the section on Data Encodings.
If the -h option is specified then the histogram is displayed by FastK on the standard output as 
soon as it is computed.
This option allows you to specify the interval of frequencies to display, where the lower limit
is 1 if ommitted.

One can optionally request, by specifying the -t option, that FastK produce a sorted table of
all canonical k-mers that occur -t or more times along with their counts, where the default
for the threshold is 3.
The output is placed in a single *stub* file with the name <code>\<source>.ktab</code> and N
roughly equal-sized *hidden* files with the names <code>.\<base>.ktab.#</code> assuming
\<source> = \<dir>/\<base> and
where # is a thread number between 1 and N where N is the number of threads used by FastK (4 by default).
The exact format of the N-part table is described in the section on Data Encodings.

One can also ask FastK to produce a k-mer count profile of each sequence in the input data set
by specifying the -p option.  A single *stub* file with the name <code>\<source>.prof</code> is output
along with 2N roughly equal-sized pairs of *hidden* files with the names
<code>.\<base>.pidx.#</code> and <code>.\<base>.prof.#</code> in the order of the sequences in the input.  Here \<base> is the base name part of \<source> = \<dir>/\<base>.
The profiles are individually compressed and the exact format of these
files is described in the section on Data Encodings.

The -c option asks FastK to first homopolymer compress the input sequences before analyzing
the k-mer content.  In a homopolymer compressed sequence, every substring of 2 or more a's
is replaced with a single a, and similarly for runs of c's, g's, and t's.  This is particularly useful for Pacbio data where homopolymer errors are five-fold more frequent than other
errors and thus the error rate of these "hoco" k-mers is five-fold less.

The -v option asks FastK to output information about its ongoing operation to standard error.
The -bc option allows you to ignore the prefix of each read of the indicated length, e.g. when
the reads have a bar code at the start of each read.
The -P option specifies where FastK should place all the numerous temporary files it creates, if not `/tmp` by default.
The -M option specifies the maximum amount of memory, in GB, FastK should use at any given
moment.
FastK by design uses a modest amount of memory, the default 12GB should generally
be more than enough.
Lastly, the -T option allows the user to specify the number of threads to use.
Generally, this is ideally set to the actual number of physical cores in one's machine.
            
```
2a. Fastrm [-i] <source>[.hist|.ktab|.prof] ...
2b. Fastmv [-in] <source>[.hist|.ktab|.prof] <dest>
2c. Fastcp [-in] <source>[.hist|.ktab|.prof] <dest>
```

As described above FastK produces hidden files whose names begin with a . for the -t and -p
options in order to avoid clutter when listing a directory's contents.
An issue with this approach is that it is inconvenient for the user to remove, rename, or copy these files
and often a user will forget they are there, potentially wasting disk space.
We therefore provide Fastrm, Fastmv, and Fastcp that remove, rename, and copy FastK .hist, .ktab, and .prof output files as a single unit.

If \<source> does not end with a FastK extenion then the command operates on any histogram, k-mer table, and profile files with \<source> as its prefix.  Otherwise the command operates on the file with the given extension and its hidden files.  Fastrm removes the relevant stub and
hidden files, Fastmv renames all the relevant files as if FastK had been called with option -N\<dest>, and Fastcp makes a copy of all associated files with the path name \<dest>.  If \<dest> is a directory than
the base name of source is used to form a complete destination path for both Fastmv and Fastcp.

As for the UNIX rm, mv, and cp commands, the -i option asks the command to query each file as to whether you want to delete (rm) or overwrite (mv,cp) it, but only for the stubs and not the hidden files corresponding to each stub, which share the same fate as their stub file.
The -n option tells Fastmv and Fastcp to not overwrite any files.

```
3. Histex [-h[<int(1)>:]<int(100)>] <source>[.hist]
```

This command and also Tabex and Profex are presented specifically to
give a user an example of how to use the C-interface, <code>libfastk.c</code>, module
to manipulate the histogram, k-mer count tables, and profiles produced by FastK.

Given a histogram file \<source>.hist produced by FastK,
one can view the histogram of k-mer counts with **Histex** where the -h specifies the 
interval of frequencies to be displayed where 1 is assumed if the lower bound is not given.

```
4. Tabex [-t<int>] <source>[.ktab]  (LIST|CHECK|(<k-mer:string>) ...
```

Given that a set of k-mer counter table files have been generated represented by stub file
\<source>.ktab, ***Tabex*** opens the corresponding hidden table files (one per thread) and then performs the sequence of actions specified by the
remaining arguments on the command line.  The literal argument LIST lists the contents
of the table in radix order.  CHECK checks that the table is indeed sorted.  Otherwise the
argument is interpreted as a k-mer and it is looked up in the table and its count returned
if found.  If the -t option is given than only those k-mers with counts greater or equal to the given value are operated upon.

```
5. Profex <source>[.prof] <read:int> ...
```
Given that a set of profile files have been generated represented by stub file
\<source>.prof, ***Profex*** opens the corresonding hiddent profile files (two per thread)
and gives a display of each sequence profile whose ordinal id is given on
the remainder of the command line.  The index of the first read is 1 (not 0).

## Current Limitations

Currently if multiple input files are given they must all be of the same type, e.g. fasta
or cram.  This restrictionis is not fundamental and could be removed with some coding effort.

## The FastK C-interface to K-mer Histograms, Tables, and Profiles

For each of the 3 distinct outputs of FastK, we have suppled a simple C library
that gives a user access to the data therein.  The library is simply embodied in
the C-file, <code>libfastk.c</code>, and associdated include file <code>libfastk.h</code>.
The makefile commands for building Histex, Tabex, and Profex illustrate how to
easily incorporate the library into your C or C++ code.

### K-mer Histogram

```
typedef struct
  { int    kmer;  //  Histogram is for k-mers of this length
    int    low;   //  Histogram is for range [low,hgh]
    int    high;
    int64 *hist;  //  hist[i] for i in [low,high] = # of k-mers occuring i times
  } Histogram;

Histogram *Load_Histogram(char *name, int low, int high);
void       Free_Histogram(Histogram *H);
```

### K-mer Table

```
typedef struct
  { int     kmer;    //  Kmer length
    int     kbyte;   //  Kmer encoding in bytes
    int     tbyte;   //  Kmer+count entry in bytes
    int64   nels;    //  # of unique, sorted k-mers in the table
    uint8  *table;   //  The (huge) table in memory
    int64  *index;   //  Search accelerator if needed
  } Kmer_Table;

Kmer_Table *Load_Kmer_Table(char *name, int cut_freq, int smer, int nthreads);
void        Free_Kmer_Table(Kmer_Table *T);

char       *Fetch_Kmer(Kmer_Table *T, int i);
int         Fetch_Count(Kmer_Table *T, int i);

void        List_Kmer_Table(Kmer_Table *T);
void        Check_Kmer_Table(Kmer_Table *T);
int         Find_Kmer(Kmer_Table *T, char *kseq);
```

### Sequence Profiles

```
typedef struct
  { int    kmer;     //  Kmer length
    int    nparts;   //  # of threads/parts for the profiles
    int    nreads;   //  # of threads/parts for the profiles
    int64 *nbase;    //  nbase[i] = id of last read in part i
    FILE **nfile;    //  nfile[i] = stream for "P" file of part i
    int64 *index;    //  Kmer+count entry in bytes
  } Profile_Index;

Profile_Index *Open_Profiles(char *name, int kmer, int nthreads);

void Free_Profiles(Profile_Index *P);

int Fetch_Profile(Profile_Index *P, int64 id, int plen, uint16 *profile);
```

## Data Encodings

### K-mer Histogram

The histogram file has a name of the form <code>\<source>.hist</code> where \<source> is the
output path name used by FastK.  It contains an initial integer
giving the k-mer size <k>, followed by two integers giving the range \[l,h] (inclusive) of frequencies in the histogram, followed by (h-l)+1 64-bit counts for frequencies l, l+1, ..., h.
Formally,

```
    < kmer size(k)  : int >
    < 1st freq.(l)  : int >
    < last freq.(h) : int >
    ( < k-mer count for f=l, l+1, ... h : int64 > ) ^ (h-l)+1
```

FastK always outputs a histogram with l=1 and h=32,637, and so the file is exactly 262,148
bytes in size.  Other auxiliary programs can produce histograms over a subrange of this range.
Note carefully that the count for the entry h, is actually the count of all k-mers that occur
h-or-more times, and when l>1, the entry l, is the count of all k-mers that occur l-or-fewer
times.

### K-mer Table

A table of canonical k-mers and their counts is produced in N parts, where N is the number of threads FastK was run with.
A single *stub* file <code>\<source>.ktab</code> where <source> is the output path name used by
FastK.
This stub file contains just the k-mer length followed by the number of threads FastK was run
with as two integers.
The table file parts are in N hidden files in the same directory as the stub
file with the names <code>.\<base>.ktab.[1,N]</code> assuming that \<source> = \<dir>/\<base>.
The k-mers in each part are lexicographically ordered and the k-mers in Ti are all less than the k-mers in T(i+1), i.e. the concatention of the N files in order of thread index is sorted.  

The data in each table file is as follows:

```
    < kmer size(k)   : int   >
    < # of k-mers(n) : int64 >
    ( < bit-encoded k-mer : uint8 ^ (k+3)/4 > < count : uint16 ) ^ n
```
    
In words, an initial integer gives the kmer size that FastK was run with followed by
a 64-bit int, n, that gives the # of entries in the file.  The remainder of the file is then n
k-mer,cnt pairs.  The k-mer is encoded in (k+3)/4 bytes where each base is compressed into 2-bits so that each byte contains up to four bases, in order of high bits to low bits.  The
bases a,c,g,t are assigned to the values 0,1,2,3, respectively.  As an example, 0xc6 encodes
tacg.  The last byte is partially filled if k is not a multiple of 4, and the remainder is
guaranteed to be zeroed.  The byte sequence for a k-mer is then followed by a 2-byte 
unsigned integer count with a maximum value of 32,767.

### Sequence Profiles

The read profiles are stored in N pairs of file, an index and a data pair, that are hidden
and identified by a single *stub* file <code>\<source>.prof</code>.
This stub file contains just the k-mer length followed by the number of threads FastK was
run with as two integers.
The hidden data files, <code>.\<base>.prof.[1,N]</code>, contain the compressed profiles for
each read
in their order in the input data set, and the hidden index files,
<code>.\<base>.pidx.[1,N]</code>, 
contain arrays of offsets into the P-files giving the start of each compressed profile,
assuming the path name \<source> = \<dir>/\<base>.
An A-file contains a brief header followed by an array of offsets.
The number of offsets in the A-file is one more than the number of profiles.
This last offset is to the end of the P-file so that the profile for
sequence b+i is the bytes off[i] to off[i+1]-1 where off[i] is the i'th offset in the A-file.

```
      < kmer size(k)                                     : int   >
      < index of sequence of 1st profile in this file(b) : int64 >
      < # of profile offsets in this file(n)             : int64 >
      ( < profile offset for sequence i in [b,b+n) : int64 > ) ^ n+1
```

A P-file contains compressed profiles.
The sequence of a profile is given by the first count followed by the first forward difference to
each successive count as these are expected to be small integers that will compress well.
The the first count is encoded in the first one or two bytes, depending on
its value, as follows:

```
  0x   => x+1 in [1,128]
  1x,y => x.y in [129,32767]
```
That is, if the high order bit is set then the count is the unsigned integer encoded in the
remaining 15-bits of the firsts two bytes.  If it is not set, then it is one more than the
unsigned integer encoded in the remaining 7-bits of the first byte.  The one byte encoding is
used whenever possible.

All subsequent counts are given as first forward differences to the immediately preceeding
count, and are encoded in one or two bytes.  Often there is a run of 0 differences in a
profile and a special one byte form codes this as a run length.  Formally the forward differences are encoded as follows:

```
  00x  => x+1 in [1,64] 0 diffs
  010x =>   x+1  in [1,32]
  011x => -(x+1) in [-1,-32]
  1x,y => x.y in 15 bit 2's complement in [-16384,-33] U [33,16383] (modulo 32768 arithmetic)
```

If the high order bit of the current byte is set, then the remaining 15-bits of the current
byte and the one following it are interpreted as a *signed* 2's complement integer.
Moreover, the forward difference is taken module 32768 (2^15), e.g. 32766 + 3 = 1.
If the 2 highest order bits of the current byte are zero, then the next x+1 counts are the
same as the most recent count, where x is the lower 6-bits interpreted as an unsigned integer.
If the 2 highest order bit of the current byte are 01, then the remaining 6 bits are interpreted
as a *1's complement* integer and the difference is one more or less than said value depending
on the sign.  The single byte encoding is used whenever possible.