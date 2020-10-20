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
The histogram is always produced whereas the production of 2. and 3. are controlled by command
line options.  The table (2.) is over just the canonical k-mers present in the data set.  Producing profiles (3.) as part of the underlying sort is much more efficient than producing them after the fact using a table or hash of all k-mers such as is necessitated when using other k-mer counter programs.  The profiles are recorded in a space-efficient compressed form, e.g.
about 4.7-bits per base for a recent 50X HiFi asssembly data set.

```
FastK [-k<int(40)] [-h[<int(1)>:]<int>] [-t<int(3)>] [-p]
        [-v] [-P<dir(/tmp)>] [-M<int(12)>] [-T<int(4)>]
          <source:cram|[bs]am|f[ast][aq][.gz]|db|dam> ...
```

FastK counts the number of k-mers in a corpus of DNA sequences, \<source> ,over the alphabet {a,c,g,t} for a specified k-mer size, 40 by default.  The input data can be in one or more CRAM, BAM, SAM,
fasta, or fastq files, where the later two can be gzip'd, preferrably
by [VGPzip](https://github.com/VGP/vgp-tools/tree/master/VGP). The data can also be in [Dazzler databases](https://github.com/thegenemyers/DAZZ_DB).  The type of the
file is determined by its extension on the command line (and not its
contents).  The extension need not be given if the root name suffices
to uniquely identify a file.  If more than one source file is given
they must all be of the same type in the current implementation.

One can select any value of k &ge; 5 with the -k option.
FastK always outputs a file \<source_root>.K\<k> that contains a histogram of the k-mer frequency
distribution where the highest possible count is 2<sup>15</sup>-1 = 32,767 -- FastK clips all higher values to this upper limit.  The string \<source\_root> is the source file name of the
first file given without any extensions.  The file is placed in the same directory as the input
file and its exact format is described in the section on Data Encodings.
If the -h option is specified then histogram is displayed by FastK on the standard output as 
soon as it is computed.
The option allows you to specify the interval of frequencies to display, where the lower end
is 1 if ommitted.

One can optionally request, by specifying the -t option, that FastK produce a sorted table of
all canonical k-mers that occur -t or more times along with their counts, where the default
for the threshold is 3.
The output is placed in N equal sized files with the names \<source\_root>.K\<k>.T# where \<k> is the
specified k-mer size and # is a thread
number between 1 and N where N is the number of threads used by FastK (4 by default).
The files are place in the same directory as the input.
The exact format of the N-part table is described in the section on Data Encodings.

One can also ask FastK to produce a k-mer count profile of each sequence in the input data set
by specifying the -p option.  The profiles are individually compressed and placed in N roughly
equal sized pairs of files with the names
\<source\_root>.K\<k>.A# and \<source\_root>.K\<k>.P# in the order of the sequences in the input.  The files are placed in the same directory as the input.  The exact format of these
files is described in the section on Data Encodings.

The -v option asks FastK to output information about its ongoing operation to standard error.
The -P option specifies where FastK should place all the numerous temporary files it creates, if not `/tmp` by default.
The -M options specifies the maximum amount of memory, in GB, FastK should use at any given
moment.
FastK by design uses a modest amount of memory, the default 12GB should generally
be more than enough.
Lastly, the -T option allows the user to specify the number of threads to use.
Generally, this is ideally set to the actual number of physical cores in one's machine.

```
2. Histex [-h[<int(1)>:]<int(100)>] <source_root>.K<k>
```

This command and also Tabex and Profex are presented specifically to
give a user an example and potentially code they can use to manipulate
the histogram, k-mer count tables, and profiles produced by FastK.

Given a histogram file \<data>.K# produced by FastK,
one can view the histogram of k-mer counts with **Histex** where the -h specifies the 
interval of frequencies do be displayed where 1 is assumed if the lower bound is not given.

```
3. Tabex [-t<int>] <source_root>.K<k>  (LIST|CHECK|(<k-mer:string>) ...
```

Given that a set of k-mer counter table files for k-mer size \<k> have been generated, ***Tabex*** opens the tables (one per thread) and then performs the sequence of actions specified by the
remaining arguments on the command line.  The literal argument LIST lists the contents
of the table in radix order.  CHECK checks that the table is indeed sorted.  Otherwise the
argument is interpreted as a k-mer and it is looked up in the table and its count returned
if found.  If the -t option is given than only those k-mers with counts greater or equal to the given value are operated upon.

```
4. Profex <source_root>.K<k> <read:int> ...
```
Given that a set of profile files for k = \<k> have been generated, ***Profex*** opens the profile filles (one per thread) and gives a
display of each sequence profile whose ordinal id is given on
the remainder of the command line.  The index of the first read is 1 (not 0).

## Current Limitations

Currently if multiple input files are given they must all be of the same type, e.g. fasta
or cram.  This restriction could be removed with some more code.

Currently the maximum size of any read/sequence is 1Mbp.  This limit will be problematic
if one is for example counting k-mers in a high-quality assembled genome where the contigs
exceed this limit.  Again the restriction could be removed but is quite complex if only a limited amount of memory is to be consumed.

Lastly, very small k-mer sizes may not work for large data sets.  The absolute minimum is 5,
but the core prefix trie used for distributing k-mers may want to use a larger minimizer
length than 5, and the minimizer length must be less than the k-mer length.  The logic
here needs to be investigated so that one has at least a dependable lower bound on k-mer
length.

## Data Encodings

### K-mer Histogram

The histogram file has a name of the form \<source\_root>.K\<k> where \<source\_root> is the root name of the input.  It contains an initial integer
giving the k-mer size <k>, followed by 32,767 64-bit counts for frequencies 1, 2, 3 ... so it
is precisely 262,140 bytes in size.
Formally,

```
    < kmer size(k) : int >
    ( < k-mer count for f=1,2,3, ... : int64 > ) ^ 32,767
```

### K-mer Table

A table of canonical k-mers and their counts is produced in N parts, where N is the number of threads FastK was run with.  The files are placed in the same directory as the input and named
\<source\_root>.K\<k>.T[1,N] where \<source\_root> is the root name of the input.  The k-mers in each part are lexicographically ordered and the k-mers in Ti are all less than the k-mers in T(i+1), i.e. the concatention of the N files in order of thread index is sorted.  

The data in each table file is as follows:

```
    < kmer size(k)   : int   >
    < # of k-mers(n) : int64 >
    ( < bit-encoded k-mer : uint8 ^ (k+3)/4 > < cnt : int16> ) ^ n
```
    
In words, an intial 64-bit int, n, gives the # of entries in the file that is followed by an int, k, giving the k-mer size that FastK was run with.  The remainder of the file is then n
k-mer,cnt pairs.  The k-mer is encoded in (k+3)/4 bytes where each base is compressed into 2-bits so thata each byte contains up to four bases, in order of high bits to low bits.  The
bases a,c,g,t are assigned to the values 0,1,2,3, respectively.  As an example, 0xc6 encodes
tacg.  The last byte is partially filled if k is not a multiple of 4.  The byte sequence for a
k-mer is then followed by a 2-byte integer count with a maximumb value of 32,767.

### Sequence Profiles

The files with suffixes P# are compressed profiles for each sequence in the input data set,
and the files with the suffixes A# are arrays of offsets into the P-files giving the start of
each compressed profile.  The number of offsets in the A-file is one more than the number of
profiles in the corresponding P-fille.  The A-file contains a brief header followed by an
array of offsets.  The last offset is to the end of the P-file so that the profile for
sequence b+i is the bytes off[i] to off[i+1]-1 where off[i] is the i'th offset in the A-file.

```
      < kmer size(k)                                     : int   >
      < index of sequence of 1st profile in this file(b) : int64 >
      < # of profile offsets in this file(n)             : int64 >
      ( < profile offset for sequence i in [b,b+n) : int64 > ) ^ n+1
```

The P-file contains compressed profiles.  The sequence of a profile is given by the first count followed by the first forward difference to each successive count as these are expected to be
small integers that will compress well.  The the first count is encoded in the first one or two bytes, depending on
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