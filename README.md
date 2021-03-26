# FastK: A K-mer counter (for HQ assembly data sets)
  
<font size ="4">**_Author:  Gene Myers_**<br>
**_First:   July 22, 2020_**<br>
**_Current: November 18, 2020_**</font>

- [Command Line](#command-line)
  - [FastK](#fastk)
  - [Fastrm, Fastcp, & Fastmv](#fastrm)

- [Sample Applications](#sample-applications)
  - [Histex](#histex): Display a FastK histogram
  - [Tabex](#tabex): List, Check, or find a k&#8209;mer in a FastK table
  - [Profex](#profex): Display a FastK profile
  - [Logex](#logex): Combine kmer,count tables with logical expressions & filter with count cutoffs
  - [Vennex](#vennex): Produce histograms for the Venn diagram of 2 or more tables
  - [Haplex](#haplex): Find k&#8209;mer pairs with a SNP in the middle
  - [Homex](#homex): Estimate homopolymer error rates

- [C-Library Interface](#c-library-interface)
  - [K-mer Histogram Class](#k-mer-histogram-class)
  - [K-mer Table Class](#k-mer-table-class)
  - [K-mer Stream Class](#k-mer-stream-class)
  - [K-mer Profile Class](#k-mer-profile-class)
 
- [File Encodings](#file-encodings)
  - [`.hist`: K-mer Histogram File](#k-mer-histogram-file)
  - [`.ktab`: K-mer Table Files](#k-mer-table-files)
  - [`.prof`: K-mer Profile Files](#k-mer-profile-files)


## Command Line

FastK is a k&#8209;mer counter that is optimized for processing high quality DNA assembly data
sets such as those produced with an Illumina instrument or a PacBio run in HiFi mode.
For example it is about 2 times faster than KMC3 when counting 40-mers in a 50X HiFi data
set.  Its relative speedup decreases with increasing error rate or increasing values of k,
but regardless is a general program that works for any DNA sequence data set and choice of k.
It is further designed to handle data sets of arbitrarily large size, e.g. a 100X data
set of a 32GB Axolotl genome (3.2Tbp) can be performed on a machine with just 12GB of memory
provided it has ~6.5TB of disk space available.

FastK can produce the following outputs:

1. a histogram of the frequency with which each k&#8209;mer in the data set occurs.
2. a table of k&#8209;mer/count pairs sorted lexicographically on the k&#8209;mer where a < c < g < t.
3. a k&#8209;mer count profile of every sequence in the data set.  A **profile** is the sequence of counts of the n-(k-1) consecutive k&#8209;mers of a sequence of length n.
4. a **relative profile** of every sequence in the data set against a FastK table produced for another data set.

Note carefully, that in order to accommodate the unknown orientation of a sequencing read,
a k&#8209;mer and its Watson Crick complement are considered to be the same k&#8209;mer by FastK, where the
lexicograpahically smaller of the two alternatives is termed **canonical**.
The histogram is always produced whereas the production of a k&#8209;mer table (2.) and profiles (3.&4.)
are controlled by command
line options.  The table (2.) is over just the *canonical* k&#8209;mers present in the data set.  Producing profiles (3.&4.) as part of the underlying sort is much more efficient than producing them after the fact using a table or hash of all k&#8209;mers such as is necessitated when using other k&#8209;mer counter programs.  The profiles are recorded in a space-efficient compressed form, e.g.
about 4.7-bits per base for a recent 50X HiFi asssembly data set.

<a name="fastk"></a>

```
1. FastK [-k<int(40)>] [-t[<int(4)>]] [-p[:<table>[.ktab]]] [-c] [-bc<int>]
          [-v] [-N<path_name>] [-P<dir(/tmp)>] [-M<int(12)>] [-T<int(4)>]
            <source>[.cram|.[bs]am|.db|.dam|.f[ast][aq][.gz]] ...
```

FastK counts the number of k&#8209;mers in a corpus of DNA sequences over the alphabet {a,c,g,t} for a specified k&#8209;mer size, 40 by default.
The input data can be in one or more CRAM, BAM, SAM, fasta, or fastq files, where the later two can be gzip'd. The data can also be in [Dazzler databases](https://github.com/thegenemyers/DAZZ_DB).  The type of the file is determined by its extension (and not its contents).  The extension need not be given if the root name suffices to uniquely identify a file.  If more than one source file is given
they must all be of the same type in the current implementation.

FastK produces a number of outputs depending on the setting of its options.  By default, the
outputs will be placed in the same directory as that of the first input and begin with the
prefix \<source> which is the first path name absent any suffix extensions.  For example, if
the input is <code>../BLUE/foo.fastq</code> then \<source> is <code>../BLUE/foo</code>, the
outputs will be placed in directory <code>../BLUE</code>, and all result file names will begin with <code>foo</code>.  If the &#8209;N option is specified then the path name specified is used
as \<source>.

One can select any value of k &ge; 5 with the &#8209;k option.
FastK always outputs a file with path name `<source>.hist` that contains a histogram of the k&#8209;mer frequency
distribution where the highest possible count is 2<sup>15</sup>-1 = 32,767 -- FastK clips all higher values to this upper limit.  Its exact format is described in the section on Data Encodings.

One can optionally request, by specifying the &#8209;t option, that FastK produce a sorted table of
all canonical k&#8209;mers that occur &#8209;t or more times along with their counts, where the default
for the threshold is 4.
The output is placed in a single *stub* file with path name `<source>.ktab` and N
roughly equal-sized *hidden* files with the path names `<dir>/.<base>.ktab.#` assuming
\<source> = \<dir>/\<base> and
where # is a thread number between 1 and N where N is the number of threads used by FastK (4 by default).
The exact format of the N&#8209;part table is described in the section on Data Encodings.

One can also ask FastK to produce a k&#8209;mer count profile of each sequence in the input data set
by specifying the &#8209;p option.  A single *stub* file with path name `<source>.prof` is output
along with 2N roughly equal-sized pairs of *hidden* files with path names
`<dir>/.<base>.pidx.#` and `<dir>/.<base>.prof.#` in the order of the sequences in the input assuming \<source> = \<dir>/\<base>.
The profiles are individually compressed and the exact format of these
files is described in the section on Data Encodings.

If the data file contains sequences with letters other than upper or lower case a, c, g, or t,
then all k-mers involving these letters are considered *invalid* and they are not counted.
Specifically, the do not occur in the k-mer table and in profiles they are regions of 2k-1
or more 0's.  So for example, if one passes a fasta "assembly" file to FastK wherein gaps
between contigs are indicated by runs of N's, then the profile of a scaffold "sequence" will
contain a corresponding run of 0's where the contig gaps are.

The -p option can contain an optional reference to a k&#8209;mer table such as produced by the
-t option.  If so, then FastK produces profiles of every read where the k&#8209;mer counts
are those found in the referenced table, or zero if a k&#8209;mer in a read is not in the
table.  This *relative* profile is often useful to see how the k&#8209;mers from one source
are reflected in another by tools such as [merfin](https://github.com/arangrhie/merfin).
They could also be used to distinguish haplotypes in a trio-based project, by producing
relative profiles with respect to the k&#8209;mers of the father and mother sequencing data sets.
If this version of the -p option is specified then only profiles are produced -- the
-t option is ignored and the default histogram is not produced.

The &#8209;c option asks FastK to first homopolymer compress the input sequences before analyzing
the k&#8209;mer content.  In a homopolymer compressed sequence, every substring of 2 or more a's
is replaced with a single a, and similarly for runs of c's, g's, and t's.  This is particularly useful for Pacbio data where homopolymer errors are five&#8209;fold more frequent than other
errors and thus the error rate of such "hoco" k&#8209;mers is five&#8209;fold less.

The &#8209;v option asks FastK to output information about its ongoing operation to standard error including a time and resource summary at completion.
The &#8209;bc option allows you to ignore the prefix of each read of the indicated length, e.g. when
the reads have a bar code at the start of each read.
The &#8209;P option specifies where FastK should place all the numerous temporary files it creates, if not `/tmp` by default.
The &#8209;M option specifies the maximum amount of memory, in GB, FastK should use at any given
moment.
FastK by design uses a modest amount of memory, the default 12GB should generally
be more than enough.
Lastly, the &#8209;T option allows the user to specify the number of threads to use.
Generally, this is ideally set to the actual number of physical cores in one's machine.
            
<a name="fastrm"></a>

```
2a. Fastrm [-if] <source>[.hist|.ktab|.prof] ...
2b. Fastmv [-inf] <source>[.hist|.ktab|.prof] ( <target> | ... <directory> )
2c. Fastcp [-inf] <source>[.hist|.ktab|.prof] ( <target> | ... <directory> )
```

As described above FastK produces hidden files whose names begin with a . for the &#8209;t and &#8209;p
options in order to avoid clutter when listing a directory's contents.
An issue with this approach is that it is inconvenient for the user to remove, rename, or copy these files
and often a user will forget the hidden files are there, potentially wasting disk space.
We therefore provide Fastrm, Fastmv, and Fastcp that remove, rename, and copy FastK .hist, .ktab, and .prof output files as a single unit.

If \<source> does not end with a FastK extenion then the command operates on any histogram, k&#8209;mer table, and profile files with \<source> as its prefix.  Otherwise the command operates on the file with the given extension and its hidden files.  Fastrm removes the relevant stubs and
hidden files, Fastmv renames all the relevant files as if FastK had been called with option &#8209;N\<target>, and Fastcp makes a copy of all associated files with the path name \<target>.  For
Fastmv and Fastcp, if the last argument is a directory (as opposed to a path name), then any number of stub files can be moved or copied where
the base name of each source is combined with \<directory> to form complete target path names.

As for the UNIX rm, mv, and cp commands, the &#8209;i option asks the command to ask the user
about each file as to whether you want to delete it (rm) or overwrite an existing file (mv,cp), but only for the stubs and not the hidden files corresponding to each stub, which share the same fate as their stub file.
The &#8209;n option tells Fastmv and Fastcp to not overwrite any existing files.
Finally, the &#8209;f option forces the creation of the new files and overides both the
&#8209;i and &#8209;n options.
            
### Current Limitations & Known Bugs

Currently if multiple input files are given they must all be of the same type, e.g. fasta
or cram.  This restriction is not fundamental and could be removed with some coding effort.

FastK is not working when memory exceeds 128GB.  This should generally not be an issue as it is designed specifically to not require large memory, 16GB should always be enough.  It does operate a bit faster with a lot of memory though, so we will track down the 32-bit
integers(s) that need to be 64-bit.

FastK is not working for k greater than roughly 128.  Again this is an unusually large k for a practical application but in principle it should work for unlimited k and we will address this problem shortly.

&nbsp;

&nbsp;

## Sample Applications

<a name="histex"></a>
```
1. Histex [-kAG] [-h[<int(1)>:]<int(100)>] <source>[.hist]
```

This command and also Tabex and Profex are presented specifically to
give a user simple examples of how to use the C&#8209;interface modules,
<code>libfastk.c</code>,
to manipulate the histogram, k&#8209;mer count tables, and profiles produced by FastK.

Given a histogram file \<source>.hist produced by FastK,
one can view the histogram of k&#8209;mer counts with **Histex** where the &#8209;h specifies the 
interval of frequencies to be displayed where 1 is assumed if the lower bound is not given.

If the &#8209;k option is *not* set then the counts of the histogram are of the # of unique
k-mers U(f) that occur with the given frequency f in the input data set, but if it is set then
the counts are of the # of k-mer instances I(f) that occur with the given frequency.  I(f) = f&middot;U(f)
except for possibly the lowest and highest frequency of the histogram as these entries
include the counts for the frequencies below and above them, respectively.

If the -A option is set, then Histex outputs a simple ASCII listing suitable for
input to other programs.  If the -G option is set, then Histex produces a special ASCII
histogram where the last entry is specifically adjusted so that **GeneScope.FK**
outputs a correct estimate of genome size and % repetitiveness.  (**GenomeScope2.0**
requires an unbounded histogram whereas GeneScope.FK typically only requires the
first 1000 or so frequencies.)

<a name="tabex"></a>
```
2. Tabex [-t<int>] <source>[.ktab] (LIST|CHECK|(<k-mer:string>) ...
```

Given that a set of k&#8209;mer counter table files have been generated represented by stub file
\<source>.ktab, ***Tabex*** opens the corresponding hidden table files (one per thread) and then performs the sequence of actions specified by the
remaining arguments on the command line.  The literal argument LIST lists the contents
of the table in radix order.  CHECK checks that the table is indeed sorted.  Otherwise the
argument is interpreted as a k&#8209;mer and it is looked up in the table and its count returned
if found.  If the &#8209;t option is given than only those k&#8209;mers with counts greater or equal to the given value are operated upon.

<a name="profex"></a>
```
3. Profex <source>[.prof] <read:int> ...
```

Given that a set of profile files have been generated and are represented by stub file
\<source>.prof, ***Profex*** opens the corresonding hidden profile files (two per thread)
and gives a display of each sequence profile whose ordinal id is given on
the remainder of the command line.  The index of the first read is 1 (not 0).

<a name="logex"></a>
```
4. Logex [-T<int(4)>] [-[hH][<int(1)>:]<int>] <name=expr> ... <source>[.ktab] ...
```

*UNDER CONSTRUCTION*

Logex takes one or more k&#8209;mer table "assignments" as its initial arguments and applies these to the ordered merge of the k&#8209;mer count tables that follow, each yielding a new k&#8209;mer tables with the assigned names, of the k&#8209;mers satisfying the logic of the associated expression along with counts computed per the "modulators" of the expression.  For example,
`Logex 'AnB = A &. B' Tab1 Tab2` would produce a new table stub file AnB.ktab
and associated hidden files, of the k&#8209;mers common to the tables represented by the
stub files Tab1.ktab and Tab2.ktab.  If the &#8209;h option is given then a histogram over
the given range is generated for each asssignment, and if the &#8209;H option is given then
only the histograms are generated and not the tables.  The &#8209;T option can be used to
specify the number of threads used.

Each assignment arguments is a path name followed by an =-sign and then a "k&#8209;mer&#8209;count"
expression.  The path name specifies the location and name of the table that will be
produced in response to the application of the k&#8209;mer&#8209;count expresssion to the input
tables.

A k&#8209;mer&#8209;count expression has as its basis a logical predicate made up from the binary
operators '|' (or), '&' (and), '^' (xor), and '-' (minus) over arguments that are alphabetic letters from a-h or A-H where case does not matter.
So for example, the logical predicate `(A^B)-C` would select those k&#8209;mers that occur in either the first or second table, but not both, and that do not occur in the third table.  The order of precedence of the operators is '&' (highest), then '^' then '-' then '|' (lowest).  Parenthesis can be used to override precedence and spaces may be freely interspersed in the expression.
If there are k &le; 8 table arguments after the assignments, then the assignment expressions in toto are expected to involve the k-consecutive letters starting with 'a'.

FastK tables are not just ordered lists of k&#8209;mers, but ordered lists of k&#8209;mers *with a
count for each*, i.e. k&#8209;mer,count pairs.  So a k&#8209;mer&#8209;count expression must also specify
how to combine the counts of the k&#8209;mers that satisfy the logical operations.
This is accomplished by adding two new unary operators, '[]' and '#', and adding a
"modulator" character to the logical operators.  Specifically:

* Any sub-expression can be followed by a post-fix modulation operator consisting of a
comma separated list of integer ranges in square brackets, i.e.
'[' \<range> \( ',' \<range> )\* ']'  where \<range> is an integer range where either
the upper or lower bound (or both) can be missing, i.e. [\<int>] '-' [\<int>], or an
integer, in which case the range is just the specific integer.
Only those k&#8209;mer,count pairs produced by the filter's sub-expression whose counts are
in one of the supplied ranges is accepted by this modulator expression.  For example,
`A[5-10]` accepts all k&#8209;mers in the first table with count between 5 and
10 (inclusive), and `(A-B)[7-]` accepts all k&#8209;mer that are in the first
table but not the second and have a count of 7 or more.  With regard to precedence,
this operator binds more tightly than any of the logical operators.

* Any sub-expression can be preceded by the prefix modulation operator #, which
returns the k&#8209;mer of its subexpression with a count of 1 (if the sub-expression
produces a k&#8209;mer).  For example, `#A |+ #B |+ #C` will produce the union of
all the k&#8209;mers in the tables where the count will be the number of tables the k&#8209;mer
occurred in.  This operator binds the most tightly with regard to precedence.

* When the same k&#8209;mer is in several of the tables and so accepted by a logical
expression, e.g. it is in both the first and second table and the operator is & or |,
then the question arises as to what the count of the accepted k&#8209;mer should be.
At this time we provide 6 "modulators" that immediately follow the logical operator
as follows:

	* '+' takes the sum of the k&#8209;mers
	* '-' subtracts the count of the left-kmer,counter pair from the right.
	* '<' takes the minimum count
	* '>' takes the maximum count
	* '*' takes the average of the two counts
	* '.' takes the count of the left&#8209;kmer whenever it is available, the right otherwise

So for example, A &+ B will produce a k&#8209;mer, count pair when a k&#8209;mer is in both the first and second
tables and give the k&#8209;mer the sum of the counts of the two instances.
A |+ B will produce a k&#8209;mer, count pair when a k&#8209;mer is in one or both of the first
and second tables and give the k&#8209;mer the sum of the counts of the instances available.
The operators ^ and - do not require modulators as only one k&#8209;mer satisfies the operator
which then uses that count as the count for its result.  As a final example,
`(A |> B |> C |> D)[-3]` will output any k&#8209;mer that has a count of 3 or less in the
first four tables along with its smallest count.

In summary, k&#8209;mer&#8209;count expressions permit all the typical filtration and logical combination operators provided in the post&#8209;count framework of most other k&#8209;mer counter software suites.  Some efficiency may be lost due to the interpretive realization
of the expressions but this is hopefully compensated for by the expressiveness of
the concept which unifies most of the desired manipulations in a single program.

<a name="vennex"></a>
```
5. Vennex [-h[<int(1)>:]<int(100)>] <source_1>[.ktab] <source_2>[.ktab] ...
```

*UNDER CONSTRUCTION*

Vennex takes two or more, say k, tables, and produces histograms of the counts for each
region in the k-way Venn diagram.  That is `Vennex Alpha Beta` where
`Alpha` and `Beta` are .ktab's will produce histograms of the
counts of:

* the k&#8209;mers in both Alpha and Beta, i.e. Alpha & Beta, in file `ALPHA.BETA.hist`
* the k&#8209;mers in Alpha but not Beta, i.e. Alpha-Beta, in file `ALPHA.beta.hist`, and
* the k&#8209;mers in Beta but not Alpha, i.e. Beta-Alpha, in file `alpha.BETA.hist`.

Generalizing,
`Vennex A B C`, produces 7 ( = 2<sup>k</sup>-1) histograms with the names, a.b.C, a.B.c, a.B.C.,
A.b.c, A.b.C, A.B.c, and A.B.C where the convention is that a table name is in upper case if it is in, and the name is in
lower case if it is out.  For example, a.B.c is a histogram of the counts of the k&#8209;mers  that are in B but not A and not C, i.e. B-A-C.  The range of the histograms is 1 to 100 (inclusive) by
default but may be specified with the -h option.

It may interest one to observe that the command `Vennex Alpha Beta` is equivalent to the command
`Logex -H100 'ALPHA.BETA=#A&B' 'ALPHA.beta=#A-B' 'alpha.BETA=#B-A' Alpha Beta` further illustrating the flexibility of the Logex command.

<a name="haplex"></a>
```
6. Haplex [-g<int>:<int>] <source>[.ktab]
```

In a scan of \<source> identify all bundles of 2&#8209;4 k&#8209;mers that differ only in their
middle base, i.e. the &lfloor;k/2&rfloor;<sup>th</sup> base.  If the &#8209;g option is
given then only bundles where the count of each k&#8209;mer is in the specified integer
range (inclusively) are reported.  Each bundle is output to the standard output with
each k&#8209;mer followed by its count on a line and a blank line between bundles.  For example,

```
...

cgatcctatcacttctaggaCccccatatgaatatagata 21
cgatcctatcacttctaggaTccccatatgaatatagata 12

cgatcctatctgtgcagattCccagcagcaccaataagaa 7
cgatcctatctgtgcagattTccagcagcaccaataagaa 19

cgatcctcaaccccggtgtgAgggtttgtttggccccgca 17
cgatcctcaaccccggtgtgGgggtttgtttggccccgca 19

cgatcctcacacttattcgaAcgctttttcggtactcgcc 18
cgatcctcacacttattcgaCcgctttttcggtactcgcc 21
cgatcctcacacttattcgaGcgctttttcggtactcgcc 8

cgatcctcacactttttcgaCgctttttcggtactcgccc 30
cgatcctcacactttttcgaTgctttttcggtactcgccc 26

...
```
in response to <code>Haplex -h7:36 CB.ktab</code> where CB is a 50X HiFi data set of
Cabernet Sauvignon.

<a name="homex"></a>
```
7. Homex -e<int> -g<int>:<int> <source_root>[.ktab]
```
In a scan of \<source> identify all k&#8209;mers that contain a homopolymer straddling the
mid-point with count in the range given by the &#8209;g parameter.  Consider the k&#8209;mers with
one extra or one less homopolymer base aded to the right end of the homopolymer.  If
those k&#8209;mers have count no more than &#8209;e they they are considered homopolymer errors.
The total # of correct and erroneous homopolymer instances for each puridine and pyrimidine and of each length up to 10 is collected and reported.  For example,
<code>Homex -e6 -g10:60 CB.ktab</code> where CB is a 50X HiFi data set of
Cabernet Sauvignon, the following table is output (in about a minute).

```
              -1      Good          +1      Error Rate
  2 at:    1083037 1444004873    1354886 -> 0.2%
  3 at:    2280582  544204150    1342588 -> 0.7%
  4 at:    2967419  231095068    1234288 -> 1.8%
  5 at:    2376973   87922494    1025822 -> 3.7%
  6 at:    1290169   29863257     587888 -> 5.9%
  7 at:     704839   11839930     368393 -> 8.3%
  8 at:     411238    5250289     217606 -> 10.7%
  9 at:     244298    2493850     129673 -> 13.0%
 10 at:     189475    1603996      98766 -> 15.2%

  2 cg:    1023595  764923391     626469 -> 0.2%
  3 cg:    1092236  168190300     307076 -> 0.8%
  4 cg:     608471   33793608     115702 -> 2.1%
  5 cg:     257885    7252509      45142 -> 4.0%
  6 cg:      84391    1570406      16852 -> 6.1%
  7 cg:      22940     351055       6788 -> 7.8%
  8 cg:       3655      49731       1563 -> 9.5%
  9 cg:       1369      15226        719 -> 12.1%
 10 cg:       1167      11200        751 -> 14.6%
```

&nbsp;

&nbsp;


## C-Library Interface

For each of the 3 distinct outputs of FastK, we have suppled a simple C library
that gives a user access to the data therein.  The library is simply embodied in
the C&#8209;file, `libfastk.c`, and associated include file `libfastk.h`.
The makefile commands for building Histex, Tabex, and Profex illustrate how to
easily incorporate the library into your C or C++ code.


### K-mer Histogram Class

A Histogram object is a record with 5 fields
as described in the comments of the declaration below:

```
typedef struct
  { int    kmer;    //  Histogram is for k-mers of this length
    int    unique;  //  1 => count  of unique k-mers, 0 => count of k-mer instances
    int    low;     //  Histogram is for range [low,hgh]
    int    high;
    int64 *hist;    //  hist[i] for i in [low,high] = # of k-mers occuring i times
  } Histogram;
```

The frequencies are stored in the array pointed at by the field `hist`
where indexing said with any value between `low` and `high`,
inclusive will deliver a valid count.  But caution: indexing with any frequency outside this range may result in an out-of-bounds memory access and possible segfault.
If unique is set then the counts are of the number of **unique** k&#8209;mers U(f) that occur with
a given frequency f, and if not set then the counts are the number of k&#8209;mer **instances** I(f) with the given frequency.  Note that I(f) = f&middot;U(f) in general.
By convention, the lowest and highest frequencies always contain the number of k&#8209;mers with the given frequency **plus** the number of k&#8209;mers with lower or higher frequencies, respectively.  This is to ensure the convention that the total sum of the counts in a histogram is equal to the total number of k&#8209;mer instances in the originating source sequence data set (if unique = 0) or the total number of unique k&#8209;mers in the data set (if unique = 1).
This boundary convention does imply that I(low) &ne; low&middot;U(low) and I(high) &ne; high&middot;U(high) but otherwise I(f) = f&middot;U(f) for all f in (low,high).

```
Histogram *Load_Histogram(char *name);
void       Modify_Histogram(Histogram *H, int low, int high, int unique);
int        Write_Histogram(char *name, Histogram *H);
void       Free_Histogram(Histogram *H);
```
`Load_Histogram` opens the FastK histogram at path name `name`
adding the .hist extension if it is not present.  It returns a pointer to a
newly allocated `Histogram` object for the data encoded in the specified
file, where the histogram is of unique k&#8209;mer counts (use `Modify_Histogram` to
change this).  The routine returns NULL if it cannot open the file, and if
there is insufficient memory available (very unlikely given its size), it prints a message
to standard error and exits.

`Modify_Histogram` modifies a given histogram so it is over the given
range and so that its counts reflect the boolean `unique`, i.e. unique or instance
counts.  The routine does nothing if the supplied subrange is not a subset of the range
of the supplied histogram.  The lowest and highest frequencies have the cumulative
counts of the frequencies below and above them, per our convention.

`Write_Histogram` writes the given histogram to the file with path `name` as a FastK
histogram file, adding a .hist extension if not given.  It returns a non-zero value
if it cannot create and write the named file.

`Free_Histogram` removes all memory encoding the input histogram.

&nbsp;

### K-mer Table Class

The Kmer\_Table class offers a simple, basic interface to access FastK table objects.
It loads entire tables into memory and shields the user from the encoding details,
so optimization for say serial access are not possible.  Some parts of a Kmer\_Table
record are visible for convenience but should never be modified by a user:

```
typedef struct
  { int     kmer;       //  Kmer length
    int     minval;     //  The minimum count of a k-mer in the table
    int64   nels;       //  # of unique, sorted k-mers in the table

    void   *private[7]; //  Private fields
  } Kmer_Table;
```

The table is conceptually an array of <code>nels</code> entries 
where each entry encodes a k&#8209;mer, count pair that
are **sorted** in lexicographical order of the k&#8209;mers. 

```
Kmer_Table *Load_Kmer_Table(char *name, int cut_off);
void        Free_Kmer_Table(Kmer_Table *T);

char       *Fetch_Kmer(Kmer_Table *T, int64 i, char *seq);
int         Fetch_Count(Kmer_Table *T, int64 i);

int64       Find_Kmer(Kmer_Table *T, char *seq);
```

`Load_Kmer_Table` opens the FastK k&#8209;mer table represented by the stub file
at path name `name`, adding the .ktab extension if it is not present.  It returns a pointer to a newly allocated `Kmer_Table` object for
the data encoded in the relevant files.  The routine returns NULL if it cannot open the stub file.  If there is insufficient memory available or the hidden files are inconsistent with
the stub file, it prints an informative message to standard error and exits.  This routine attempts to load the entire table into memory and so may fail as these tables can be very large.  For example, if FastK is run on a human genome data set with -t4, the table can require as much as 40-50GB.  As
a result, in the cases where one wants the reduced table of only those k&#8209;mers
whose counts are not less than `cut_off`, then the load actually reads the table
twice with a `Kmer_Stream` to use only the memory required for exactly those
k&#8209;mers.  This can save significant space at the expense of taking more time to load.
`Free_Kmer_Table` removes all memory encoding the table object.

The two `Fetch` routines return the k&#8209;mer and count, respectively, of the
`i`<sup>th</sup> entry in the given table.  `Fetch_Kmer` in particular returns a pointer to an ascii, 0-terminated string giving the k&#8209;mer in lower-case
a, c, g, t.  If the parameter `seq` is not NULL then the string is placed there and
the pointer returned is to `seq` which much be of length at least `kmer+3`.  If `seq`
is NULL then an array of the appropriate size is allocated and returned containing the requested string. 

`Find_Kmer` searches the table for the supplied k&#8209;mer string and returns the
index of the k&#8209;mer if found, or -1 if not found.  The string `seq` must be
at least `kmer` bases long, and if longer, the trailing bases are ignored.  The string
may use either upper- or lower-case Ascii letters.  The input k&#8209;mer need not be
canonical, `Find_Kmer` will automatically search for the canonical form.

The sample code below opens a table for "foo.ktab", prints out the contents of the table, and ends by freeing all memory involved.

```
Kmer_Table *T = Open_Kmer_Table("foo",0);
char       *s = Fetch_Kmer(T,0,NULL);       //  Create buffer s, value ignored
for (int i = 0; i < T->nels; i++)
  printf("%s : %d\n",Fetch_Kmer(T,i,s),Fetch_Count(T,i));
free(s);
Free_Kmer_Table(T);
```

&nbsp;

### K-mer Stream Class

The Kmer\_Stream class realizes a more complex interface to FastK tables that
entails saving memory by buffering and direct access to the raw encoding
if desired (see [K-mer Table Files](#k-mer-table-files) below).

K&#8209;mer tables can be truly large so that when loaded in memory 10's of gigabytes of main memory are required.  On the other hand most operations can be arranged as a scan of
one or more tables especially given that they are sorted, e.g finding the k&#8209;mers common to two tables.  The Kmer\_Stream class is designed for efficient scanning of the
table, trading off speed of random accesss for efficient memory utilization.
As such, it has the concept of a current position and a small, several KB, buffer to
efficiently move sequentially through consecutive positions.  It is possible to jump
to a specific position, but less efficiently as then one must externally seek the table on disk and reload the buffer.

All fields that might be useful are visible to the user as
follows, with the provisio that they are read-only:

```
typedef struct
  { int     kmer;        //  Kmer length
    int     minval;      //  The minimum count of a k-mer in the stream
    int64   nels;        //  # of unique, sorted k-mers in the stream
                     //  Current position
    int64   cidx;        //  Index of current entry (in table as a whole)
    uint8  *csuf;        //  current entry suffix
    int     cpre;        //  current entry prefix
                     // Other useful parameters
    int     ibyte;       //  # of bytes in prefix
    int     kbyte;       //  Kmer encoding in bytes (= ceiling(kmer/4))
    int     tbyte;       //  Kmer+count entry in bytes (= kbyte + 2)
    int     hbyte;       //  Kmer suffix in bytes (= kbyte - ibyte)
    int     pbyte;       //  Kmer,count suffix in bytes (= tbyte - ibyte)
    
    void   *private[9];  //  Private fields
  } Kmer_Stream;
```
A Kmer\_Stream has a current position that is initialized to the first entry in the
table and that is typically then advanced sequentially through the table.
The current position is directly available
in the fields (1) `cidx`, the ordinal index in the table of the current entry,
(2) `cpre`, the first `ibyte` bytes of the bit compressed k-mer encoded as an integer, and (3) `csuf`, a pointer to the remaining `pbyte` bytes of the k-mer,count encoding.  When the current position is at the end of the table `cidx` will equal `nels` and `csuf` will 
be NULL.  The operators for manipulating a table are as as follows:

```
Kmer_Stream *Open_Kmer_Stream(char *name);
Kmer_Stream *Clone_Kmer_Stream(Kmer_Stream *S);
void         Free_Kmer_Stream(Kmer_Stream *S);

void         First_Kmer_Entry(Kmer_Stream *S);
void         Next_Kmer_Entry(Kmer_Stream *S);

char        *Current_Kmer(Kmer_Streaam *S, char *seq);
int          Current_Count(Kmer_Streaam *S);
uint8       *Current_Entry(Kmer_Streaam *S, uint8 *entry);

void         GoTo_Kmer_Index(Kmer_Stream *S, int64 i);
int          GoTo_Kmer_String(Kmer_Stream *S, char *seq);
int          GoTo_Kmer_Entry(Kmer_Stream *S, uint8 *entry);
```

`Open_Kmer_Stream` opens a k&#8209;mer table as a streamable table object.  Note carefully that the routine conceptually **opens** the table for reading, but does not **load** it (into memory).  The routine returns NULL if it cannot open the stub file.  If there is insufficient memory available or the hidden files are inconsistent with the stub file, it prints an informative message to standard error and exits.  The current position or cursor is set to be the start of the table.
`Free_Kmer_Stream` removes all memory encoding the stream object and closes any open
files associated with it.

`Clone_Kmer_Stream` creates a stream object that shares its read-only indexing tables with the input
stream `S`.  This provides space efficiency when opening a table with multiple threads.  One must
take care to free all clones, prior to freeing the stream the clones were spawned from.

`First_Kmer_Entry` sets the position/entry for the stream to the first entry of
the table and `Next_Kmer_Entry` advance the current position to the next entry.
One needs to check if `csuf` is NULL to determine if the position has advanced to the
end of the table.

One can work with the k-mer,count pair for the current entry directly through the fields
`cpre` and `csuf` permiting optimizations in some situations, but for the most part
the three `Current` routines will produce components of an entry in a more convenient
form.
 
`Current_Kmer` returns a pointer to an ascii, 0-terminated string giving the k&#8209;mer at the current position
in lower-case a, c, g, t.  If the parameter `seq` is not NULL then the string is placed there and the pointer returned is to `seq` which much be of length at least `kmer+3`.
If `seq` is NULL then an array of the appropriate size is allocated and returned containing the requested string.  `Current_Count` returns the count of the kmer,count
pair at the current position.

`Current_Entry` has the same calling conventions as `Current_Kmer`, but returns the `tbyte` bit-compressed encoding of a k-mer,count pair.
The k&#8209;mer is encoded in the first `kbyte` = (kmer+3)/4 bytes where each base is compressed into 2&#8209;bits so that each byte contains up to four bases, in order of high bits to low bits.  The
bases a,c,g,t are assigned to the values 0,1,2,3, respectively.  As an example, 0xc6 encodes
tacg.  The last byte is partially filled if `kmer` is not a multiple of 4, and the remainder is guaranteed to be zeroed.  The byte sequence for the k&#8209;mer is then followed by a 2-byte 
unsigned integer count (implying tbytes = kbytes+2) with a maximum value of 32,767 and
a minimum value of `minval`.  So the pointer `entry` if non-NULL should point at an array of at least
ceiling(k/4)+2 (= `tbyte`) bytes.

The three `GoTo` routines allow one to jump to a specific position.
`GoTo_Kmer_Index` sets the current position to the `i`<sup>th</sup> entry of the
table.  `GoTo_Kmer_String` sets the position to the first entry in the table whose
k&#8209;mer is not less the `seq` in canonical form.
`GoTo_Kmer_Entry` sets the position to the first entry in the table whose value is not less than that of the 2-bit compressed k&#8209;mer encoding pointed at by `entry`.
`GoTo_Kmer_String` searches for the k&#8209;mer in canonical form, whereas
`GoTo_Kmer_Entry` searches for the 2-bit compressed k&#8209;mer as given.
Both routines return 1 if the k-mer at the new position matches the search argument, and 0
otherwise.
The GoTo-routines are not efficient, especially `GoTo_Kmer_String` and `GoTo_Kmer_Entry`
which must binary search for the desired position.  They are intended for the expert who
wishes to use them for tasks like partitioning a table for simultaneous processing by multiple threads.

As an example, the code below opens a stream for "foo.ktab", prints out the contents of the table, and ends by freeing all memory involved.

```
Kmer_Stream *S = Open_Kmer_Stream("foo",1);
char        *s = Current_Kmer(S,NULL);       //  Create buffer s, value ignored
for (First_Kmer_Entry(S); S->csuf != NULL; Next_Kmer_Entry(S))
  printf("%s : %d\n",Current_Kmer(S,s),Current_Count(S));
free(s);
Free_Kmer_Stream(S);
```

&nbsp;

### K-mer Profile Class

A Profile\_Index object is a record with 6 fields as described in the comments of the declaration below:

```
typedef struct
  { int    kmer;     //  Kmer length
    int    nparts;   //  # of threads/parts for the profiles
    int    nreads;   //  total # of reads in data set
    int64 *nbase;    //  nbase[i] for i in [0,nparts) = id of last read in part i + 1
    int   *nfile;    //  nfile[i] for i in [0,nparts) = stream for ".prof" file of part i
    int64 *index;    //  index[i] for i in [0,nreads] = offset in relevant part file of
                     //      compressed profile for read i.
  } Profile_Index;
```

Like the k&#8209;mer stream class, the set of all profiles is not **loaded** into memory,
but rather only **opened** so that individual profiles for a sequence
can be read in and uncompressed on demand.  So `nparts`> indicates how many
hidden part files constitute the set of all profiles and `>nfile<` is an open
file stream to each of the npart hidden .prof files containing the compressed profiles.
On the otherhand, all the .pidx files are loaded into an array of <code>nreads+1</code> offsets into the hidden files pointed at by `index` where the small nparts element table `nbase<` is used to resolve which part file a read is in.
Specificaly, the reads whose compressed profile are found in part p, are those in [x,nbase[p]] where x is 0 if p = 0 and nbase[p-1] otherwise.
For those reads whose compressed profile is in part p, the profile is at [y,index[i])
in the stream nfile[p] where y is 0 if i = x and index[i-1] otherwise.

```
Profile_Index *Open_Profiles(char *name);

void Free_Profiles(Profile_Index *P);

int Fetch_Profile(Profile_Index *P, int64 id, int plen, uint16 *profile);
```

`Open_Profiles` opens the FastK profile files represented by the stub file
at path name `name`, adding the .prof extension if it is not present.  It returns a pointer to a newly allocated `Profile_Index` object that facilitates access
to individual compressed profiles in the hidden .prof data files.
The routine returns NULL if it cannot open the stub file.  If there is insufficient memory available or the hidden files are inconsistent with
the stub file, it prints an informative message to standard error and exits.

`Free_Profiles` removes all memory encoding the profile index object.

`Fetch_Profile` uses the index `P` to fetch the compressed profile
for the sequence with ordinal index `id` in the data set (starting with 0),
and attempts to decompress it into the 2-byte integer array presumed to be pointed at
by `profile` of presumed length `plen`.  It always returns the
length of the indicated profile.  If this is greater than plen, then only the first
plen values of the profile are placed in profile, otherwise the entire profile of the
given length is placed at the start of the array.

&nbsp;

&nbsp;

## File Encodings

### K-mer Histogram File

The histogram file has a name of the form `<source>.hist` where \<source> is the
output path name used by FastK.  It contains an initial integer
giving the k&#8209;mer size <k>, followed by two integers giving the range \[l,h] (inclusive) of frequencies in the histogram, followed by two integers giving the count of  k&#8209;mer instances at the boundaries of the histogram, followed by (h-l)+1 64-bit unique k&#8209;mer counts for frequencies l, l+1, ..., h.
Formally,

```
    < kmer size(k)        : int >
    < 1st freq.(l)        : int >
    < last freq.(h)       : int >
    < 1st instance count  : int64 >
    < last instance count : int64 >
    ( <unique k-mer counts for f=l, l+1, ... h : int64 > ) ^ (h-l)+1
```

FastK always outputs a histogram with l=1 and h=32,637, so the file is exactly 262,164
bytes in size.  Other auxiliary programs can produce histograms over a subrange of this range.
Note carefully that the count for the entry h, is actually the count of all k&#8209;mers that occur
h-or-more times, and when l>1, the entry l, is the count of all k&#8209;mers that occur l-or-fewer
times.  This guaranties that the sum of all counts is the number of k&#8209;mers in the input.

Note carefully that the FastK histogram counts the # of **unique** k&#8209;mers that have a given frequency.  So for example if there is a particular k&#8209;mer that occurs exactly 40 times in the data set, then only 1 count is accumulated at frequency 40.  But one may
also want to count the number of k&#8209;mers instances I(f), i.e. add a count for each instance of the k&#8209;mer or 40 counts in our example,
rather than the number of unique k&#8209;mers U(f).  Observe that I(f) = f&middot;U(f) for any given frequency.  However this is not true for f = l or h as these counts are the sum of the counts above and below l and h, respectively.  Therefore, to be able to produce either
type of histogram counts, the encoding contains explicitly I(l) and I(h) as 4th and 5th integers in its header, while the histogram gives the U-counts from l to h, inclusive.

&nbsp;

### K-mer Table Files

A table of canonical k&#8209;mers and their counts is produced in N hidden parts, where N is the number of threads FastK was run with.  These hidden files are identified by
a single *stub* file `<source>.ktab` where \<source> is the output path name used by
FastK.
The information in the stub file is as follows:

```
   < kmer size(k)    : int >
   < # of parts(N)   : int >
   < min count(m)    : int >
   < prefix bytes(p) : int >
   < 1st index to entries with prefix = i+1 : int64 >, i = 0, ... 4^(4p)-1
```
The first 4 integers of the stub file give (1) the k&#8209;mer length, (2) the number of threads FastK was run with, (3) the frequency cutoff (&#8209;t option) used to prune the table, and (4) the number of prefix bytes of each k-mer (when encoded as a 2-bit
compressed byte array) that are indexed by the 4<sup>4p</sup>+1 table, call it IDX, that constitutes the remainder of the stub file.  The i<sup>th</sup> element, IDX[i],
gives the ordinal index of the first element in the sorted table for which its first
4p bases have the value i+1.  Thus the entries in the table whose first 4p bases have
value i can be found in the interval [&nbsp;IDX[i-1],&nbsp;IDX[i]&nbsp;) assuming IDX[-1] = 0.
Note carefully that these intervals are guaranteed not to span table parts.
   
The table file parts are in N hidden files in the same directory as the stub
file with the names `.<base>.ktab.[1,N]` assuming that \<source> = \<dir>/\<base>.
The k&#8209;mers in each part are lexicographically ordered and the k&#8209;mers in part i are all less than the k&#8209;mers in part i+1, i.e. the concatention of the N files in order of thread index is sorted.  The information in each table file is as follows:

```
    < kmer size(k)   : int   >
    < # of k-mers(n) : int64 >
    ( < bit-encoded k-mer : uint8 ^ (ceiling(k/4)-p) > < count : uint16 ) ^ n
```
    
In words, an initial integer gives the kmer size that FastK was run with followed by
a 64&#8209;bit int, n, that gives the # of entries in this part file.  The remainder of the file is then n
k&#8209;mer,cnt pairs save that the first p bytes of the k&#8209;mer encodings are excised
as they can be obtained from the index in the stub file above.  Doing so saves a great
deal of disk space at the expense of a bit more encoding complexity.
The k&#8209;mer bases are encoded in 2-bits with 4 to a byte, from high bits to low bits,
where the bases a,c,g,t are assigned to the values 0,1,2,3, respectively.
For example, 0xc6 encodes tacg.
The last byte is partially filled if k is not a multiple of 4, and the remainder is
guaranteed to be zeroed.  The truncated ceiling(k/4)-p &ge; 0 byte encoding for a k&#8209;mer is then followed by a 2-byte 
unsigned integer count with a maximum value of 32,767.

&nbsp;

### K-mer Profile Files

The read profiles are stored in N pairs of file, an index and a data pair, that are hidden
and identified by a single *stub* file `<source>.prof`.
This stub file contains just the k&#8209;mer length followed by the number of threads FastK was
run with as two integers.
The hidden data files, `.\<base>.prof.[1,N]`, contain the compressed profiles for
each read
in their order in the input data set, and the hidden index files,
`.\<base>.pidx.[1,N]`, 
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
  0x   => x in [0,127]
  1x,y => x.y in [128,32767]
```
That is, if the high order bit is set then the count is the unsigned integer encoded in the
remaining 15-bits of the firsts two bytes.  If it is not set, then it is the
unsigned integer encoded in the remaining 7-bits of the first byte.  The one byte encoding is
used whenever possible.

All subsequent counts are given as first forward differences to the immediately preceeding
count, and are encoded in one or two bytes.  Often there is a run of 0 differences in a
profile and a special one byte form codes this as a run length.  Formally the forward differences are encoded as follows:

```
  00x  =>  x in [1,63] 0 diffs
  010x =>  x  in [1,31]
  011x => -x in [-1,-31]
  1x,y =>  x.y in 15 bit 2's complement in [-16384,-32] U [32,16383]
                                             (modulo 32768 arithmetic)
```

If the high order bit of the current byte is set, then the remaining 15-bits of the current
byte and the one following it are interpreted as a *signed* 2's complement integer.
Moreover, the forward difference is taken module 32768 (2^15), e.g. 32766 + 3 = 1.
If the 2 highest order bits of the current byte are zero, then the next x+1 counts are the
same as the most recent count, where x is the lower 6-bits interpreted as an unsigned integer.
If the 2 highest order bit of the current byte are 01, then the remaining 6 bits are interpreted
as a *1's complement* integer and the difference is one more or less than said value depending
on the sign.  The single byte encoding is used whenever possible.