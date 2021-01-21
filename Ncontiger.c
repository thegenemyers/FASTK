/*******************************************************************************************
 *
 *  scaff2contig is a pipe that reads in a fasta or fastq file from stdin presumed to contain
 *    an assembly where each sequence is a collection of contigs separated by N's to denote
 *    gaps of a size equal to the estimated length of the gap.  It outputs a fasta file
 *    of the contigs, where the header of each contig is the # of the scaffold, the # of
 *    the contig within the scaffold and then the original header line for the relevant
 *    scaffold.
 *
 *  Author:  Gene Myers
 *  Date  :  January 2021
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <sys/stat.h>
#include <unistd.h>

#include "libfastk.h"

static char *Usage = "< <in:fast[aq]> > <out:fasta>";

typedef struct
  { char *header;  // header line (excluding 1st char)
    char *seq;     // DNA sequence
    int   len;     // length of DNA sequence
  } Entry;

static FILE *INPUT;
static FILE *OUTPUT;

//  Read next line into a buffer and return a pointer to the buffer and set *plen
//    the length of the line.  NB: replaces '\n' with '\0'.

static char *read_line(int *plen)
{ static char *buffer;
  static int   bmax = 0;
  int len;
  
  if (bmax == 0)
    { bmax = 500;
      buffer = (char *) Malloc(bmax,"Allocating line buffer");
      if (buffer == NULL)
        exit (1);
    }
  
  if (fgets(buffer,bmax,INPUT) == NULL)
    return (NULL);
  
  len = strlen(buffer);
  while (buffer[len-1] != '\n') 
    { bmax = ((int) (1.4*bmax)) + 100;
      buffer = (char *) Realloc(buffer,bmax,"Reallocating line buffer");
      if (buffer == NULL)
        exit (1);
      if (fgets(buffer+len,bmax-len,INPUT) == NULL)
        { fprintf(stderr,"%s: Last line of file does not end with new-line\n",Prog_Name);
          exit (1);
        } 
      len += strlen(buffer+len);
    }
  buffer[--len] = '\0';
  
  *plen = len;
  return (buffer);
}

static Entry *get_fastq_entry()
{ static Entry entry;
  static int   hmax = 0;
  static int   smax = 0;
  static char *seq = NULL;
  static char *header = NULL;

  int   hlen, slen, qlen;
  char *line;

  line = read_line(&hlen);
  if (line == NULL)
     return (NULL);

  if (line[0] != '@')
    { fprintf(stderr,"%s: Entry header does not start with an @-sign\n",Prog_Name);
      exit (1);
    }
  if (hlen > hmax)
    { hmax = 1.2*hlen + 1000; 
      header = Realloc(header,hmax+1,"Allocating header buffer");
      if (header == NULL)
        exit (1);
    }
  memcpy(header,line,hlen+1);

  line = read_line(&slen);
  if (line == NULL)
    { fprintf(stderr,"%s: Incomplete fastq entry\n",Prog_Name);
      exit (1);
    }
  if (slen > smax)
    { smax = 1.2*slen + 10000; 
      seq = Realloc(seq,smax+1,"Allocating sequence buffer");
      if (seq == NULL)
        exit (1);
    }
  memcpy(seq,line,slen+1);

  line = read_line(&qlen);
  if (line == NULL)
    { fprintf(stderr,"scaff2contig: Incomplete fastq entry\n");
      exit (1);
    }
  if (line[0] != '+')
    { fprintf(stderr,"scaff2contig: Divider line does not start with a +-sign\n");
      exit (1);
    }

  line = read_line(&qlen);
  if (line == NULL)
    { fprintf(stderr,"scaff2contig: Incomplete fastq entry\n");
      exit (1);
    }
  if (slen != qlen)
    { fprintf(stderr,"scaff2contig: QV line does not have the same length as sequence line\n");
      exit (1);
    }

  entry.header = header;
  entry.seq    = seq;
  entry.len    = slen;
  return (&entry);
}


static Entry *get_fasta_entry()
{ static Entry entry;
  static int   hmax = 0;
  static int   smax = 0;
  static char *seq = NULL;
  static char *header = NULL;
  static char *next = NULL;
  static int   empty = 1;

  int   hlen, slen, m;
  char *line;

  if (empty)
    { empty = 0;
      next = read_line(&hlen);
      if (next == NULL)
        { fprintf(stderr,"%s: INPUT is empty!\n",Prog_Name);
          exit (1);
        }
      if (next[0] != '>')
        { fprintf(stderr,"%s: Entry header does not start with a >-sign\n",Prog_Name);
          exit (1);
        }
    }
  else
    { if (next == NULL)
        return (NULL);
      hlen = strlen(next);
    }

  if (hlen > hmax)
    { hmax = 1.2*hlen + 1000; 
      header = Realloc(header,hmax+1,"Allocating header buffer");
      if (header == NULL)
        exit (1);
    }
  memcpy(header,next,hlen+1);

  slen = 0;
  while (1)
    { line = read_line(&m);
      if (line == NULL || line[0] == '>')
        break;

      if (slen+m > smax)
        { smax = 1.2*(slen+m) + 1000000; 
          seq = Realloc(seq,smax+1,"Allocating sequence buffer");
          if (seq == NULL)
            exit (1);
        }

      memcpy(seq+slen,line,m);
      slen += m;
    }
  if (slen == 0)
    { if (line == NULL)
        fprintf(stderr,"%s: Sequence missing after last header\n",Prog_Name);
      else
        fprintf(stderr,"%s: Sequence missing between two headers\n",Prog_Name);
      exit (1);
    }

  next = line;
  entry.header = header;
  entry.seq    = seq;
  entry.len    = slen;
  return (&entry);
}


int main(int argc, char *argv[])
{
  INPUT   = stdin;
  OUTPUT  = stdout;
  
  { int    i, j, k;
    int    flags[128];

    ARG_INIT("Ncontiger")

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

    if (argc != 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  { Entry *entry, *(*fetch)();
    int    len;
    char  *seq;
    int    s, c;
    int     i, j, k;

    c = fgetc(INPUT);
    if (c == '@')
      fetch = get_fastq_entry;
    else
      fetch = get_fasta_entry;
    ungetc(c,INPUT);

    for (s = 0; (entry = fetch()) != NULL; s++)
      { seq = entry->seq;
        len = entry->len;
        for (k = 0; k < len; k++)
          if (seq[k] != 'N' && seq[k] != 'n')
            break;
        for (c = 0, i = k; i < len; i = k, c++)
          { for (j = i+1; j < len; j++)
              if (seq[j] == 'N' || seq[j] == 'n')
                break;
            if (j < len)
              { for (k = j+1; k < len; k++)
                  if (seq[k] != 'N' && seq[k] != 'n')
                    break;
              }
            else
              k = j;

            fprintf(OUTPUT,"> %d %d %d %d :: %s\n",s,c,i,j,entry->header);
            for ( ; i+100 < j; i += 100)
              fprintf(OUTPUT,"%.100s\n",seq+i);
            fprintf(OUTPUT,"%.*s\n",j-i,seq+i);
          }
      }
  }

  exit (0);
}
