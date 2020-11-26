/*******************************************************************************************
 *
 *  Logical Expression Parser & Evaluator
 *
 *  Author:  Gene Myers
 *  Date  :  Oct. 31, 2016
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#define PRINT_TREE

#include "libfastk.h"

static char *Usage[] = { " [-T<int(4)>] [-[hH]{<int(1)>:]<int>]",
                         "   <output:name=expr> ... <source_root>[.ktab] ..." };


/****************************************************************************************
 *
 *  Logical expression parser
 *
 *****************************************************************************************/

#define OP_OR  0
#define OP_AND 1
#define OP_MIN 2
#define OP_XOR 3
#define OP_CNT 4
#define OP_ARG 5

#define MOD_AVE 0
#define MOD_SUM 1
#define MOD_SUB 2
#define MOD_MIN 3
#define MOD_MAX 4
#define MOD_LFT 5

#ifdef PRINT_TREE

static char *Operator[] =
  { "OR", "AND", "MIN", "XOR", "CNT", "ARG" };

static char *Modulator[] =
  { "AVE", "SUM", "SUB", "MIN", "MAX", "LFT" };

#endif

typedef struct _node
  { uint16        op;
    uint16        mode;
    struct _node *lft;
    struct _node *rgt;
  } Node;

static char *Scan;
static int   Error;
static int   VarV;
static Node *or();

#define ERROR(msg)	\
{ Error = msg;		\
  return (NULL);	\
}

static char *Error_Messages[] =
  { "Out of memory",
    "Expecting an argument or (",
    "Expecting )",
    "Expecting a count modulator",
    "Premature end of expression",
    "Expecting a -",
    "Expecting a , or ]",
    "Do not recognize this operator"
    "Largest argument possible is H/h"
  };

static Node *node(int op, int mode, Node *lft, Node *rgt)
{ Node *v;

  v = (Node *) malloc(sizeof(Node));
  if (v == NULL)
    ERROR(0)
  v->op   = op;
  v->mode = mode;
  v->lft  = lft;
  v->rgt  = rgt;
  return (v);
}

static void free_tree(Node *v)
{ if (v->op < OP_ARG)
    { if (v->lft != NULL)
        free_tree(v->lft);
      if (v->op == OP_CNT)
        free(v->rgt);
      else if (v->rgt != NULL)
        free_tree(v->rgt);
    }
  free(v);
}

static Node *terminal()
{ while (isspace(*Scan))
    Scan += 1;

  if (*Scan == '(')
    { Node *v;

      Scan += 1;
      v = or();
      if (v == NULL)
        return (NULL);
      while (isspace(*Scan))
        Scan += 1;
      if (*Scan != ')')
        ERROR(2)
      Scan += 1;
      return (v);
    }
  else if (isalpha(*Scan))
    { int64 x;

      if (islower(*Scan))
        x = *Scan-'a';
      else
        x = *Scan-'A';
      if (x > 8)
        ERROR(8)
      VarV |= (1 << x);
      Scan += 1;
      return (node(OP_ARG,0,(Node *) x,NULL));
    }
  else
    if (*Scan == '\0')
      ERROR(4)
    else
      ERROR(1)
}

static int get_number()
{ int x;

  x = *Scan++ - '0';
  while (isdigit(*Scan))
    x = 10*x + (*Scan++ - '0');
  while (isspace(*Scan))
    Scan += 1;
  return (x);
}

static Node *filter()
{ Node *v;

  v = terminal();
  if (v == NULL)
    return (NULL);

  while (isspace(*Scan))
    Scan += 1;
  while (*Scan == '[')
    { int   len, *f;
      char *beg;

      Scan += 1;
      while (isspace(*Scan))
        Scan += 1;

      beg = Scan;
      len = 0;
      while (1)
        { if (isdigit(*Scan))
            get_number();
          if (*Scan != '-')
            { free_tree(v);
              ERROR(5)
            }
          Scan += 1;
          while (isspace(*Scan))
            Scan += 1;
          if (isdigit(*Scan))
            get_number();
          len += 2;
          if (*Scan == ']')
            break;
          if (*Scan != ',')
            { free_tree(v);
              ERROR(6)
            }
          Scan += 1;
          while (isspace(*Scan))
            Scan += 1;
        }

      f = malloc(sizeof(int)*len);
      if (f == NULL)
        { free_tree(v);
          ERROR(0)
        }

      Scan = beg;
      len  = 0;
      while (1)
        { if (isdigit(*Scan))
            f[len++] = get_number();
          else
            f[len++] = 1;
          Scan += 1;
          while (isspace(*Scan))
            Scan += 1;
          if (isdigit(*Scan))
            f[len++] = get_number();
          else
            f[len++] = 0x7fff;
          if (*Scan == ']')
            break;
          Scan += 1;
          while (isspace(*Scan))
            Scan += 1;
        }
      Scan += 1;

      v = node(OP_CNT,len,v,(Node *) f);
    } 
  return (v);
}

static int get_mode()
{ switch (*Scan)
  { case '*':
      return (MOD_AVE);
    case '+':
      return (MOD_SUM);
    case '-':
      return (MOD_SUB);
    case '<':
      return (MOD_MIN);
    case '>':
      return (MOD_MAX);
    case '.':
      return (MOD_LFT);
    default:
      return (-1);
  }
}

static Node *and()
{ Node *v, *w; 
  int   m;

  v = filter();
  if (v == NULL)
    return (NULL);
  while (1)
    { while (isspace(*Scan))
        Scan += 1;
      if (*Scan != '&')
        return (v);
      Scan += 1;
      m = get_mode();
      if (m < 0)
        { free_tree(v);
          ERROR(3)
        }
      Scan += 1;
      w = filter();
      if (w == NULL)
        { free_tree(v);
          return (NULL);
        }
      v = node(OP_AND,m,v,w);
    }
}

static Node *xor()
{ Node *v, *w; 

  v = and();
  if (v == NULL)
    return (NULL);
  while (1)
    { while (isspace(*Scan))
        Scan += 1;
      if (*Scan != '^')
        return (v);
      Scan += 1;
      w = and();
      if (w == NULL)
        { free_tree(v);
          return (NULL);
        }
      v = node(OP_XOR,0,v,w);
    }
}

static Node *minus()
{ Node *v, *w; 

  v = xor();
  if (v == NULL)
    return (NULL);
  while (1)
    { while (isspace(*Scan))
        Scan += 1;
      if (*Scan != '-')
        return (v);
      Scan += 1;
      w = xor();
      if (w == NULL)
        { free_tree(v);
          return (NULL);
        }
      v = node(OP_MIN,0,v,w);
    }
}

static Node *or()
{ Node *v, *w; 
  int   m;

  v = minus();
  if (v == NULL)
    return (NULL);
  while (1)
    { while (isspace(*Scan))
        Scan += 1;
      if (*Scan == '\0' || *Scan == ')')
        return (v);
      if (*Scan != '|')
        return (v);
      Scan += 1;
      m = get_mode();
      if (m < 0)
        { free_tree(v);
          ERROR(3)
        }
      Scan += 1;
      w = minus();
      if (w == NULL)
        { free_tree(v);
          return (NULL);
        }
      v = node(OP_OR,m,v,w);
    }
}

#ifdef PRINT_TREE

static void print_tree(Node *v, int level)
{ if (v->op == OP_CNT)
    { int *f, i;

      printf("%*s%s",level,"",Operator[v->op]);
      f = (int *) (v->rgt);
      for (i = 0; i < v->mode; i+= 2)
        printf(" %d-%d",f[i],f[i+1]);
      printf("\n");
      print_tree(v->lft,level+2);
    }
  else if (v->op == OP_ARG)
    printf("%*s%s %lld\n",level,"",Operator[v->op],(int64) (v->lft));
  else
    { printf("%*s%s",level,"",Operator[v->op]);
      if (v->op <= OP_AND)
        printf(" %s",Modulator[v->mode]);
      printf("\n");
      print_tree(v->lft,level+2);
      print_tree(v->rgt,level+2);
    }
}

#endif

static Node *parse_expression(char *expr, int *pargs)
{ Node *v;
  int   i, nargs;

  VarV = 0;
  Scan = expr;
  v    = or();
  if (v != NULL && *Scan != '\0')
    { free_tree(v);
      Error = 7;
      v = NULL;
    }

  if (v == NULL)
    { if (Error == 0)
        fprintf(stderr,"%s: Out of memory parsing expression\n",Prog_Name);
      else
        { fprintf(stderr,"%s: Expression syntax error (%d):\n\n",Prog_Name,Error);
          fprintf(stderr,"    %s\n",expr);
          fprintf(stderr,"%*s^ %s\n",(int) ((Scan-expr)+4),"",Error_Messages[Error]);
        }
      return (NULL);
    }

  nargs = 0;
  for (i = 0; i < 8; i++)
    if ((VarV & (1 << i)) != 0)
      { if (nargs < i)
          { fprintf(stderr,"%s: Expression refers to '%c' but not '%c'?\n",
                           Prog_Name,'a'+i,'a'+(i-1));
            return (NULL);
          }
        nargs = i+1;
      }

  *pargs = nargs;
  return (v);
}

static int eval_logic(Node *t, int i)
{ switch (t->op)
  { case OP_OR:
      return (eval_logic(t->lft,i) || eval_logic(t->rgt,i));
    case OP_MIN:
      return (eval_logic(t->lft,i) && (1-eval_logic(t->rgt,i)));
    case OP_XOR:
      return (eval_logic(t->lft,i) != eval_logic(t->rgt,i));
    case OP_AND:
      return (eval_logic(t->lft,i) && eval_logic(t->rgt,i));
    case OP_CNT:
      return (eval_logic(t->lft,i));
    case OP_ARG:
      { int64 x = (int64) (t->lft);
        if ( (i & (1<<x)) != 0)
          return (1);
        else
          return (0);
      }
    default:
      return (0);
  }
}

static int *compile_expression(Node *t, int ntabs)
{ int *table;
  int  i;

  table = Malloc(sizeof(int)*(1<<ntabs),"Allocating logic table");
  for (i = (1 << ntabs)-1; i >= 0; i--)
    table[i] = eval_logic(t,i);
  return (table);
}

static int eval_expression(Node *t, int *cnts)
{ switch (t->op)
  { case OP_OR:
    case OP_AND:
      { int x, y;

        x = eval_expression(t->lft,cnts);
        y = eval_expression(t->rgt,cnts);
        switch (t->mode)
        { case '*':
            return ((x+y) >> 2);
          case '+':
            return (x+y);
          case '-':
            x -= y;
            if (x < 0)
              x = 0;
            return (x);
          case '<':
            if (x < y)
              return (x);
            else
              return (y);
          case '>':
            if (x > y)
              return (x);
            else
              return (y);
          case '.':
            if (x == 0)
              return (y);
            return (x);
          default:
            return (0);
        }
      }

    case OP_MIN:
    case OP_XOR:
      return (eval_expression(t->lft,cnts) + eval_expression(t->rgt,cnts));

    case OP_CNT:
      { int x, *r;

        x = eval_expression(t->lft,cnts);
        r = (int *) (t->rgt);
        for (int i = 0; i < t->mode; i += 2)
          if (r[i] <= x && x <= r[i+1])
            return (x);
        return (0);
      }

    case OP_ARG:
      return (cnts[(int64) (t->lft)]);

    default:
      return (0);
  }
}


/****************************************************************************************
 *
 *  Streaming eval
 *
 *****************************************************************************************/

#define  COUNT_OF(p) (*((uint16 *) (p+kbyte)))

static inline int mycmp(uint8 *a, uint8 *b, int n)
{ while (n-- > 0)
    { if (*a++ != *b++)
        return (a[-1] < b[-1] ? -1 : 1);
    }
  return (0);
}

static void Merge(Kmer_Stream **T, int ntabs, Node *E)
{ int kbyte = T[0]->kbyte;

  uint8 *ptr[ntabs];
  int    itop, in[ntabs], cnt[ntabs], imin;
  int    c, v, x;
  int   *logic;

  logic = compile_expression(E,ntabs);

  for (c = 0; c < ntabs; c++)
    { ptr[c] = First_Kmer_Entry(T[c]);
      cnt[c] = 0;
    }
  while (1)
    { for (c = 0; c < ntabs; c++)
        if (ptr[c] != NULL)
          break;
      if (c >= ntabs)
        break;
      itop = 0;
      in[itop++] = imin = c;
      v = (1 << c);
      for (c++; c < ntabs; c++)
        { if (ptr[c] == NULL)
            continue;
          v = mycmp(ptr[c],ptr[imin],kbyte);
          if (v == 0)
            { in[itop++] = c;
              v |= (1 << c);
            }
          else if (v < 0)
            { itop = 0;
              in[itop++] = imin = c;
              v = (1 << c);
            }
        }

      if (logic[v])
        { for (c = 0; c < itop; c++)
            { x = in[c];
              cnt[x] = COUNT_OF(ptr[x]);
            }

          v = eval_expression(E,cnt);
          if (v > 0)
            { // write(f,ptr[x],kbyte);
              // write(f,&v,sizeof(short);
            }

          for (c = 0; c < itop; c++)
            { x = in[c];
              cnt[x] = 0;
              ptr[x] = Next_Kmer_Entry(T[x]);
            }
        }
      else
        for (c = 0; c < itop; c++)
          { x = in[c];
            ptr[x] = Next_Kmer_Entry(T[x]);
          }
    }
}


/****************************************************************************************
 *
 *  Main
 *
 *****************************************************************************************/

int main(int argc, char *argv[])
{ Kmer_Stream **T;
  int           ntabs;
  Node         *E;

  int DO_TABLE;
  int NTHREADS;
  int HIST_LOW, HIST_HGH;

  { int    i, j, k;
    int    flags[128];
    char  *eptr, *fptr;

    (void) flags;

    ARG_INIT("Logex");

    NTHREADS = 4;
    HIST_LOW = 0;
    DO_TABLE = 1;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
            break;
          case 'H':
            DO_TABLE = 0;
          case 'h':
            HIST_LOW = strtol(argv[i]+2,&eptr,10);
            if (eptr > argv[i]+2)
              { if (HIST_LOW < 1 || HIST_LOW > 0x7fff)
                  { fprintf(stderr,"\n%s: Histogram count %d is out of range\n",
                                   Prog_Name,HIST_LOW);
                    exit (1);
                  }
                if (*eptr == ':')
                  { HIST_HGH = strtol(eptr+1,&fptr,10);
                    if (fptr > eptr+1 && *fptr == '\0')
                      { if (HIST_LOW > HIST_HGH)
                          { fprintf(stderr,"\n%s: Histogram range is invalid\n",Prog_Name);
                            exit (1);
                          }
                        break;
                      }
                  }
                else if (*eptr == '\0')
                  { HIST_HGH = HIST_LOW;
                    HIST_LOW = 1;
                    break;
                  }
              }
            fprintf(stderr,"\n%s: Syntax of -h option invalid -h[<int(1)>:]<int>\n",Prog_Name);
            exit (1);
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
        }   
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc < 2)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -T: Use -T threads.\n");
        fprintf(stderr,"      -h: Generate histograms.\n");
        fprintf(stderr,"      -H: Generate histograms only, no tables.\n");
        exit (1);
      } 
  }   

  E = parse_expression(argv[1],&ntabs);
  if (E == NULL)
    exit (1);

  print_tree(E,0);

  { int *logic = compile_expression(E,ntabs);
    int  i;
    for (i = 0; i < (1<<ntabs); i++)
      printf(" %0*x: %d\n",(ntabs-1)/4+1,i,logic[i]);
  }

  if (ntabs != argc-2)
    { fprintf(stderr,"%s: # of arguments (%d) != # of supplied tables (%d)\n",
                     Prog_Name,ntabs,argc-2);
      exit (1);
    }

  { int c;

    T = Malloc(sizeof(Kmer_Stream *)*ntabs,"Allocating streams");
    if (T == NULL)
      exit (1);

    for (c = 0; c < ntabs; c++)
      T[c] = Open_Kmer_Stream(argv[c+2],1);

    Merge(T,ntabs,E);

    for (c = 0; c < ntabs; c++)
      Free_Kmer_Stream(T[c]);
    free(T);
  }

  free_tree(E);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
