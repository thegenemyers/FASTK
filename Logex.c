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

#define DEBUG

#include "libfastk.h"

static char *Usage[] = { " [-T<int(4)>] [-[hH]{<int(1)>:]<int>]",
                         "   <output:name=expr> ... <source_root>[.ktab] ..." };

#define MAX_TABS 8

static int DO_TABLE;
static int NTHREADS;
static int HIST_LOW, HIST_HGH;

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
#define OP_NUM 5
#define OP_ARG 6

#define MOD_AVE 0
#define MOD_SUM 1
#define MOD_SUB 2
#define MOD_MIN 3
#define MOD_MAX 4
#define MOD_LFT 5

#ifdef DEBUG

static char *Operator[] =
  { "OR", "AND", "MIN", "XOR", "CNT", "NUM", "ARG" };

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
      if (x > MAX_TABS)
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

static Node *num()
{ Node *v;

  while (isspace(*Scan))
    Scan += 1;
  if (*Scan == '#')
    { Scan += 1;
      v = terminal();
      if (v == NULL)
        return (NULL);
      else
        return (node(OP_NUM,0,v,NULL));
    }
  else
    return (terminal());
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

static int RSORT(const void *l, const void *r)
{ int *x = (int *) l;
  int *y = (int *) r;
  return (*x-*y);
}

static Node *filter()
{ Node *v;

  v = num();
  if (v == NULL)
    return (NULL);

  while (isspace(*Scan))
    Scan += 1;
  while (*Scan == '[')
    { int   len, *f;
      char *beg;
      int   i, p;

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

      qsort(f,len/2,2*sizeof(int),RSORT);

      p = 0;
      for (i = 2; i < len; i += 2)
        if (f[i] <= f[p+1])
          { if (f[i+1] > f[p+1])
              f[p+1] = f[i+1];
          }  
        else
          { p += 2;
            f[p] = f[i];
            f[p+1] = f[i+1];
          }
      len = p+2;

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
      if (m >= 0)
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
      if (m >= 0)
        Scan += 1;
      w = minus();
      if (w == NULL)
        { free_tree(v);
          return (NULL);
        }
      v = node(OP_OR,m,v,w);
    }
}

#ifdef DEBUG

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
  else if (v->op == OP_NUM)
    { printf("%*s%s\n",level,"",Operator[v->op]);
      print_tree(v->lft,level+2);
    }
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

static Node *parse_expression(char *expr, int *varg)
{ Node *v;

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
      exit (1);
    }

  *varg = VarV;
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
    case OP_NUM:
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
  { case OP_NUM:
      if (eval_expression(t->lft,cnts) > 0)
        return (1);
      else
        return (0);

    case OP_OR:
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
 *  Output Assignment
 *
 *****************************************************************************************/

typedef struct
  { Node *expr;
    int   varg;
    char *root;
    char *path;
    int  *logic;
  } Assignment;


static Assignment *parse_assignment(char *ass)
{ char       *expr;
  Assignment *A;

  A = Malloc(sizeof(Assignment),"Allocating assignment\n");
  if (A == NULL)
    exit (1);

  expr = index(ass,'=');
  while (isspace(expr[-1]))
    expr -= 1;
  *expr++ = '\0';

  A->root  = Root(ass,".ktab");
  A->path  = PathTo(ass);
  A->expr  = parse_expression(expr+1,&(A->varg));
  A->logic = NULL;

  return (A);
}

static void print_assignment(Assignment *A)
{ printf("'%s' '%s' %02x:\n",A->path,A->root,A->varg);
  print_tree(A->expr,0);
}

static void free_assignment(Assignment *A)
{ free(A->logic);
  free(A->expr);
  free(A->path);
  free(A->root);
  free(A);
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

static void Merge(Kmer_Stream **T, int ntabs, Assignment **A, int nass)
{ static int one = 1;

  int kbyte = T[0]->kbyte;
  int kmer  = T[0]->kmer;
  int hgram = (HIST_LOW > 0);

  uint8 *ptr[ntabs];
  FILE  *out[nass];
  int64 *hist[nass];
  int    any_logic[256], any_count[256];
  int    itop, in[ntabs], cnt[ntabs], imin;
  int    c, v, x, i;

  for (c = 0; c < ntabs; c++)
    { ptr[c] = First_Kmer_Entry(T[c]);
      cnt[c] = 0;
    }

  for (v = 0; v < (1 << ntabs); v++)
    any_logic[v] = any_count[v] = 0;
  for (i = 0; i < nass; i++)
    { int *logic = A[i]->logic;
      for (v = 0; v < (1 << ntabs); v++)
        any_logic[v] |= logic[v];
      if (A[i]->expr->op != OP_NUM)
        for (v = 0; v < (1 << ntabs); v++)
          any_count[v] |= logic[v];
    }

  for (i = 0; i < nass; i++)
    { if (DO_TABLE)
        { out[i] = fopen(Catenate(A[i]->path,"/.",A[i]->root,".ktab"),"w");
          fwrite(&kmer,sizeof(int),1,out[i]);
          fwrite(&one,sizeof(int),1,out[i]);
          fwrite(&one,sizeof(int),1,out[i]);
        }
      if (hgram)
        { hist[i] = Malloc(sizeof(int64)*((HIST_HGH-HIST_LOW)+1),"Allocating histogram");
          hist[i] -= HIST_LOW;
        }
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

      if (any_logic[v])
        { if (any_count[v])
            { for (c = 0; c < itop; c++)
                { x = in[c];
                  cnt[x] = COUNT_OF(ptr[x]);
                }

              for (i = 0; i < nass; i++)
                if (A[i]->logic[v])
                  { v = eval_expression(A[0]->expr,cnt);
                    if (v > 0)
                      { fwrite(ptr[x],kbyte,1,out[i]);
                        fwrite(&v,sizeof(short),1,out[i]);
                        if (hgram)
                          { if (v < HIST_LOW)
                              hist[i][HIST_LOW] += v;
                            else if (v > HIST_HGH)
                              hist[i][HIST_HGH] += v;
                            else
                              hist[i][v] += v;
                          }
                      }
                  }

              for (c = 0; c < itop; c++)
                { x = in[c];
                  cnt[x] = 0;
                  ptr[x] = Next_Kmer_Entry(T[x]);
                }
            }
          else
            { for (i = 0; i < nass; i++)
                if (A[i]->logic[v])
                  { fwrite(ptr[x],kbyte,1,out[i]);
                    fwrite(&v,sizeof(short),1,out[i]);
                    if (hgram)
                      hist[i][HIST_LOW] += 1;
                  }
            }
        }
      else
        for (c = 0; c < itop; c++)
          { x = in[c];
            ptr[x] = Next_Kmer_Entry(T[x]);
          }
    }

  if (DO_TABLE)
    for (i = 0; i < nass; i++)
      fclose(out[i]);
}


/****************************************************************************************
 *
 *  Main
 *
 *****************************************************************************************/

int main(int argc, char *argv[])
{ int           nass, narg;
  Assignment  **A;
  Kmer_Stream **S;

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
  
  { int c, kmer;

    for (c = 1; c < argc; c++)
      if (index(argv[c],'=') == NULL)
        break;
    nass = c-1;
    narg = argc-c;
    if (nass == 0)
      { fprintf(stderr,"%s: There must be at least one assignment argument\n",Prog_Name);
        exit (1);
      }
    if (narg == 0)
      { fprintf(stderr,"%s: There must be at least one table argument\n",Prog_Name);
        exit (1);
      }
    if (narg > MAX_TABS)
      { fprintf(stderr,"%s: So sorry, but at most %d tables are possible.\n",Prog_Name,MAX_TABS);
        exit (1);
      }

    A = Malloc(sizeof(Kmer_Stream *)*nass,"Allocating assignment pointers");
    S = Malloc(sizeof(Kmer_Stream *)*narg,"Allocating table pointers");
    if (A == NULL || S == NULL)
      exit (1);

    for (c = 1; c <= nass; c++)
      { Assignment *a = parse_assignment(argv[c]);
        if (a == NULL)
          exit (1);
        if (DO_TABLE)
          { FILE *f;
            int   x, yes;

            f = fopen(Catenate(a->path,"/",a->root,".ktab"),"r");
            if (f != NULL)
              { fclose(f);
                printf("Output table %s already exists, continue? ",
                       Catenate(a->path,"/",a->root,".ktab"));
                fflush(stdout);
                yes = 0;
                while ((x = getc(stdin)) != '\n')
                  if (x == 'y' || x == 'Y')
                    yes = 1;
                  else if (x == 'n' || x == 'N')
                    yes = 0;
                if (!yes)
                  exit (1);
              }
          }
        A[c-1] = a;
#ifdef DEBUG
        print_assignment(a);
#endif
      }

    kmer = 0;
    for (c = 1; c <= narg; c++)
      { Kmer_Stream *s = Open_Kmer_Stream(argv[c+nass]);
        if (s == NULL)
          { fprintf(stderr,"%s: Cannot open table %s\n",Prog_Name,argv[c+nass]);
            exit (1);
          }
        if (c == 1)
          kmer = s->kmer;
        else
          { if (s->kmer != kmer)
              { fprintf(stderr,"%s: K-mer tables do not involve the same K\n",Prog_Name);
                exit (1);
              }
          }
        S[c-1] = s;
      }
  }

  { int c, varg;

    varg = 0;
    for (c = 0; c < nass; c++) 
      varg |= A[c]->varg;

    if (varg >= (1 << narg))
      { if (nass == 1)
          fprintf(stderr,"%s: Expression refers to tables not given\n",Prog_Name);
        else
          fprintf(stderr,"%s: Expressions refer to tables not given\n",Prog_Name);
        exit (1);
      }
    if (varg < (1 << (narg-1)))
      { fprintf(stderr,"%s: There are tables not referred to by an expression\n",Prog_Name);
        exit (1);
      }
    if (varg != (1 << narg)-1)
      { if (nass == 1)
          fprintf(stderr,"%s: Expression does not refer ta all the tables\n",Prog_Name);
        else
          fprintf(stderr,"%s: Expressions do not refer ta all the tables\n",Prog_Name);
        exit (1);
      }

    for (c = 0; c < nass; c++)
      { int *logic = compile_expression(A[c]->expr,narg);
        A[c]->logic = logic;
#ifdef DEBUG
        for (int i = 0; i < (1<<narg); i++)
          printf(" %0*x: %d\n",(narg-1)/4+1,i,logic[i]);
#endif
       }
  }

  { int   i;
    int64 p;

    for (i = 1; i < NTHREADS; i++)
      { p = (S[0]->nels*i)/NTHREADS; 
        GoTo_Kmer_Index(S[0],p);
        printf(" %lld: %s\n",p,Current_Kmer(S[0]));
      }
    First_Kmer_Entry(S[0]);
    printf(" %dll: %s\n",0,Current_Kmer(S[0]));

    if (narg > 1)
      { for (i = 1; i < NTHREADS; i++)
          { p = (S[0]->nels*i)/NTHREADS; 
            GoTo_Kmer_Index(S[0],p);
            GoTo_Kmer_String(S[1],S[0]->celm);
            printf(" %lld: %s\n",S[1]->cidx,Current_Kmer(S[1]));
          }
      }
  }

  Merge(S,narg,A,nass);

  { int c;

    for (c = 0; c < narg; c++)
      Free_Kmer_Stream(S[c]);
    free(S);

    for (c = 0; c < nass; c++)
      free_assignment(A[c]);
    free(A);
  }

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
