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
#include <pthread.h>

#undef   DEBUG
#undef   DEBUG_THREADS

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
#define MOD_ONE 6

#ifdef DEBUG

static char *Operator[] =
  { "OR", "AND", "MIN", "XOR", "CNT", "NUM", "ARG" };

static char *Modulator[] =
  { "AVE", "SUM", "SUB", "MIN", "MAX", "LFT", "1" };

#endif

typedef struct _node
  { int16         op;
    int16         mode;
    struct _node *lft;
    struct _node *rgt;
  } Node;

static char *Scan;
static int   Error;
static int   VarV;
static int   hasFilter;
static char *hasNoMode;
static int   Narg;
static Node *or();

#define ERROR(msg)	\
{ Error = msg;		\
  return (NULL);	\
}

static char *Error_Messages[] =
  { "Out of memory",					// 0
    "Expecting an argument or (",			// 1
    "Expecting )",					// 2
    "Modeless operator iniside argument of []",		// 3
    "Premature end of expression",			// 4
    "Expecting a -",					// 5
    "Expecting a , or ]",				// 6
    "Do not recognize this operator",			// 7
    "Largest argument possible is H/h",			// 8
    "Modeless operator not in # argument",		// 9
    "Invalid modulator"	                		// 10
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
      else if (v->op == OP_NUM && v->mode)
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
    { int   hF;
      char *hM;

      Scan += 1;
      hF = hasFilter;
      hM = hasNoMode;
      hasFilter = 0;
      v = terminal();
      if (v == NULL)
        return (NULL);
      if (hasFilter)
        { hasFilter = 1;
          hasNoMode = hM;
          return (node(OP_NUM,0,v,NULL));
        }
      else
        { int *tab = Malloc(sizeof(int)*(1<<Narg),"Allocating logic table");
          if (tab == NULL)
            { free_tree(v);
              ERROR(0);
            }
          hasFilter = hF;
          hasNoMode = hM;
          return (node(OP_NUM,1,v,(Node *) tab));
        }
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
  char *hM;

  hM = hasNoMode;
  hasNoMode = NULL;
  v = num();
  if (v == NULL)
    return (NULL);

  while (isspace(*Scan))
    Scan += 1;
  while (*Scan == '[')
    { int   len, *f;
      char *beg;
      int   i, p;

      if (hasNoMode != NULL)
        { free_tree(v);
          Scan = hasNoMode;
          ERROR(3);
        }
      hasFilter = 1;

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

  if (v->op == OP_CNT || hM != NULL)
    hasNoMode = hM;
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
      if (*Scan == '(' || isalpha(*Scan) || *Scan == '#' || isspace(*Scan))
        return (MOD_ONE);
      else
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
      if (m == MOD_ONE)
        hasNoMode = Scan-1;
      else if (m >= 0)
        Scan += 1;
      else
        { free_tree(v);
          ERROR(10);
        }
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
      if (m == MOD_ONE)
        hasNoMode = Scan-1;
      else if (m >= 0)
        Scan += 1;
      else
        { free_tree(v);
          ERROR(10);
        }
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
      fflush(stdout);
      print_tree(v->lft,level+2);
    }
  else if (v->op == OP_ARG)
    { printf("%*s%s %lld\n",level,"",Operator[v->op],(int64) (v->lft));
      fflush(stdout);
    }
  else if (v->op == OP_NUM)
    { printf("%*s%s (%d)\n",level,"",Operator[v->op],v->mode);
      if (v->mode)
        for (int i = 0; i < (1<<Narg); i++)
          printf("%*s  %0*x: %d\n",level,"",(Narg-1)/4+1,i,((int *) (v->rgt))[i]);
      fflush(stdout);
      print_tree(v->lft,level+2);
    }
  else
    { printf("%*s%s",level,"",Operator[v->op]);
      if (v->op <= OP_AND)
        printf(" %s",Modulator[v->mode]);
      printf("\n");
      fflush(stdout);
      print_tree(v->lft,level+2);
      print_tree(v->rgt,level+2);
    }
}

#endif

static Node *parse_expression(char *expr, int *varg, int ntabs)
{ Node *v;

  hasNoMode = NULL;
  hasFilter = 0;
  Narg = ntabs;
  VarV = 0;
  Scan = expr;
  v    = or();
  if (v != NULL)
    { if (*Scan != '\0')
        { free_tree(v);
          Error = 7;
          v = NULL;
        }
      else if (hasNoMode)
        { free_tree(v);
          Scan = hasNoMode;
          Error = 9;
          v = NULL;
        }
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
{ int x, y;
  switch (t->op)
  { case OP_OR:
      x = eval_logic(t->lft,i);
      y = eval_logic(t->rgt,i);
      return (x || y);
    case OP_MIN:
      x = eval_logic(t->lft,i);
      y = eval_logic(t->rgt,i);
      return (x && (1-y));
    case OP_XOR:
      x = eval_logic(t->lft,i);
      y = eval_logic(t->rgt,i);
      return (x != y);
    case OP_AND:
      x = eval_logic(t->lft,i);
      y = eval_logic(t->rgt,i);
      return (x && y);
    case OP_CNT:
      return (eval_logic(t->lft,i));
    case OP_NUM:
      x = eval_logic(t->lft,i);
      if (t->mode)
        ((int *) (t->rgt))[i] = x;
      return (x);
    case OP_ARG:
      { x = (int64) (t->lft);
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

  if (t->op == OP_NUM)
    { for (i = (1 << ntabs)-1; i >= 0; i--)
        eval_logic(t,i);
      return ((int *) (t->rgt));
    }
  else
    { table = Malloc(sizeof(int)*(1<<ntabs),"Allocating logic table");
      for (i = (1 << ntabs)-1; i >= 0; i--)
        table[i] = eval_logic(t,i);
      return (table);
    }
}

static inline int modulate(int x, int y, int mode)
{ switch (mode)
  { case MOD_AVE:
      return ((x+y) >> 1);
    case MOD_SUM:
      return (x+y);
    case MOD_SUB:
      x -= y;
      if (x < 0)
        x = 0;
      return (x);
    case MOD_MIN:
      if (x < y)
        return (x);
      else
        return (y);
    case MOD_MAX:
      if (x > y)
        return (x);
      else
        return (y);
    case MOD_LFT:
      if (x == 0)
        return (y);
      return (x);
    case MOD_ONE:
    default:
      return (1);
  }
}

static int eval_expression(Node *t, int *cnts)
{ switch (t->op)
  { case OP_NUM:
      if (eval_expression(t->lft,cnts) > 0)
        return (1);
      else
        return (0);

    case OP_OR:
      { int x, y;

        x = eval_expression(t->lft,cnts);
        y = eval_expression(t->rgt,cnts);
        if (x == 0)
          return (y);
        if (y == 0)
          return (x);
        return (modulate(x,y,t->mode));
      }
 
    case OP_AND:
      { int x, y;

        x = eval_expression(t->lft,cnts);
        y = eval_expression(t->rgt,cnts);
        if (x == 0)
          return (0);
        if (y == 0)
          return (0);
        return (modulate(x,y,t->mode));
      }

    case OP_XOR:
      { int x, y;

        x = eval_expression(t->lft,cnts);
        y = eval_expression(t->rgt,cnts);
        if (x == 0 && y != 0)
          return (y);
        if (y == 0 && x != 0)
          return (x);
        return (0);
      }

    case OP_MIN:
      { int x, y;

        x = eval_expression(t->lft,cnts);
        y = eval_expression(t->rgt,cnts);
        if (x != 0 && y == 0)
          return (x);
        return (0);
      }

    case OP_CNT:
      { int x, *r;

        x = eval_expression(t->lft,cnts);
        r = (int *) (t->rgt);
        for (int i = 0; i < t->mode; i += 2)
          { if (x < r[i])
              return (0);
            else if (x <= r[i+1])
              return (x);
          }
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
    int  *filter;
    int   logical;
    int   ntabs;
  } Assignment;

static Assignment *parse_assignment(char *ass, int ntabs)
{ char       *expr;
  Assignment *A;

  A = Malloc(sizeof(Assignment),"Allocating assignment\n");
  if (A == NULL)
    exit (1);

  expr = index(ass,'=');
  while (isspace(expr[-1]))
    expr -= 1;
  *expr++ = '\0';

  A->ntabs   = ntabs;
  A->root    = Root(ass,".ktab");
  A->path    = PathTo(ass);
  A->expr    = parse_expression(expr+1,&(A->varg),ntabs);
  A->filter  = compile_expression(A->expr,ntabs);
  A->logical = (A->expr->op == OP_NUM);
  return (A);
}

#ifdef DEBUG

static void print_assignment(Assignment *A)
{ printf("'%s' '%s' %02x:\n",A->path,A->root,A->varg);
  Narg = A->ntabs;
  print_tree(A->expr,0);
  if ( ! A->logical)
    for (int i = 0; i < (1 << Narg); i++)
      printf(" %0*x: %d\n",(Narg-1)/4+1,i,A->filter[i]);
}

#endif

static void free_assignment(Assignment *A)
{ if ( ! A->logical)
    free(A->filter);
  free_tree(A->expr);
  free(A->path);
  free(A->root);
  free(A);
}

/****************************************************************************************
 *
 *  Streaming eval
 *
 *****************************************************************************************/

typedef struct
  { int           tid;
    char        **arg;
    Kmer_Stream **S;
    int           narg;
    Assignment  **A;
    int           nass;
    int64        *begs;
    int64        *ends;
    int64       **hist;
  } TP;

#define  COUNT_OF(p) (*((uint16 *) (p+kbyte)))

static inline int mycmp(uint8 *a, uint8 *b, int n)
{ while (n-- > 0)
    { if (*a++ != *b++)
        return (a[-1] < b[-1] ? -1 : 1);
    }
  return (0);
}

#ifdef DEBUG

static char dna[4] = { 'a', 'c', 'g', 't' };

static char *fmer[256], _fmer[1280];

static void setup_fmer_table()
{ char *t;
  int   i, l3, l2, l1, l0;

  i = 0;
  t = _fmer;
  for (l3 = 0; l3 < 4; l3++)
   for (l2 = 0; l2 < 4; l2++)
    for (l1 = 0; l1 < 4; l1++)
     for (l0 = 0; l0 < 4; l0++)
       { fmer[i] = t;
         *t++ = dna[l3];
         *t++ = dna[l2];
         *t++ = dna[l1];
         *t++ = dna[l0];
         *t++ = 0;
         i += 1;
       }
}

static void print_seq(FILE *out, uint8 *seq, int len)
{ int i, b, k;

  b = len >> 2;
  for (i = 0; i < b; i++)
    fprintf(out,"%s",fmer[seq[i]]);
  k = 6;
  for (i = b << 2; i < len; i++)
    { fprintf(out,"%c",dna[seq[b] >> k]);
      k -= 2;
    }
}

#endif

static void *merge_thread(void *args)
{ TP *parm = (TP *) args;
  int           tid   = parm->tid;
  char        **arg   = parm->arg;
  Assignment  **A     = parm->A;
  int           ntabs = parm->narg;
  int           nass  = parm->nass;
  Kmer_Stream **T     = parm->S;
  int64        *begs  = parm->begs;
  int64        *ends  = parm->ends;

  int one   = 1;
  int kbyte = T[0]->kbyte;
  int kmer  = T[0]->kmer;
  int hgram = (HIST_LOW > 0);

  uint8 **ptr, *bp;
  FILE  **out;
  int64 **hist, *nels;
  int    *filter, need_counts;
  int    itop, *in, *cnt;
  int    c, v, x, i;

#ifdef DEBUG
  setup_fmer_table();
#endif

  ptr    = Malloc(sizeof(uint8 *)*ntabs,"Allocating thread working memory");
  in     = Malloc(sizeof(int)*ntabs,"Allocating thread working memory");
  cnt    = Malloc(sizeof(int)*ntabs,"Allocating thread working memory");
  nels   = Malloc(sizeof(int64)*nass,"Allocating thread working memory");
  filter = Malloc(sizeof(int)*(1<<ntabs),"Allocating thread working memory");

#ifdef DEBUG_THREADS
  printf("Doing %d:",tid);
  for (c = 0; c < ntabs; c++)
    printf(" [%lld-%lld]",begs[c],ends[c]);
  printf("\n");
#endif

  if (tid != 0)
    { T = Malloc(sizeof(Kmer_Stream *)*ntabs,"Allocating thread working memory");
      for (c = 1; c <= ntabs; c++)
        T[c-1] = Open_Kmer_Stream(arg[c]);
    }
  else if (DO_TABLE)
    { for (i = 0; i < nass; i++)
        { FILE *f = fopen(Catenate(A[i]->path,"/",A[i]->root,".ktab"),"w");
          fwrite(&kmer,sizeof(int),1,f);
          fwrite(&NTHREADS,sizeof(int),1,f);
          fwrite(&one,sizeof(int),1,f);
          fclose(f);
        }
    }

  need_counts = 0;
  for (v = 0; v < (1 << ntabs); v++)
    filter[v] = 0;
  for (i = 0; i < nass; i++)
    { int *table = A[i]->filter;
      for (v = 0; v < (1 << ntabs); v++)
        filter[v] |= table[v];
      if ( ! A[i]->logical)
        need_counts = 1;
    }

  if (DO_TABLE)
    { out = Malloc(sizeof(FILE *)*nass,"Allocating thread working memory");
      for (i = 0; i < nass; i++)
        { out[i] = fopen(Catenate(A[i]->path,"/.",A[i]->root,
                                        Numbered_Suffix(".ktab.",tid+1,"")),"w");
          nels[i] = 0;
          fwrite(&kmer,sizeof(int),1,out[i]);
          fwrite(nels+i,sizeof(int64),1,out[i]);
        }
    }

  if (hgram)
    { hist = Malloc(sizeof(int64 *)*nass,"Allocating thread working memory");
      for (i = 0; i < nass; i++)
        { hist[i] = Malloc(sizeof(int64)*((HIST_HGH-HIST_LOW)+1),"Allocating histogram");
          hist[i] -= HIST_LOW;
        }
    }

  for (c = 0; c < ntabs; c++)
    { ptr[c] = GoTo_Kmer_Index(T[c],begs[c]);
      cnt[c] = 0;
    }

  while (1)
    { for (c = 0; c < ntabs; c++)
        if (T[c]->cidx < ends[c])
          break;
      if (c >= ntabs)
        break;
      itop  = 1;
      in[0] = c;
      bp    = ptr[c];
      v = (1 << c);
      for (c++; c < ntabs; c++)
        { if (T[c]->cidx >= ends[c])
            continue;
          x = mycmp(ptr[c],bp,kbyte);
          if (x == 0)
            { in[itop++] = c;
              v |= (1 << c);
            }
          else if (x < 0)
            { itop  = 1;
              in[0] = c;
              bp    = ptr[c];
              v = (1 << c);
            }
        }

      if (filter[v])
        { if (need_counts)
            { for (c = 0; c < itop; c++)
                { x = in[c];
                  cnt[x] = COUNT_OF(ptr[x]);
                }
            }

          for (i = 0; i < nass; i++)
            if (A[i]->filter[v])
              { if (A[i]->logical)
                  { if (DO_TABLE)
                      { fwrite(bp,kbyte,1,out[i]);
                        fwrite(&one,sizeof(short),1,out[i]);
                        nels[i] += 1;
                      }
                    if (hgram)
                      hist[i][HIST_LOW] += 1;
                  }
                else
                  { v = eval_expression(A[i]->expr,cnt);
                    if (v > 0)
                      { if (DO_TABLE)
                          { fwrite(bp,kbyte,1,out[i]);
                            fwrite(&v,sizeof(short),1,out[i]);
// print_seq(stdout,bp,kmer);
// printf(" %d %d %d\n",v,i,kbyte);
                            nels[i] += 1;
                          }
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
              }

          if (need_counts)
            for (c = 0; c < itop; c++)
              { x = in[c];
                cnt[x] = 0;
                ptr[x] = Next_Kmer_Entry(T[x]);
              }
          else
            for (c = 0; c < itop; c++)
              { x = in[c];
                ptr[x] = Next_Kmer_Entry(T[x]);
              }
        }
      else
        for (c = 0; c < itop; c++)
          { x = in[c];
            ptr[x] = Next_Kmer_Entry(T[x]);
          }
    }

  if (DO_TABLE)
    { for (i = 0; i < nass; i++)
        { rewind(out[i]);
          fwrite(&kmer,sizeof(int),1,out[i]);
          fwrite(nels+i,sizeof(int64),1,out[i]);
          fclose(out[i]);
        }
      free(out);
    }

  if (tid != 0)
    { for (c = 0; c < ntabs; c++) 
        Free_Kmer_Stream(T[c]);
      free(T);
    }

  free(filter);
  free(nels);
  free(cnt);
  free(in);
  free(ptr);

  parm->hist = hist;
  return (NULL);
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
      { Assignment *a = parse_assignment(argv[c],narg);
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
  }

  { int64     range[NTHREADS+1][narg];
    pthread_t threads[NTHREADS];
    TP        parm[NTHREADS];
    int       t, a, i;
    int64     p;

    for (a = 0; a < narg; a++)
      { range[0][a] = 0;
        range[NTHREADS][a] = S[a]->nels;
      }
    for (t = 1; t < NTHREADS; t++)
      { p = (S[0]->nels*t)/NTHREADS; 
        GoTo_Kmer_Index(S[0],p);
        printf(" %lld: %s\n",p,Current_Kmer(S[0]));
        range[t][0] = p;
        for (a = 1; a < narg; a++)
          { GoTo_Kmer_String(S[a],S[0]->celm);
            printf(" %lld: %s\n",S[a]->cidx,Current_Kmer(S[a]));
            range[t][a] = S[a]->cidx;
          }
      }

    for (t = 0; t < NTHREADS; t++)
      { parm[t].tid  = t;
        parm[t].arg  = argv+nass;
        parm[t].S    = S;
        parm[t].narg = narg;
        parm[t].A    = A;
        parm[t].nass = nass;
        parm[t].begs = range[t];
        parm[t].ends = range[t+1];
      }

#ifdef DEBUG_THREADS
    for (t = 0; t < NTHREADS; t++)
      merge_thread(parm+t);
#else
    for (t = 1; t < NTHREADS; t++)
      pthread_create(threads+t,NULL,merge_thread,parm+t);
    merge_thread(parm);
    for (t = 1; t < NTHREADS; t++)
      pthread_join(threads[t],NULL);
#endif

    if (HIST_LOW > 0)
      { int64 *histi, *histt;
        FILE  *f;

        for (i = 0; i < nass; i++)
          { histi = parm[0].hist[i];
            for (t = 1; t < NTHREADS; t++)
              { histt = parm[t].hist[i];
                for (a = HIST_LOW; a <= HIST_HGH; a++)
                  histi[a] += histt[a];
                free(histt);
              }

            f = fopen(Catenate(A[i]->path,"/",A[i]->root,".hist"),"w");
            fwrite(&(S[0]->kmer),sizeof(int),1,f);
            fwrite(&HIST_LOW,sizeof(int),1,f);
            fwrite(&HIST_HGH,sizeof(int),1,f);
            fwrite(histi,sizeof(int64),(HIST_HGH-HIST_LOW)+1,f);
            fclose(f);

            free(histi);
          }

        for (t = 0; t < NTHREADS; t++)
          free(parm[t].hist);
      }
  }

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
