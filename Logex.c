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

#undef  DEBUG
#undef  DEBUG_THREADS
#undef  DEBUG_TRACE

#include "libfastk.h"

static char *Usage[] = { " [-T<int(4)>] [-[hH][<int(1)>:]<int>]",
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
#define OP_GC  5
#define OP_NUM 6
#define OP_ARG 7

#define MOD_AVE 0
#define MOD_SUM 1
#define MOD_SUB 2
#define MOD_MIN 3
#define MOD_MAX 4
#define MOD_LFT 5
#define MOD_ONE 6

#ifdef DEBUG

static char *Operator[] =
  { "OR", "AND", "MIN", "XOR", "CNT", "GC", "NUM", "ARG" };

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
static int   hasFilter;   //  There is a filter operator in the expression
static char *hasNoMode;   //  Ptr @ modeless op in current #-rooted subtree if one, NULL otherwise
static int   needGC;      //  There is a {} filter operator in the expression
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
    "Modeless operator inside argument of [] or {}",	// 3
    "Premature end of expression",			// 4
    "Expecting a -",					// 5
    "Expecting a , or ]",				// 6
    "Do not recognize this operator",			// 7
    "Largest argument possible is H/h",			// 8
    "Modeless operator not in # argument",		// 9
    "Invalid modulator"	                		// 10
    "Expecting a , or }",				// 11
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
      if (v->op == OP_CNT || v->op == OP_GC)
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
  char *hM, Close;

  hM = hasNoMode;
  hasNoMode = NULL;
  v = num();
  if (v == NULL)
    return (NULL);

  while (isspace(*Scan))
    Scan += 1;
  while (*Scan == '[' || *Scan == '{')
    { int   len, *f;
      char *beg;
      int   i, p;

      if (*Scan == '[')
        Close = ']';
      else
        { Close = '}';
          needGC = 1;
        }

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
        { if (*Scan != '-')
            { if (isdigit(*Scan))
                get_number();
              if (*Scan == Close)
                { len += 2;
                  break;
                }
              if (*Scan == ',')
                { len += 2;
                  Scan += 1;
                  while (isspace(*Scan))
                    Scan += 1;
                  continue;
                }
              if (*Scan != '-')
                { free_tree(v);
                  ERROR(5)
                }
            }
          Scan += 1;
          while (isspace(*Scan))
            Scan += 1;
          if (isdigit(*Scan))
            get_number();
          len += 2;
          if (*Scan == Close)
            break;
          if (*Scan != ',')
            { free_tree(v);
              if (Close == ']')
                ERROR(6)
              else
                ERROR(11)
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
        { if (*Scan != '-')
            { if (isdigit(*Scan))
                { f[len] = f[len+1] = get_number();
                  len += 2;
                }
              if (*Scan == Close)
                break;
              if (*Scan == ',')
                { Scan += 1;
                  while (isspace(*Scan))
                    Scan += 1;
                  continue;
                }
            }
          else
            { f[len] = f[len+1] = 1;
              len += 2;
            }
          Scan += 1;
          while (isspace(*Scan))
            Scan += 1;
          if (isdigit(*Scan))
            f[len-1] = get_number();
          else
            f[len-1] = 0x7fff;
          if (*Scan == Close)
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

      if (Close == '}')
        v = node(OP_GC,len,v,(Node *) f);
      else
        v = node(OP_CNT,len,v,(Node *) f);
    } 

  if (v->op == OP_CNT || v->op == OP_GC || hM != NULL)
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
{ if (v->op == OP_CNT || v->op == OP_GC)
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
    { int i;

      printf("%*s%s (%d)\n",level,"",Operator[v->op],v->mode);
      if (v->mode)
        for (i = 0; i < (1<<Narg); i++)
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
    case OP_GC:
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
      if (t->rgt == NULL)
        { if (eval_expression(t->lft,cnts) > 0)
            return (1);
          else
            return (0);
        }
      else
        return (((int *) (t->rgt))[cnts[-1]]);

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
      { int i, x, *r;

        x = eval_expression(t->lft,cnts);
        r = (int *) (t->rgt);
        for (i = 0; i < t->mode; i += 2)
          { if (x < r[i])
              return (0);
            else if (x <= r[i+1])
              return (x);
          }
        return (0);
      }

    case OP_GC:
      { int i, x, *r;

        x = cnts[-2];
        r = (int *) (t->rgt);
        for (i = 0; i < t->mode; i += 2)
          { if (x < r[i])
              return (0);
            else if (x <= r[i+1])
              return (eval_expression(t->lft,cnts));
          }
        return (0);
      }
     
    case OP_ARG:
      return (cnts[(int64) (t->lft)]);

    default:
      return (0);
  }
}

static inline int modulate_mins(int x, int y, int mode)
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
      return (x);
    case MOD_ONE:
    default:
      return (1);
  }
}

static int eval_minimums(Node *t, int *mins)
{ switch (t->op)
  { case OP_NUM:
      return (1);

    case OP_OR:
      { int x, y;

        x = eval_minimums(t->lft,mins);
        y = eval_minimums(t->rgt,mins);
        if (y < x)
          x = y;
        y = modulate(x,y,t->mode);
        if (y < x)
          x = y;
        return (x);
      }
 
    case OP_AND:
      { int x, y;

        x = eval_minimums(t->lft,mins);
        y = eval_minimums(t->rgt,mins);
        return (modulate_mins(x,y,t->mode));
      }

    case OP_XOR:
      { int x, y;

        x = eval_minimums(t->lft,mins);
        y = eval_minimums(t->rgt,mins);
        if (x < y)
          return (x);
        else
          return (y);
      }

    case OP_MIN:
      return (eval_minimums(t->lft,mins));

    case OP_CNT:
      { int x, *r;

        x = eval_minimums(t->lft,mins);
        r = (int *) (t->rgt);
        if (x < r[0])
          return (r[0]);
        else
          return (x);
      }

    case OP_GC:
      return (eval_minimums(t->lft,mins));
     
    case OP_ARG:
      return (mins[(int64) (t->lft)]);

    default:
      return (1);
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
    int   needGC;
    int   ntabs;
  } Assignment;

static Assignment *parse_assignment(char *ass, int ntabs)
{ char       *expr, *s;
  Assignment *A;

  A = Malloc(sizeof(Assignment),"Allocating assignment\n");
  if (A == NULL)
    exit (1);

  expr = s = index(ass,'=');
  while (isspace(s[-1]))
    s -= 1;
  *s++ = '\0';

  A->ntabs   = ntabs;
  A->root    = Root(ass,".ktab");
  A->path    = PathTo(ass);

  needGC = 0;
  A->expr    = parse_expression(expr+1,&(A->varg),ntabs);
  A->needGC  = needGC;
  A->filter  = compile_expression(A->expr,ntabs);
  A->logical = (A->expr->op == OP_NUM && A->expr->rgt != NULL);
  return (A);
}

#ifdef DEBUG

static void print_assignment(Assignment *A)
{ int i;

  printf("'%s' '%s' %02x:\n",A->path,A->root,A->varg);
  Narg = A->ntabs;
  print_tree(A->expr,0);
  if ( ! A->logical)
    for (i = 0; i < (1 << Narg); i++)
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
    Kmer_Stream **S;
    int           narg;
    Assignment  **A;
    int           nass;
    FILE        **out;
    int64       **prefx;
    int64        *begs;
    int64        *ends;
    int64       **hist;
  } TP;

static int    GC[256];
static int    GCR[256];

#define IB_OUT  3

static void gc_setup(int kmer)
{ static int isgc[4] = { 0, 100, 100, 0 };
  uint32 x;

  for (x = 0; x < 256; x++)
    { GC[x] = isgc[(x>>6)&0x3] + isgc[(x>>4)&0x3] + isgc[(x>>2)&0x3] + isgc[x&0x3];
      switch( kmer % 4)
      { case 0:
          GCR[x] = GC[x];
          break;
        case 1:
          GCR[x] = isgc[(x>>6)&0x3];
          break;
        case 2:
          GCR[x] = isgc[(x>>6)&0x3] + isgc[(x>>4)&0x3];
          break;
        case 3:
          GCR[x] = isgc[(x>>6)&0x3] + isgc[(x>>4)&0x3] + isgc[(x>>2)&0x3];
          break;
      }
    }
}

static inline int gcontent(uint8 *a, int kbyte)
{ int i, cnt;

  cnt = 0;
  for (i = 1; i < kbyte; i++)
    cnt += GC[*a++];
  return (cnt+GCR[*a]);
}

static inline int mycmp(uint8 *a, uint8 *b, int n)
{ while (n-- > 0)
    { if (*a++ != *b++)
        return (a[-1] < b[-1] ? -1 : 1);
    }
  return (0);
}

static void *merge_thread(void *args)
{ TP *parm = (TP *) args;
  int           tid   = parm->tid;
  Assignment  **A     = parm->A;
  int           ntabs = parm->narg;
  int           nass  = parm->nass;
  Kmer_Stream **S     = parm->S;
  Kmer_Stream **T;
  int64        *begs  = parm->begs;
  int64        *ends  = parm->ends;
  int64       **prefx = parm->prefx;
  FILE        **out   = parm->out;

  int one   = 1;
  int hbyte = S[0]->kbyte-IB_OUT;
  int kbyte = S[0]->kbyte;
  int kmer  = S[0]->kmer;
  int hgram = (HIST_LOW > 0);

  int64 **hist = NULL;
  int64  *nels;
  uint8 **ent, *bst;
  uint16  sho;
  int    *filter, need_counts, need_GC;
  int    itop, *in, *cnt;
  int    c, v, x, i;

#ifdef DEBUG_TRACE
  char *buffer;
#endif

  if (hgram)
    { hist = Malloc(sizeof(int64 *)*nass,"Allocating thread working memory");
      hist[0] = Malloc(sizeof(int64)*((HIST_HGH-HIST_LOW)+3)*nass,"Allocating histogram");
      bzero(hist[0],sizeof(int64)*((HIST_HGH-HIST_LOW)+3)*nass);
      hist[0] -= HIST_LOW;
      for (i = 1; i < nass; i++)
        hist[i] = hist[i-1] + ((HIST_HGH-HIST_LOW)+3);
    }

  in     = Malloc(sizeof(int)*ntabs,"Allocating thread working memory");
  cnt    = Malloc(sizeof(int)*(ntabs+2),"Allocating thread working memory");
  nels   = Malloc(sizeof(int64)*nass,"Allocating thread working memory");
  filter = Malloc(sizeof(int)*(1<<ntabs),"Allocating thread working memory");
  ent    = Malloc(sizeof(uint8 *)*ntabs,"Allocating thread working memory");

  cnt += 2;  // cnt[-1] = v cnt[-2] = gc percent (if needed)

#ifdef DEBUG_THREADS
  printf("Doing %d:",tid);
  for (c = 0; c < ntabs; c++)
    printf(" [%lld-%lld]",begs[c],ends[c]);
  printf("\n");
#endif

  if (tid != 0)
    { T = Malloc(sizeof(Kmer_Stream *)*ntabs,"Allocating thread working memory");
      for (c = 0; c < ntabs; c++)
        T[c] = Clone_Kmer_Stream(S[c]);
    }
  else
    T = S;

#ifdef DEBUG_TRACE
  buffer = Current_Kmer(T[0],NULL);
#endif

  need_counts = 0;
  need_GC     = 0;
  for (v = 0; v < (1 << ntabs); v++)
    filter[v] = 0;
  for (i = 0; i < nass; i++)
    { int *table = A[i]->filter;
      for (v = 0; v < (1 << ntabs); v++)
        filter[v] |= table[v];
      if ( ! A[i]->logical)
        { need_counts = 1;
          if (A[i]->needGC)
            need_GC = 1;
        }
    }

  if (DO_TABLE)
    { for (i = 0; i < nass; i++)
        { nels[i] = 0;
          fwrite(&kmer,sizeof(int),1,out[i]);
          fwrite(nels+i,sizeof(int64),1,out[i]);
        }
    }

  for (c = 0; c < ntabs; c++)
    { GoTo_Kmer_Index(T[c],begs[c]);
      cnt[c] = 0;
    }

  for (c = 0; c < ntabs; c++)
    ent[c] = Current_Entry(T[c],NULL);

  while (1)
    { for (c = 0; c < ntabs; c++)
        if (T[c]->cidx < ends[c])
          break;
      if (c >= ntabs)
        break;
      itop  = 1;
      in[0] = c;
      bst = ent[c];
      v = (1 << c);
      for (c++; c < ntabs; c++)
        { if (T[c]->cidx >= ends[c])
            continue;
          x = mycmp(ent[c],bst,kbyte);
          if (x == 0)
            { in[itop++] = c;
              v |= (1 << c);
            }
          else if (x < 0)
            { itop  = 1;
              in[0] = c;
              bst = ent[c];
              v = (1 << c);
            }
        }

#ifdef DEBUG_TRACE
      for (c = 0; c < itop; c++)
        { x = in[c];
          printf(" %d: %s %5d",x,Current_Kmer(T[x],buffer),Current_Count(T[x]));
        }
      printf(" %x %d\n",v,filter[v]);
#endif

      if (filter[v])
        { if (need_counts)
            { for (c = 0; c < itop; c++)
                { x = in[c];
                  cnt[x] = Current_Count(T[x]);
                }
              cnt[-1] = v;
              if (need_GC)
                cnt[-2] = gcontent(bst,kbyte)/kmer;
            }

          for (i = 0; i < nass; i++)
            if (A[i]->filter[v])
              { if (A[i]->logical)
                  { if (DO_TABLE)
                      { fwrite(bst+IB_OUT,hbyte,1,out[i]);
                        fwrite(&one,sizeof(short),1,out[i]);
                        x = (bst[0] << 16) | (bst[1] << 8) | bst[2];
                        prefx[i][x] += 1;
                        nels[i] += 1;
                      }
                    if (hgram)
                      { hist[i][HIST_LOW] += 1;
                        hist[i][HIST_HGH+1] += 1;
                      }
#ifdef DEBUG_TRACE
                    printf("   %d: *\n",i);
#endif
                  }
                else
                  { c = eval_expression(A[i]->expr,cnt);
                    if (c > 0)
                      { if (DO_TABLE)
                          { fwrite(bst+IB_OUT,hbyte,1,out[i]);
                            if (c > 32767)
                              sho = 32767;
                            else
                              sho = c;
                            fwrite(&sho,sizeof(short),1,out[i]);
                            x = (bst[0] << 16) | (bst[1] << 8) | bst[2];
                            prefx[i][x] += 1;
                            nels[i] += 1;
                          }
                        if (hgram)
                          { if (c >= HIST_HGH)
                              { hist[i][HIST_HGH] += 1;
                                hist[i][HIST_HGH+2] += c;
                              }
                            else if (c <= HIST_LOW)
                              { hist[i][HIST_LOW] += 1;
                                hist[i][HIST_HGH+1] += c;
                              }
                            else
                              hist[i][c] += 1;
                          }
#ifdef DEBUG_TRACE
                        printf("   %d: %5d\n",i,c);
#endif
                      }
                  }
              }
        }

      for (c = 0; c < itop; c++)
        { Kmer_Stream *t;

          x = in[c];
          cnt[x] = 0;
          t = T[x];
          Next_Kmer_Entry(t);
          if (t->csuf != NULL)
            Current_Entry(t,ent[x]);
        }
    }

  if (DO_TABLE)
    { for (i = 0; i < nass; i++)
        { rewind(out[i]);
          fwrite(&kmer,sizeof(int),1,out[i]);
          fwrite(nels+i,sizeof(int64),1,out[i]);
          fclose(out[i]);
        }
    }

  for (c = 0; c < ntabs; c++)
    free(ent[c]);

  if (tid != 0)
    { for (c = 0; c < ntabs; c++) 
        Free_Kmer_Stream(T[c]);
      free(T);
    }

  free(ent);
  free(filter);
  free(nels);
  free(cnt-2);
  free(in);

  parm->hist = hist;
  return (NULL);
}


/****************************************************************************************
 *
 *  Main
 *
 *****************************************************************************************/

int main(int argc, char *argv[])
{ int           nass, narg, kmer;
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
          case 'h':
            if (argv[i][1] == 'H')
              DO_TABLE = 0;
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
            else if (*eptr == 0)
              { HIST_LOW = 1;
                HIST_HGH = 0x7fff;
                break;
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
  
  { int c;

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

    gc_setup(kmer);
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
#ifndef DEBUG_THREADS
    pthread_t threads[NTHREADS];
#endif
    TP        parm[NTHREADS];
    FILE    **out[NTHREADS];
    int64    *prefx[nass];
    int       ixlen = 0;
    int       pivot;
    char     *seq;
    uint8    *ent;
    int       t, a, i;
    int64     p;

    if (DO_TABLE)
      { ixlen = 0x1000000;
        prefx[0] = Malloc(sizeof(int64)*ixlen*nass,"Allocating prefix tables");
        bzero(prefx[0],sizeof(int64)*ixlen*nass);
        for (a = 1; a < nass; a++)
          prefx[a] = prefx[a-1] + ixlen;

        out[0] = Malloc(sizeof(FILE *)*nass*NTHREADS,"Allocating thread working memory");
        for (t = 1; t < NTHREADS; t++)
          out[t] = out[t-1] + nass;
        for (t = 0; t < NTHREADS; t++)
          for (a = 0; a < nass; a++)
            out[t][a] = fopen(Catenate(A[a]->path,"/.",A[a]->root,
                                        Numbered_Suffix(".ktab.",t+1,"")),"w");
      }

    for (a = 0; a < narg; a++)
      { range[0][a] = 0;
        range[NTHREADS][a] = S[a]->nels;
      }

    pivot = 0;
    for (a = 1; a < narg; a++)
      if (S[a]->nels > S[pivot]->nels)
        pivot = a;

    seq = Current_Kmer(S[0],NULL);
    ent = Current_Entry(S[0],NULL);
    for (t = 1; t < NTHREADS; t++)
      { p = (S[pivot]->nels*t)/NTHREADS; 
        GoTo_Kmer_Index(S[pivot],p);

#ifdef DEBUG
        printf("\n %lld:",p);
        if (p < S[pivot]->nels)
          printf(" %s\n",Current_Kmer(S[pivot],seq));
        else
          printf(" EOT\n");
#endif

        if (p >= S[pivot]->nels)
          for (a = 0; a < narg; a++)
            range[t][a] = S[a]->nels;
        else
          { ent = Current_Entry(S[pivot],ent);                //  Break at prefix boundaries
            for (i = IB_OUT; i < S[0]->kbyte; i++)
              ent[i] = 0;
            for (a = 0; a < narg; a++)
              { GoTo_Kmer_Entry(S[a],ent);
#ifdef DEBUG
                printf(" %lld:",S[a]->cidx);
                if (S[a]->cidx < S[a]->nels)
                  printf("  %s\n",Current_Kmer(S[a],seq));
                else
                  printf(" EOT\n");
#endif
                range[t][a] = S[a]->cidx;
              }
          }
      }
    free(seq);

    for (t = 0; t < NTHREADS; t++)
      { parm[t].tid   = t;
        parm[t].S     = S;
        parm[t].narg  = narg;
        parm[t].A     = A;
        parm[t].nass  = nass;
        parm[t].begs  = range[t];
        parm[t].ends  = range[t+1];
        if (DO_TABLE)
          { parm[t].prefx = prefx;
            parm[t].out   = out[t];
          }
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

    if (DO_TABLE)
      { int minval;
        int three = 3;
        int mins[narg];

        for (a = 0; a < narg; a++)
          mins[a] = S[a]->minval;

        for (a = 0; a < nass; a++)
          { FILE  *f   = fopen(Catenate(A[a]->path,"/",A[a]->root,".ktab"),"w");
            int64 *prf = prefx[a];

            minval = eval_minimums(A[a]->expr,mins);

            fwrite(&kmer,sizeof(int),1,f);
            fwrite(&NTHREADS,sizeof(int),1,f);
            fwrite(&minval,sizeof(int),1,f);
            fwrite(&three,sizeof(int),1,f);

            for (i = 1; i < ixlen; i++)
              prf[i] += prf[i-1];

            fwrite(prf,sizeof(int64),ixlen,f);
            fclose(f);
          }

        free(out[0]);
        free(prefx[0]);
      }

    if (HIST_LOW > 0)
      { int64 *hist0, *histt;
        FILE  *f;

        for (i = 0; i < nass; i++)
          { hist0 = parm[0].hist[i];
            for (t = 1; t < NTHREADS; t++)
              { histt = parm[t].hist[i];
                for (a = HIST_LOW; a <= HIST_HGH+2; a++)
                  hist0[a] += histt[a];
              }

            f = fopen(Catenate(A[i]->path,"/",A[i]->root,".hist"),"w");
            fwrite(&(S[0]->kmer),sizeof(int),1,f);
            fwrite(&HIST_LOW,sizeof(int),1,f);
            fwrite(&HIST_HGH,sizeof(int),1,f);
            fwrite(hist0+HIST_HGH+1,sizeof(int64),1,f);
            fwrite(hist0+HIST_HGH+2,sizeof(int64),1,f);
            fwrite(hist0+HIST_LOW,sizeof(int64),(HIST_HGH-HIST_LOW)+1,f);
            fclose(f);
          }

        for (t = 0; t < NTHREADS; t++)
          { free(parm[t].hist[0] + HIST_LOW);
            free(parm[t].hist);
          }
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
