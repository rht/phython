:Begin:
:Function: readform_file
:Pattern: ReadForm[filename_String]
:Arguments: {filename}
:ArgumentTypes: {Manual}
:ReturnType: Manual
:End:

:Begin:
:Function: readform_exec
:Pattern: ReadForm[formcmd_String, incpath_String, filename_String]
:Arguments: {formcmd, incpath, filename}
:ArgumentTypes: {Manual}
:ReturnType: Manual
:End:

:Evaluate: _ReadForm := (Message[ReadForm::syntax]; Abort[])

:Begin:
:Function: readformclear
:Pattern: ReadFormClear[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Begin:
:Function: readformdebug
:Pattern: ReadFormDebug[debug_Integer, filename_:""]
:Arguments: {debug, filename}
:ArgumentTypes: {Integer, Manual}
:ReturnType: Manual
:End:

:Evaluate: ReadForm::syntax = "Bad syntax."

:Evaluate: ReadForm::noopen = "Cannot open `1`."

:Evaluate: ReadForm::nooutput =
  "Something went wrong, there was no output from FORM."

:Evaluate: ReadForm::toomany =
  "Too many expressions.  Increase MAXEXPR in ReadForm.tm."

:Evaluate: ReadForm::formerror = "`1`"


/*
	ReadForm.tm
		reads FORM output back into Mathematica
		this file is part of FormCalc
		last modified 23 Dec 11 th

Note: FORM code must have
	1. #- (no listing),
	2. off stats,
	3. should produce output with print (not print +s).

Debug:
	bit 0   = stderr listing of file output
	bit 1   = oob communication
	bit 1+2 = + FORM -> Mma verbose
	bit 1+3 = + Mma -> FORM raw verbose
	bit 1+4 = + Mma -> FORM verbose
	bit 5   = final result transfer
	bit 6   = internal chains before OrderChain
	bit 7   = internal chains after OrderChain
*/

#include "mathlink.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <signal.h>
#include <sys/wait.h>

#define MAXEXPR 5000
#define TERMBUF 500000
#define STRINGSIZE 32767

#ifndef WITHEXTERNALCHANNEL
#define WITHEXTERNALCHANNEL 1
#endif

#define Abort(s) abort1(s, __LINE__)
#define abort1(s, line) abort2(s, line)
#define abort2(s, line) { perror(s " " __FILE__ "(" #line ")"); exit(1); }

#define Die(p) if( (p) == NULL ) Abort("malloc")

#define DEBUG "\e[31m"
#define RESET "\e[0m\n"


/*

A term in the FORM output is organized into the TERM structure
in the following way:

 ____4_____     __3___     __2___     ___0____     _____1_____
/          \   /      \   /      \   /        \   /           \
SumOver(...) * Mat(...) * Den(...) * paveM(...) * ..... * (...)

Hierarchy of collecting:
4. SumOver
3. Mat
2. Den
1. [coefficient]
0. paveM

*/

#define LEVEL_PAVE 0
#define LEVEL_COEFF 1
#define LEVEL_DEN 2
#define LEVEL_MAT 3
#define LEVEL_SUMOVER 4
#define NLEVELS 5

typedef const int cint;
typedef unsigned char byte;
typedef byte *string;
typedef MLCONST byte *cstring;

typedef struct {
  cstring name;
  int level;
} FUN;

static const FUN funtab[] = {
  {(cstring)"SumOver", LEVEL_SUMOVER},
  {(cstring)"Mat",     LEVEL_MAT},
  {(cstring)"Den",     LEVEL_DEN},
  {(cstring)"paveM",   LEVEL_PAVE},
  {(cstring)"A0i",     LEVEL_PAVE},
  {(cstring)"B0i",     LEVEL_PAVE},
  {(cstring)"C0i",     LEVEL_PAVE},
  {(cstring)"D0i",     LEVEL_PAVE},
  {(cstring)"E0i",     LEVEL_PAVE},
  {(cstring)"F0i",     LEVEL_PAVE},
  {(cstring)"Acut",    LEVEL_PAVE},
  {(cstring)"Bcut",    LEVEL_PAVE},
  {(cstring)"Ccut",    LEVEL_PAVE},
  {(cstring)"Dcut",    LEVEL_PAVE},
  {(cstring)"Ecut",    LEVEL_PAVE},
  {(cstring)"Fcut",    LEVEL_PAVE},
  {(cstring)"A0",      LEVEL_PAVE},
  {(cstring)"A00",     LEVEL_PAVE},
  {(cstring)"B0",      LEVEL_PAVE},
  {(cstring)"B1",      LEVEL_PAVE},
  {(cstring)"B00",     LEVEL_PAVE},
  {(cstring)"B11",     LEVEL_PAVE},
  {(cstring)"B001",    LEVEL_PAVE},
  {(cstring)"B111",    LEVEL_PAVE},
  {(cstring)"DB0",     LEVEL_PAVE},
  {(cstring)"DB1",     LEVEL_PAVE},
  {(cstring)"DB00",    LEVEL_PAVE},
  {(cstring)"C0",      LEVEL_PAVE},
  {(cstring)"D0",      LEVEL_PAVE},
  {(cstring)"E0",      LEVEL_PAVE},
  {(cstring)"F0",      LEVEL_PAVE}
};

typedef struct term {
  struct term *last;
  string f[NLEVELS];
  int nterms[NLEVELS], coll;
  byte expr[];
} TERM;

typedef struct btree {
  struct btree *lt, *gt;
  string abbr;
  byte expr[];
} BTREE;

static byte zero[] = "";
static BTREE *root = NULL;
static int debug = 0;
static FILE *stddeb;

/******************************************************************/

#define Strlen(s) strlen((const char *)s)
#define Strchr(s, c) (string)strchr((const char *)s, (char)c)
#define Strstr(s1, s2) (string)strstr((const char *)s1, (const char *)s2)
#define Strcmp(s1, s2) strcmp((const char *)s1, (const char *)s2)
#define Strncpy(s1, s2, n) strncpy((char *)s1, (const char *)s2, n)
#define MLPutStr(mlp, s) MLPutByteString(mlp, s, Strlen(s))

#if 1
static inline long Strtol(cstring s, string *e, int base) {
  return strtol((const char *)s, (char **)e, base);
}
#else
#define Strtol(s, e, base) strtol((const char *)s, (char **)e, base)
#endif

/******************************************************************/

static void PrintPointer(TERM *tp)
{
#if 1
	/* for tough cases w/segfault on access */
  fprintf(stddeb, "address: %p\n", tp);
  fprintf(stddeb, "PAVE:    %s\n", tp->f[LEVEL_PAVE]);
  fprintf(stddeb, "COEFF:   %s\n", tp->f[LEVEL_COEFF]);
  fprintf(stddeb, "DEN:     %s\n", tp->f[LEVEL_DEN]);
  fprintf(stddeb, "MAT:     %s\n", tp->f[LEVEL_MAT]);
  fprintf(stddeb, "SUMOVER: %s\n", tp->f[LEVEL_SUMOVER]);
  fprintf(stddeb, "coll:    %d\n\n", tp->coll);
#else
  fprintf(stddeb,
    "address: %p\n"
    "PAVE:    %s\n"
    "COEFF:   %s\n"
    "DEN:     %s\n"
    "MAT:     %s\n"
    "SUMOVER: %s\n"
    "coll:    %d\n\n",
    tp,
    tp->f[LEVEL_PAVE],
    tp->f[LEVEL_COEFF],
    tp->f[LEVEL_DEN],
    tp->f[LEVEL_MAT],
    tp->f[LEVEL_SUMOVER],
    tp->coll);
#endif
}

static void PrintChain(TERM *tp, const char *info)
{
  int n = 0;
  while( tp ) {
    fprintf(stddeb, "\n%s term %d:\n", info, ++n);
    PrintPointer(tp);
    tp = tp->last;
  }
}

/******************************************************************/

static inline string MLString(MLINK mlp)
{
  cstring s;
  string d;
  int n;

  if( MLGetByteString(mlp, &s, &n, ' ') == 0 ) {
    MLClearError(mlp);
    MLNewPacket(mlp);
    return NULL;
  }
  d = malloc(n + 1);
  if( d ) {
    memcpy(d, s, n);
    d[n] = 0;
  }
  MLReleaseByteString(mlp, s, n);
  return d;
}

/******************************************************************/

static inline void MLSendPacket(MLINK mlp)
{
  MLEndPacket(mlp);
  while( MLNextPacket(mlp) != RETURNPKT )
    MLNewPacket(mlp);
}

/******************************************************************/

static inline void MLEmitMessage(MLINK mlp, cstring tag, cstring arg)
{
  MLPutFunction(mlp, "EvaluatePacket", 1);

  MLPutFunction(mlp, "Message", (arg) ? 2 : 1);
  MLPutFunction(mlp, "MessageName", 2);
  MLPutSymbol(mlp, "ReadForm");
  MLPutStr(mlp, tag);
  if( arg ) MLPutStr(mlp, arg);
  MLSendPacket(mlp);
  MLNewPacket(mlp);	/* discard returned Null */
}

/******************************************************************/

static inline int MLPutExpr(MLINK mlp, TERM *tp, int lev)
{
  static const char *levname[] = {"pave", "coeff", "den", "mat", "sum"};
  if( debug & 32 )
    fprintf(stddeb, "toexpr %s |%s|\n", levname[lev], tp->f[lev]);
  MLPutFunction(mlp, "ToExpression", 1);
  return MLPutStr(mlp, tp->f[lev]);
}

/******************************************************************/

static inline void InsertLor(string s)
{
  s -= 3;
  do s[3] = s[0]; while( *--s < 'A' );
  memcpy(s, "Lor[", 4);
}

/******************************************************************/

static string GetAbbr(string expr)
{
  BTREE *lp, **node = &root;
  cstring abbr;
  int exprlen, abbrlen;

  while( (lp = *node) ) {
    cint t = Strcmp(expr, lp->expr);
    if( t == 0 ) return lp->abbr;
    node = (t < 0) ? &lp->lt : &lp->gt;
  }

  MLPutFunction(stdlink, "EvaluatePacket", 1);
  MLPutFunction(stdlink, "FormEval", 1);
  exprlen = Strlen(expr);
  MLPutByteString(stdlink, expr, exprlen);
  MLSendPacket(stdlink);

  if( MLGetByteString(stdlink, &abbr, &abbrlen, ' ') == 0 ) {
	/* returned expr not a string (e.g. FormEval not defined) */
    MLClearError(stdlink);
    MLNewPacket(stdlink);
    return expr;
  }

  Die(lp = malloc(exprlen + abbrlen + 2 + sizeof *lp));
  *node = lp;
  lp->lt = lp->gt = NULL;
  memcpy(lp->expr, expr, exprlen);
  lp->expr[exprlen] = 0;
  lp->abbr = lp->expr + exprlen + 1;
  memcpy(lp->abbr, abbr, abbrlen);
  lp->abbr[abbrlen] = 0;

  MLReleaseByteString(stdlink, abbr, abbrlen);

  return lp->abbr;
}

/******************************************************************/

static inline void Shift(TERM *tp, cstring begin, cstring end, cstring new)
{
  if( new != begin ) {
    string *f;
    for( f = tp->f; f < &tp->f[NLEVELS]; ++f )
      if( (unsigned)(*f - begin) < (unsigned)(end - begin) )
        *f += new - begin;
  }
}

/******************************************************************/

static TERM *Resize(TERM *tp, cstring end, cint newsize)
{
  TERM *new;
  Die(new = realloc(tp, newsize));
  Shift(new, tp->expr, end, new->expr);
  return new;
}

/******************************************************************/

static inline void MoveToEnd(TERM *tp, cint lev, cstring end)
{
  string s = tp->f[lev];
  cint len = Strlen(s);
  cstring begin = s + len + 1;
  cint rest = end - begin;
  string tmp;

  Die(tmp = malloc(len));
  memcpy(tmp, s, len);
  memmove(s, begin, rest);
  s += rest;
  *s++ = 0;
  memcpy(s, tmp, len);
  free(tmp);

  Shift(tp, begin, end, tp->f[lev]);
  tp->f[lev] = s;
}

/******************************************************************/

static inline TERM *FinalizeTerm(TERM *tp, string di, int *maxpavesize)
{
  *di++ = 0;
  if( *tp->f[LEVEL_MAT] )
    tp->f[LEVEL_MAT] = GetAbbr(tp->f[LEVEL_MAT]);
  if( *tp->f[LEVEL_PAVE] )
    *maxpavesize += Strlen(tp->f[LEVEL_PAVE]) + 2;
  return Resize(tp, di, di - (cstring)tp);
}

/******************************************************************/

static string PutFactor(string to, cstring from)
{
  if( *from ) {
    cint len = Strlen(from) + 1;
    memcpy(to, from, len);
    return to + len;
  }
  *to++ = '1';
  *to++ = 0;
  return to;
}

/******************************************************************/

static void CollectPaVe(TERM *termp, cint maxpavesize)
{
  do {
    TERM *old = termp, *tp;
    string s = termp->f[LEVEL_PAVE];

    while( (tp = old->last) ) {
      int lev;
      for( lev = LEVEL_PAVE + 1; lev < NLEVELS; ++lev )
        if( Strcmp(termp->f[lev], tp->f[lev]) != 0 ) {
          old = tp;
          goto loop;
        }
      if( termp->coll == 0 ) {
        Die(termp->f[LEVEL_PAVE] = malloc(maxpavesize));
        s = PutFactor(termp->f[LEVEL_PAVE], s);
        termp->coll = 1;
      }
      s[-1] = '+';
      s = PutFactor(s, tp->f[LEVEL_PAVE]);
      old->last = tp->last;
      free(tp);
loop: ;
    }

    if( termp->coll ) termp->f[LEVEL_PAVE] =
      realloc(termp->f[LEVEL_PAVE], s - termp->f[LEVEL_PAVE]);
  } while( (termp = termp->last) );
}

/******************************************************************/

static TERM *OrderChain(TERM *t1p, cint level)
{
  TERM *old1;
  int *const nterms = &t1p->nterms[level];
  int c = 0;

  do {
    TERM *t2p, *old2, *ini = t1p;
    int c2 = 0;
    ++c;

    do {
      ++c2;
      old1 = t1p;
      t1p = old1->last;
      if( t1p == NULL ) goto next;
    } while( Strcmp(ini->f[level], t1p->f[level]) == 0 );

    t2p = t1p;

    do {
      old2 = t2p;
over:
      t2p = old2->last;
      if( t2p == NULL ) goto next;
    } while( Strcmp(ini->f[level], t2p->f[level]) != 0 );

    old1->last = t2p;
    old1 = t2p;
    old2->last = t2p->last;
    ++c2;
    goto over;

next:
    if( level > LEVEL_COEFF ) {
      old1->last = NULL;
      old1 = OrderChain(ini, level - 1);
    }
    else ini->nterms[level - 1] = c2;
  } while( (old1->last = t1p) );

  *nterms = c;
  return old1;
}

/******************************************************************/

static TERM *Transmit(TERM *tp, int level)
{
  cint nterms = tp->nterms[level];
  int term;

  if( level == LEVEL_SUMOVER ) MLPutFunction(stdlink, "List", nterms);
  else if( nterms > 1 ) MLPutFunction(stdlink, "Plus", nterms);

  for( term = 1; term <= nterms; ++term ) {
    int lev, ntimes = (*tp->f[level] != 0);

    if( debug & 32 )
      fprintf(stddeb, DEBUG "sending %d of %d terms" RESET,
        term, nterms);

    for( lev = level - 1; lev > LEVEL_PAVE; --lev ) {
      ++ntimes;
      if( tp->nterms[lev] > 1 ) goto sendit;
      if( *tp->f[lev] == 0 ) --ntimes;
    }
	/* OrderChain goes down only to LEVEL_COEFF, hence: */
    if( *tp->f[LEVEL_PAVE] ) ++ntimes;

sendit:
    switch( ntimes ) {
    case 0:
      MLPutInteger(stdlink, 1);
      break;

    default:
      MLPutFunction(stdlink, "Times", ntimes);
    case 1:
      if( *tp->f[level] ) MLPutExpr(stdlink, tp, level);
      for( lev = level - 1; lev > LEVEL_PAVE; --lev ) {
        if( tp->nterms[lev] > 1 ) {
          tp = Transmit(tp, lev);
          goto loop;
        }
        if( *tp->f[lev] ) MLPutExpr(stdlink, tp, lev);
      }
      if( *tp->f[LEVEL_PAVE] ) MLPutExpr(stdlink, tp, LEVEL_PAVE);
      break;
    }
    tp = tp->last;
loop: ;
  }

  return tp;
}

/******************************************************************/

static void ReadForm(FILE *file)
{
  TERM *expressions[MAXEXPR], **exprp = expressions;
  TERM *termp = NULL, *tp = NULL;
  byte br[64];
  int inexpr = 0, maxpavesize = 0, b = 0, thislev = 0;
  cstring beg = NULL;
  string delim = NULL, di = NULL;
  enum { nfun = sizeof funtab/sizeof(FUN) };
  int lineno = 0;

  maxpavesize = 0;

  for( ; ; ) {
    byte line[256];
    cstring si;
    byte errmsg[512];
    string errend = errmsg;
    int i, tpsize = sizeof *tp + TERMBUF;

    *errend = 0;

nextline:
    for( ; ; ) {
      string pos;

      *line = 0;
      si = (string)fgets((char *)line, sizeof line, file);
      if( debug & 1 ) fprintf(stddeb, "%06d %c%s",
        ++lineno,
        inexpr ? (b ? '&' : ' ') : '*',
        (const char *)line);

      if( si == NULL || MLAbort ) {
        int nexpr;
        TERM **ep;

        if( !feof(file) ) goto abort;

        if( errend > errmsg ) {
          errend[-1] = 0;	/* discard last \n */
          MLEmitMessage(stdlink, (cstring)"formerror", errmsg);
          goto abort;
        }

        nexpr = exprp - expressions;
        if( nexpr == 0 ) {
          MLEmitMessage(stdlink, (cstring)"nooutput", NULL);
          goto abort;
        }

        /* successful exit */
        MLPutFunction(stdlink, "FormExpr", nexpr);
        for( ep = expressions; ep < exprp; ++ep ) {
          if( debug & 64 ) PrintChain(*ep, "unordered");
          OrderChain(*ep, NLEVELS - 1);
          if( debug & 128 ) PrintChain(*ep, "ordered");
          Transmit(*ep, NLEVELS - 1);
        }
        goto quit;
      }

      if( (pos = Strstr(si, "-->")) ||
          (pos = Strstr(si, "==>")) ||
          (pos = Strstr(si, "===")) ) {
        pos += 4;
        if( Strstr(errmsg, pos) == NULL ) {
          Strncpy(errend, pos, errmsg + sizeof errmsg - errend);
          errend += Strlen(errend);
        }
        continue;
      }

      if( inexpr ) {
        if( *si >= ' ' ) break;  /* catch both \n and \r (on Windows) */
        if( di == tp->expr ) continue;
        termp = FinalizeTerm(termp, di, &maxpavesize);
      }
      else if( (pos = Strchr(si, '=')) == NULL ) continue;

      Die(tp = malloc(tpsize));
      beg = delim = di = tp->expr;
      tp->last = termp;
      termp = tp;
      for( i = 0; i < NLEVELS; ++i ) tp->f[i] = zero;
      tp->f[thislev = LEVEL_COEFF] = di;
      tp->coll = 0;

      if( !inexpr ) {
        inexpr = 1;
        si = pos + 1;
        break;
      }
    }

    if( di > (cstring)tp + tpsize - sizeof line )
      tp = Resize(tp, di, tpsize += TERMBUF);

    while( *si ) {
      byte c = *si++;
      if( c <= ' ' ) continue;

      switch( c ) {
      case '+':
      case '-':
      case '*':
        if( b == 0 ) delim = di + 1;
        break;

      case '(':
        if( b == 0 ) {
          int newlev = LEVEL_COEFF;
          *di = 0;
          for( i = 0; i < nfun; ++i )
            if( Strcmp(delim, funtab[i].name) == 0 ) {
              newlev = funtab[i].level;
              break;
            }

          if( thislev != newlev ) {
            if( *termp->f[newlev] )
              MoveToEnd(termp, newlev, delim - 1);
            else {
              if( delim > beg ) delim[-1] = 0;
              termp->f[newlev] = delim;
            }
            thislev = newlev;
          }
        }

        if( di == beg || Strchr("+-*/^,([", di[-1]) ) br[b++] = ')';
        else c = '[', br[b++] = ']';
        break;

      case ')':
        if( b > 0 ) c = br[--b];
        break;

      case '[':
        *di++ = '\\';
        break;

      case ';':
        *exprp++ = termp = FinalizeTerm(termp, di, &maxpavesize);
        if( exprp >= expressions + MAXEXPR ) {
          MLEmitMessage(stdlink, (cstring)"toomany", NULL);
          goto abort;
        }
        if( maxpavesize ) CollectPaVe(termp, maxpavesize);
        termp = NULL;
        maxpavesize = inexpr = 0;
        goto nextline;

      case '_':
        *di++ = c = '$';
        break;

      case '?':
        InsertLor(di++);
        c = ']';
        break;
      }

      *di++ = c;
    }
  }

abort:
  MLPutFunction(stdlink, "Abort", 0);

quit:
  MLEndPacket(stdlink);

  while( exprp > expressions ) {
    TERM *last;
    for( tp = *--exprp; tp; tp = last ) {
      if( tp->coll ) free(tp->f[LEVEL_PAVE]);
      last = tp->last;
      free(tp);
    }
  }
}

/******************************************************************/

static void readform_file(void)
{
  string filename = MLString(stdlink);
  FILE *file = fopen((const char *)filename, "r");

  if( file ) {
    ReadForm(file);
    fclose(file);
  }
  else {
    MLEmitMessage(stdlink, (cstring)"noopen", filename);
    MLPutFunction(stdlink, "Abort", 0);
    MLEndPacket(stdlink);
  }

  free(filename);
}

/******************************************************************/

static inline int writeall(cint h, cstring buf, long n)
{
  long w = 0;
  while( (w = write(h, buf, n -= w)) > 0 );
  if( w < 0 ) close(h);
  return w;
}

/******************************************************************/

static int ToMma(cint hw, string expr)
{
  int b = 0, decl = 256, verb = 0;
  cint exprlen = Strlen(expr);
  byte br[64], c;
  string result, r, s;

  if( debug & 2 ) {
    fprintf(stddeb, DEBUG "to mma (%lu bytes)" RESET,
      (unsigned long)exprlen);
    if( debug & 4 ) fprintf(stddeb, "4 |%s|\n", expr);
  }

  MLPutFunction(stdlink, "EvaluatePacket", 1);
  MLPutFunction(stdlink, "FormEvalDecl", 1);
  MLPutByteString(stdlink, expr, exprlen);
  MLSendPacket(stdlink);

  if( (result = MLString(stdlink)) == NULL ) return 0;

  if( debug & 2 ) {
    fprintf(stddeb, DEBUG "from mma raw (%lu bytes)" RESET,
      (unsigned long)Strlen(result));
    if( debug & 8 ) fprintf(stddeb, "8 |%s|\n", result);
  }

  for( r = result, s = expr; (c = *r++); ) {
    if( (c | decl) <= ' ' ) continue;
    switch( c | verb ) {
    case '$':
      if( *r == c ) ++r, c = '_';
      break;
    case '{':
      c = '(', br[b++] = ')';
      if( b == 1 ) goto send;
      *s++ = 'L';
      *s++ = 'i';
      *s++ = 's';
      *s++ = 't';
      break;
    case '[':
      if( r[-2] == '\\' ) verb = 256, br[b++] = ']';
      else c = '(', br[b++] = ')';
      break;
    case '(':
      br[b++] = ')';
      break;
    case ']' + 256:
      verb = 0;
    case ']':
    case ')':
    case '}':
      if( b > 0 ) c = br[--b];
      break;
    case '\\':
      if( *r != '0' ) continue;
      c = Strtol(r, &r, 8);
      break;
    case ',':
      if( decl | (b - 1) ) break;
      *s++ = '\n';
send:
      if( debug & 2 ) {
        fprintf(stddeb, DEBUG "from mma (%lu bytes)" RESET,
          (unsigned long)(s - expr));
        if( debug & 16 ) {
          *s = 0;
          fprintf(stddeb, "16 |%s|\n", expr);
        }
      }
      *s++ = '\n';
      if( writeall(hw, expr, s - expr) < 0 ) {
        free(result);
        return 0;
      }
      s = expr;
      decl = 0;
      continue;
    }
    *s++ = c;
  }

  sync();  /* another sync needed for Mac OS, again unclear why */

  free(result);
  return 1;
}

/******************************************************************/

static void *ExtIO(void *h)
{
  cint hw = ((int *)h)[3];
  cint hr = ((int *)h)[4];
  string expr;
  byte br[64];
  enum { blocksize = 40960, linesize = 512, ahead = 32 };
  int size = 2*blocksize, verb = 0, b, w, n;

  if( (n = read(hr, br, sizeof br) - 1) < 0 ||
      writeall(hw, br,
        n + sprintf((char *)br + n, ",%d\n", getpid())) < 0 )
    pthread_exit(NULL);

  Die(expr = malloc(size));

loop:
  w = 1;
  b = 0;

  for( ; ; ) {
    long r = w + ahead;
    string s;

    if( r + linesize > size ) {
      size += blocksize;
      Die(expr = realloc(expr, size));
    }

    sync();  /* needed for Mac OS, unclear why */

    n = read(hr, s = expr + r, size - r);
    if( n <= 0 ) {
      free(expr);
      pthread_exit(NULL);
    }

    do {
      byte c = *s++;
      if( c <= ' ' ) continue;
      switch( c | verb ) {
      case '#':
        expr[w++] = '0';
        expr[w++] = '}';
        expr[w] = 0;
        *expr = '{';
        if( ToMma(hw, expr) ) goto loop;
        free(expr);
        pthread_exit(NULL);
      case '?':
        InsertLor(expr + w++);
        c = ']';
        break;
      case '_':
        expr[w++] = c = '$';
        break;
      case '(':
        if( Strchr("+-*/^,([{", expr[w-1]) ) br[b++] = ')';
        else c = '[', br[b++] = ']';
        break;
      case ')':
        if( b > 0 ) c = br[--b];
        break;
      case '[':
        expr[w++] = '\\';
        verb = 256;
        break;
      case ']' + 256:
        verb = 0;
        break;
      }
      expr[w++] = c;
    } while( --n );
  }
}

/******************************************************************/

static void readform_exec(void)
{
  int hh[10] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
  int *h, i;
  pid_t pid = -1;
  FILE *file;
  string formcmd = MLString(stdlink);
  string incpath = MLString(stdlink);
  string filename = MLString(stdlink);
  string argv[] = {formcmd, (string)"-p", incpath,
    filename, NULL, filename, NULL};
#if WITHEXTERNALCHANNEL
  pthread_t tid;
  void *tj = NULL;
  byte arg[32];
#endif

  /* make sure we don't overlap with 0, 1, 2 */
  i = -2;
  do {
    h = &hh[i += 2];
    if( pipe(h) == -1 ) goto abort;
  } while( h[1] <= 2 );

#if WITHEXTERNALCHANNEL
  if( pipe(&h[2]) != -1 &&
      pipe(&h[4]) != -1 &&
      pthread_create(&tid, NULL, ExtIO, h) == 0 ) {
    tj = (void *)1;
    sprintf((char *)arg, "%d,%d", h[2], h[5]);
    argv[3] = (string)"-pipe";
    argv[4] = arg;
  }
#endif

  while( i > 0 ) close(hh[--i]);

  signal(SIGCHLD, SIG_IGN);
  pid = fork();
  if( pid == -1 ) goto abort;

  if( pid == 0 ) {
    usleep(1000);
    close(h[0]);
    dup2(h[1], 1);
    close(h[1]);
    if( h[3] != -1 ) close(h[3]);
    if( h[4] != -1 ) close(h[4]);
    exit(execvp((char *)argv[0], (char **)argv));
  }

  close(h[1]);
  if( h[2] != -1 ) close(h[2]);
  if( h[5] != -1 ) close(h[5]);
  h[1] = h[2] = h[5] = -1;

  file = fdopen(h[0], "r");
  if( file ) {
    ReadForm(file);
    fclose(file);
    h[0] = -1;
  }
  else {
abort:
    MLEmitMessage(stdlink, (cstring)"noopen", formcmd);
    MLPutFunction(stdlink, "Abort", 0);
    MLEndPacket(stdlink);
  }

  if( pid > 0 ) {
    kill(pid, SIGKILL);
    wait(&i);
  }

  for( i = 6; --i >= 0; ) if( h[i] != -1 ) close(h[i]);

#if WITHEXTERNALCHANNEL
  if( tj ) pthread_join(tid, &tj);
#endif

  free(filename);
  free(incpath);
  free(formcmd);
}

/******************************************************************/

static void CutBranch(BTREE *node)
{
  if( node ) {
    CutBranch(node->lt);
    CutBranch(node->gt);
    free(node);
  }
}

/******************************************************************/

static void readformclear(void)
{
  CutBranch(root);
  root = NULL;
  MLPutSymbol(stdlink, "Null");
  MLEndPacket(stdlink);
}

/******************************************************************/

static void readformdebug(cint deb)
{
  cstring filename = MLString(stdlink);
  debug = deb;

  stddeb = stderr;
  if( filename && *filename ) {
    stddeb = fopen((const char *)filename, "w");
    if( stddeb == NULL ) {
      MLEmitMessage(stdlink, (cstring)"noopen", filename);
      MLPutSymbol(stdlink, "$Failed");
      MLEndPacket(stdlink);
      debug = 0;
      return;
    }
    setbuf(stddeb, NULL);
  }

  MLPutSymbol(stdlink, "True");
  MLEndPacket(stdlink);
}

/******************************************************************/

int main(int argc, char **argv)
{
  return MLMain(argc, argv);
}

