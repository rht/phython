:Begin:
:Function: readdata
:Pattern: ReadData[filename_String, setno_Integer:1, parahead_:Para, datahead_:Data]
:Arguments: {filename, setno, ToString[parahead], ToString[datahead]}
:ArgumentTypes: {String, Integer, String, String}
:ReturnType: Manual
:End:

:Evaluate: DataRow = List
:Evaluate: ReadData::noopen = "Cannot open `1`."


/*
	ReadData.tm
		reads data files produced by num.F into Mathematica
		this file is part of FormCalc
		last modified 14 Jul 10 th

known shortcomings:
- fixed memory requirements:
  1) at most MAXPARA parameters in each set
  2) at most MAXCOLS columns
  3) at most MAXDATA data lines in each set

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mathlink.h"

#ifndef MLCONST
#define MLCONST
#endif

#define MAXCOLS 32
#define MAXPARA 100
#define MAXDATA 4000
#define MAXEVAL 1000

typedef struct {
  double val;
  char name[32];
} Para;

static Para para[MAXPARA];
static int npara;

static double data[MAXDATA][MAXCOLS];
static int ncols[MAXDATA];
static int ndata;

static int neval;


static inline void waitpacket()
{
  int pkt;

  MLEndPacket(stdlink);

  do {
    pkt = MLNextPacket(stdlink);
    MLNewPacket(stdlink);
  } while( pkt != RETURNPKT );
}


static void transmit(const char *parahead, const char *datahead, int setno)
{
  int i;

  if( ++neval == 1 ) MLPutFunction(stdlink, "EvaluatePacket", 1);

  MLPutFunction(stdlink, "List", (neval > MAXEVAL) ? 2 : 3);

  MLPutFunction(stdlink, "Set", 2);

  MLPutFunction(stdlink, parahead, 1);
  MLPutInteger(stdlink, setno);

  MLPutFunction(stdlink, "List", npara);
  for( i = 0; i < npara; ++i ) {
    Para *p = &para[i];
    char *s;

    for( s = p->name; *s; ++s )
      switch( *s ) {
      case '_':
        *s = '$';
        break;
      case '(':
        *s = '[';
        break;
      case ')':
        *s = ']';
        break;
      }

    MLPutFunction(stdlink, "Rule", 2);
    MLPutFunction(stdlink, "ToExpression", 1);
    MLPutString(stdlink, p->name);
    MLPutDouble(stdlink, p->val);
  }

  MLPutFunction(stdlink, "Set", 2);

  MLPutFunction(stdlink, datahead, 1);
  MLPutInteger(stdlink, setno);

  MLPutFunction(stdlink, "List", ndata);
  for( i = 0; i < ndata; ++i ) {
    const double *d = data[i];
    int j;

    MLPutFunction(stdlink, "DataRow", ncols[i]);
    for( j = 0; j < ncols[i]; ++j ) MLPutDouble(stdlink, d[j]);
  }

  npara = ndata = 0;

  if( neval > MAXEVAL ) {
    waitpacket();
    neval = 0;
  }
}


static void readdata(const char *filename, int setno,
  const char *parahead, const char *datahead)
{
  char line[4096];
  FILE *file = (*filename == '!') ?
    popen(filename + 1, "r") :
    fopen(filename, "r");

  npara = ndata = neval = 0;

  if( file == NULL ) {
    MLPutFunction(stdlink, "Message", 2);
    MLPutFunction(stdlink, "MessageName", 2);
    MLPutSymbol(stdlink, "ReadData");
    MLPutString(stdlink, "noopen");
    MLPutString(stdlink, filename);
    MLPutSymbol(stdlink, "$Failed");
    MLEndPacket(stdlink);
    return;
  }

  while( fgets(line, sizeof line, file) ) {
    char *s;

    if( *line == '#' ) {
      Para *p = &para[npara];
      if( sscanf(line + 1, " %[^= ] =%lg", p->name, &p->val) == 2 )
        npara++;
      continue;
    }

    s = line + strspn(line, " \t");

    if( *s != 0 && *s != '\n' ) {
      double *d = data[ndata];
      int n;

      for( n = 0; n < MAXCOLS; ++n ) {
        char *p = s;
        *d++ = strtod(s, &s);
        if( s == p ) break;
      }

      if( n ) ncols[ndata++] = n;
      continue;
    }

    if( npara || ndata ) transmit(parahead, datahead, setno++);
  }

  ((*filename == '!') ? pclose : fclose)(file);

  if( neval ) {
    MLPutFunction(stdlink, "List", 0);
    waitpacket();
  }

  MLPutInteger(stdlink, setno - 1);
  MLEndPacket(stdlink);
}


int main(int argc, char **argv)
{
  return MLMain(argc, argv);
}

