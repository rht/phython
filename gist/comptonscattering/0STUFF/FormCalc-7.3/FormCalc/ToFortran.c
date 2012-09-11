/*
	ToFortran.c
		post-processes Mathematica's FortranForm output
		- replaces all real constants by Fortran-style
		  double precision numbers (1.234D0),
		- removes " and indents lines not starting with #
		  (i.e. does not touch preprocessor statements)
		- replaces the continuation character (if any) by &
		this file is part of FormCalc
		last modified 3 Dec 10 th
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main()
{
  static const char signdigits[] = "+-0123456789";
  static const char *digits = signdigits + 2;
  char s[200], next[200];

  for( *next = 0;
       (*next) ? strcpy(s, next), *next = 0, (char *)1 :
                 fgets(s, sizeof s, stdin);
       puts(s) ) {
    char *si, *di, *pos;
    int n = strlen(s);
    char *eol = s + n - 1;
    if( *eol != '\n' ) *++eol = '\n';
    if( eol[-1] == '\\' ) {
      if( fgets(next, sizeof next, stdin) == NULL ) break;
      di = next + 6 + strspn(next + 6, " \t");
      n = strcspn(di, " */()");
      memcpy(--eol, di, n);
      eol += n;
      memmove(di, di + n, strlen(di + n) + 1);
    }
    *eol = 0;

    if( *s == '*' ) continue;

    si = s;
    while( (pos = strpbrk(si, digits)) ) {
      char term;

      si = pos + strspn(pos, digits);
      if( pos[-1] >= 'A' ) continue;  /* belongs to variable name */

      term = *si;
      if( term == '.' ) {
        if( *++si >= 'A' ) continue;
        si += strspn(si, digits);
      }

      if( (*si++ & 0xde) != 'D' ) {
        if( term != '.' ) continue;  /* is an integer */
        for( di = eol += 2; di > si; --di ) *di = di[-2];
        *si = '0';
      }
      si[-1] = 'D';
      si += strspn(si, signdigits);
    }

    if( *s != '#' ) {
      if( strncmp(s, "     ", 5) == 0 && s[5] != ' ' ) s[5] = '&';
      for( si = di = s; ; ++si ) {
        if( *si != '"' ) *di++ = *si;
        if( *si == 0 ) break;
      }
      if( *s >= 'A' ) putchar('\t');
    }
  } /* eof */

  return 0;
}

