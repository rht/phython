/*
	psripper.c
		rips comments and whitespace from a PostScript
		prologue file for faster processing
		this file is part of FeynArts
		last modified 7 May 02 th

	Usage: psripper < in.ps > out.ps
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main()
{
  char s[1024], *p, *token, c;
  int instring = 0, linewidth = 0, width;

  while( !feof(stdin) ) {
    *s = 0;
    fgets(s, sizeof(s), stdin);
    if( *s == 0 ) break;
    if( *s == '%' && (*(s + 1) == '%' || *(s + 1) == '!') ) {
      if( linewidth ) {
        putchar('\n');
        linewidth = 0;
      }
      if( strcmp(s + 2, "BeginProlog\n") == 0 ||
          strcmp(s + 2, "BeginSetup\n") == 0 ) putchar('\n');
      printf("%s", s);
      continue;
    }

    token = NULL;

    for( p = s; ; ++p ) {
      switch( *p ) {
      case ' ':
      case '\t':
      case '%':
        if( instring ) break;
      case '\n':
      case 0:
        c = *p;
        if( token ) {
          *p = 0;
          width = strlen(token) + 1;
          if( linewidth ) {
            if( width + linewidth < 75 ) putchar(' ');
            else {
              putchar('\n');
              linewidth = 0;
            }
          }
          linewidth += width;
          printf("%s", token);
          token = NULL;
        }          
        if( c == ' ' || c == '\t' ) continue;
        goto skiprest;
      case ')':
        if( *(p - 1) != '\\' ) instring = 0;
        break;
      case '(':
        if( *(p - 1) != '\\' ) instring = 1;
        break;
      }
      if( !token ) token = p;
    }
skiprest: ;
  }
  if( linewidth ) putchar('\n');
  return 0;
}

