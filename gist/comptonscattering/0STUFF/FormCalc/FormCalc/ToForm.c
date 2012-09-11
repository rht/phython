/*
	ToForm.c
		rearranges Mma's InputForm output to yield
		acceptable FORM input
		this file is part of FormCalc
		last modified 14 Jul 10 th
*/

#include <stdio.h>
#include <string.h>

int main()
{
  char in[2048], out[2048];

  while( fgets(in, sizeof in, stdin) ) {
    char *s, *d, closing = ')';

    if( *in == '#' ) {
      fputs(in, stdout);
      continue;
    }

    for( s = in, d = out; *s; ++s ) {
      switch( *s ) {
      case '\\':
        if( s[1] == '\n' ) fgets(s, in + sizeof in - s, stdin);
        break;
      case '"':
        break;
      case '[':
        *d++ = (s[-1] == '\\') ? (closing = ']', '[') : '(';
        break;
      case ']':
        *d++ = closing;
        closing = ')';
        break;
      case '{':
        *d++ = 'L';
        *d++ = 'i';
        *d++ = 's';
        *d++ = 't';
        *d++ = '(';
        break;
      case '}':
        *d++ = ')';
        break;
      case ' ':
        if( s[1] == '.' || s[-1] == '.' ) break;
        if( s - in > 75 ) {
          *d++ = '\n';
          *d = 0;
          fputs(d = out, stdout);
          s = strcpy(in, s);
        }
        *d++ = *s;
        break;
      case '*':
      case '=':
        if( s[-1] == *s ) break;
      default:
        *d++ = *s;
      }
    }

    *d = 0;
    fputs(out, stdout);
  }

  return 0;
}

