/*
	logfile.c
		I/O redirection for logfiles
		this file is part of FormCalc
		last modified 26 Feb 08 th
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>

#ifdef UNDERSCORE
#define openlog openlog_
#define closelog closelog_
#endif

static int prevstdout;


int openlog(const char *dir, const int *serial, const int len)
{
  int logfile;
  struct stat filestat;
  char filename[512];

  prevstdout = len;
  do
    if( prevstdout == 0 ) return 0;
  while( dir[--prevstdout] == ' ' );

  memcpy(filename, dir, ++prevstdout);
  filename[prevstdout] = 0;

  if( stat(filename, &filestat) == 0 ) {
    if( !S_ISDIR(filestat.st_mode) ) {
      fprintf(stderr, "%s already exists\n", filename);
      exit(1);
    }
  }
  else if( mkdir(filename, 0755) == -1 ) {
    fprintf(stderr, "Cannot create directory %s\n", dir);
    exit(1);
  }

  sprintf(filename + prevstdout, "/%07d", *serial);
  if( stat(filename, &filestat) == 0 && (filestat.st_mode & 0100) == 0 ) {
    printf("%s already complete\n", filename);
    return 1;
  }

  logfile = creat(filename, 0744);
  if( logfile == -1 ) {
    fprintf(stderr, "Cannot create %s\n", filename);
    exit(1);
  }

  puts(filename);
  fflush(stdout);

  prevstdout = dup(1);
  dup2(logfile, 1);
  close(logfile);

  return 0;
}


void closelog(void)
{
  if( prevstdout ) {
    fchmod(1, 0644);
    dup2(prevstdout, 1);
    close(prevstdout);
  }
}

