/*
	setfpu.c
		set FPU to throw an exception (terminate the
		program) rather than deliver NaNs
		necessary only for g77/gfortran
		this file is part of FormCalc
		last modified 12 Dec 08 th
*/


#include <stdlib.h>

#ifdef __GLIBC_PREREQ
#if __GLIBC_PREREQ(2,2)
#define _GNU_SOURCE
#include <fenv.h>
int feenableexcept(int excepts);
#endif
#endif

#if UNDERSCORE
#define setfpu setfpu_
#endif

void setfpu(void)
{
#ifdef _GNU_SOURCE
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
}

