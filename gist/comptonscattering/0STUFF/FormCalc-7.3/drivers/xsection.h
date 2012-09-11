* xsection.h
* common blocks for xsection.F
* this file is part of FormCalc
* last modified 30 Nov 11 th


#include "decl.h"

#ifndef SQRTS
#define SQRTS 0
#define FIXED MAXVAR+1
#define TRIVIAL MAXVAR+2
#define Var(v) var(1,v)
#define Show(v) var(2,v)
#define Lower(v) var(3,v)
#define Upper(v) var(4,v)
#define Step(v) var(5,v)
#define CutMin(v) var(6,v)
#define CutMax(v) var(7,v)

#define BIT_RESET 0
#define BIT_LOOP 1
#endif

	RealType var(8, MINVAR:TRIVIAL)
	RealType avgfac, sqrtS
	RealType mass_in, mass_out, threshold, scale
	integer*8 helicities
	integer sqrtSinvalid, flags

	common /xsection/ var, avgfac, sqrtS,
     &    mass_in, mass_out, threshold, scale,
     &    helicities, sqrtSinvalid, flags

