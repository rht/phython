* 1to2.F
* the kinematics for a 1 -> 2 process
* this file is part of FormCalc
* last modified 22 May 08 th


* some definitions for xsection.F
#define LEGS 3
#define LEGS_IN 1
#define LEGS_OUT 2
#define MASS_IN MASS1
#define MASS_OUT MASS2 + MASS3
#define TYPES TYPE1, TYPE2, TYPE3

* this defines the order of integrations:
#define XCOSTH 1
#define MAXVAR XCOSTH

#include LUMI


************************************************************************

	subroutine KinIni(fail)
	implicit none
	integer fail

#include "xsection.h"

	Lower(XCOSTH) = 0
	Upper(XCOSTH) = 1
	if( Lower(SQRTS) .gt. Upper(SQRTS) ) then
	  Step(XCOSTH) = XCOSTHSTEP
	else
	  Step(XCOSTH) = 0
	endif

	fail = 0
	end


************************************************************************

	subroutine InitialState
	implicit none

#include "xsection.h"

	double precision p

#ifdef FORCE_ONSHELL
	p = 0
#else
	p = dim(sqrtS, mass_in)
	if( p .gt. 1D-6 ) then
	  INFO "Computing off-shell decay"
	  INFO "Warning: Results correct only if decaying particle"
	  INFO "was chosen off-shell in CalcFeynAmp"
	endif
	p = sqrt(p*(sqrtS + mass_in))
#endif

	call VecSet(1, mass_in, p, 0D0, 0D0, 1D0)
	end


************************************************************************

	subroutine FinalState(fac)
	implicit none
	double precision fac

#include "xsection.h"

	double precision p, ex,ey,ez, msum, minv
	integer key

	p = 0
	ex = 0
	ey = 0
	ez = 1
	minv = sqrtS
	msum = MASS_OUT

	key = 0
#ifdef COSTHCUT
	CutMin(XCOSTH) = -(1 - COSTHCUT)
	CutMax(XCOSTH) = +(1 - COSTHCUT)
	key = key + Cut(CUT_COSTH, CUT_MIN + CUT_MAX)
#endif

	call Split(2, DBLE(MASS2), p, ex,ey,ez, minv, msum,
     &    fac, key, Var(FIXED), Var(XCOSTH), Var(TRIVIAL))

	call VecSet(3, DBLE(MASS3), p, ex,ey,ez)
	end

