* 1to2.F
* the kinematics for a 1 -> 2 process
* this file is part of FormCalc
* last modified 30 Nov 11 th


* some definitions for xsection.F
#define LEGS 3
#define LEGS_IN 1
#define LEGS_OUT 2
#define MASS_IN MASS1
#define MASS_OUT MASS2 + MASS3
#define TYPES TYPE1, TYPE2, TYPE3

* this defines the order of integrations:
#define XCOSTH1 1
#define MAXVAR XCOSTH1

#include LUMI


************************************************************************

	subroutine KinIni(fail)
	implicit none
	integer fail

#include "xsection.h"

	Lower(XCOSTH1) = 0
	Upper(XCOSTH1) = 1
	if( Lower(SQRTS) .gt. Upper(SQRTS) ) then
	  Step(XCOSTH1) = XCOSTH1STEP
	else
	  Step(XCOSTH1) = 0
	endif

	fail = 0
	end


************************************************************************

	subroutine InitialState
	implicit none

#include "xsection.h"

	RealType p

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
	RealType fac

#include "xsection.h"

	RealType p, ex,ey,ez, msum, minv
	integer key

	p = 0
	ex = 0
	ey = 0
	ez = 1
	minv = sqrtS
	msum = MASS_OUT

	key = 0
#ifdef COSTH1CUT
	CutMin(XCOSTH1) = -(1 - COSTH1CUT)
	CutMax(XCOSTH1) = +(1 - COSTH1CUT)
	key = key + Cut(CUT_COSTH, CUT_MIN + CUT_MAX)
#endif

	call Split(3, Re(MASS3), p, ex,ey,ez, minv, msum,
     &    fac, key, Var(FIXED), Var(XCOSTH1), Var(TRIVIAL))

	call VecSet(2, Re(MASS2), p, ex,ey,ez)
	end

