* 2to2.F
* the kinematics for a 2 -> 2 process
* this file is part of FormCalc
* last modified 18 Jan 08 th


* some definitions for xsection.F
#define LEGS 4
#define LEGS_IN 2
#define LEGS_OUT 2
#define MASS_IN MASS1 + MASS2
#define MASS_OUT MASS3 + MASS4
#define TYPES TYPE1, TYPE2, TYPE3, TYPE4

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

	p = ThreeMom(sqrtS, DBLE(MASS1), DBLE(MASS2))

	call VecSet(1, DBLE(MASS1), p, 0D0, 0D0, 1D0)
	call VecSet(2, DBLE(MASS2), p, 0D0, 0D0, -1D0)
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

	call Split(3, DBLE(MASS3), p, ex,ey,ez, minv, msum,
     &    fac, key, Var(FIXED), Var(XCOSTH), Var(TRIVIAL))

	call VecSet(4, DBLE(MASS4), p, ex,ey,ez)
	end

