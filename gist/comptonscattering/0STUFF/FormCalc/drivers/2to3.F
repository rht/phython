* 2to3.F
* the kinematics for a 2 -> 3 process
* this file is part of FormCalc
* last modified 18 Jan 08 th


* some definitions for xsection.F
#define LEGS 5
#define LEGS_IN 2
#define LEGS_OUT 3
#define MASS_IN MASS1 + MASS2
#define MASS_OUT MASS3 + MASS4 + MASS5
#define TYPES TYPE1, TYPE2, TYPE3, TYPE4, TYPE5

* this defines the order of integrations:
#define XMREM5 1
#define XCOSTH5 2
#define XCOSTH4 3
#define XPHI4 4
#define MAXVAR XPHI4

#include LUMI


************************************************************************

	subroutine KinIni(fail)
	implicit none
	integer fail

#include "xsection.h"

	Lower(XMREM5) = 0
	Upper(XMREM5) = 1
	if( Lower(SQRTS) .gt. Upper(SQRTS) ) then
	  Step(XMREM5) = XMREM5STEP
	else
	  Step(XMREM5) = 0
	endif

	Lower(XCOSTH5) = 0
	Upper(XCOSTH5) = 1
	Step(XCOSTH5) = 0

	Lower(XCOSTH4) = 0
	Upper(XCOSTH4) = 1
	Step(XCOSTH4) = 0

	Lower(XPHI4) = 0
	Upper(XPHI4) = 1
	Step(XPHI4) = 0

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
#ifdef E5MIN
	CutMin(XMREM5) = E5MIN
	key = key + Cut(CUT_MREM_E, CUT_MIN)
#endif
#ifdef COSTH5CUT
	CutMin(XCOSTH5) = -(1 - COSTH5CUT)
	CutMax(XCOSTH5) = +(1 - COSTH5CUT)
	key = key + Cut(CUT_COSTH, CUT_MIN + CUT_MAX)
#endif

	call Split(5, DBLE(MASS5), p, ex,ey,ez, minv, msum,
     &    fac, key, Var(XMREM5), Var(XCOSTH5), Var(TRIVIAL))
	Show(XMREM5) = momspec(SPEC_E, 5)
	if( fac .eq. 0 ) return

	call Split(4, DBLE(MASS4), p, ex,ey,ez, minv, msum,
     &    fac, 0, Var(FIXED), Var(XCOSTH4), Var(XPHI4))

	call VecSet(3, DBLE(MASS3), p, ex,ey,ez)
	end

