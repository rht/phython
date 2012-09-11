* SABcut.F
* the FormCalc wrapper for the Samurai two-point function
* this file is part of FormCalc
* last modified 31 Mar 11 th


	double complex function SABcut(rank, num, numeps, p, m1, m2)
	implicit none
	integer rank
	double complex num, numeps
	external num, numeps
	integer*8 p
	double precision m1, m2

#include "SA.h"

	double precision pj(2,4), mj(2)
	double complex tot(-2:0), totr
	logical ok

	pj(1,1) = 0
	pj(1,2) = 0
	pj(1,3) = 0
	pj(1,4) = 0
	call SAMomDecode(p, pj(2,1), 2)

	mj(1) = m1
	mj(2) = m2

	PRINT *, "SABcut:", rank
	PRINT *, "p1=", pj(2,4)
	PRINT *, "   ", pj(2,1)
	PRINT *, "   ", pj(2,2)
	PRINT *, "   ", pj(2,3)
	PRINT *, "m1=", m1
	PRINT *, "m2=", m2

	call samurai(num, tot, totr, pj, mj, 2, rank, 1, mudim, ok)
	SABcut = tot(epscoeff)
	PRINT *, "SABcut=", SABcut
	end

