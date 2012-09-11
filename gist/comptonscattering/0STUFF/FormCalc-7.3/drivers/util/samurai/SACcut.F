* SACcut.F
* the FormCalc wrapper for the Samurai three-point function
* this file is part of FormCalc
* last modified 31 Mar 11 th


	double complex function SACcut(rank, num, numeps,
     &    p1, p2, m1, m2, m3)
	implicit none
	integer rank
	double complex num, numeps
	external num, numeps
	integer*8 p1, p2
	double precision m1, m2, m3

#include "SA.h"

	double precision pj(3,4), mj(3)
	double complex tot(-2:0), totr
	logical ok

	pj(1,1) = 0
	pj(1,2) = 0
	pj(1,3) = 0
	pj(1,4) = 0
	call SAMomDecode(p1, pj(2,1), 3)
	call SAMomDecode(p2, pj(3,1), 3)

	mj(1) = m1
	mj(2) = m2
	mj(3) = m3

	PRINT *, "SACcut:", rank
	PRINT *, "p1=", pj(2,4)
	PRINT *, "   ", pj(2,1)
	PRINT *, "   ", pj(2,2)
	PRINT *, "   ", pj(2,3)
	PRINT *, "p2=", pj(3,4)
	PRINT *, "   ", pj(3,1)
	PRINT *, "   ", pj(3,2)
	PRINT *, "   ", pj(3,3)
	PRINT *, "m1=", m1
	PRINT *, "m2=", m2
	PRINT *, "m3=", m3

	call samurai(num, tot, totr, pj, mj, 3, rank, 1, mudim, ok)
	SACcut = tot(epscoeff)
	PRINT *, "SACcut=", SACcut
	end

