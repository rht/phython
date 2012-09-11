* SADcut.F
* the FormCalc wrapper for the Samurai four-point function
* this file is part of FormCalc
* last modified 31 Mar 11 th


	double complex function SADcut(rank, num, numeps,
     &    p1, p2, p3, m1, m2, m3, m4)
	implicit none
	integer rank
	double complex num, numeps
	external num, numeps
	integer*8 p1, p2, p3
	double precision m1, m2, m3, m4

#include "SA.h"

	double precision pj(4,4), mj(4)
	double complex tot(-2:0), totr
	logical ok

	pj(1,1) = 0
	pj(1,2) = 0
	pj(1,3) = 0
	pj(1,4) = 0
	call SAMomDecode(p1, pj(2,1), 4)
	call SAMomDecode(p2, pj(3,1), 4)
	call SAMomDecode(p3, pj(4,1), 4)

	mj(1) = m1
	mj(2) = m2
	mj(3) = m3
	mj(4) = m4

	PRINT *, "SADcut:", rank
	PRINT *, "p1=", pj(2,4)
	PRINT *, "   ", pj(2,1)
	PRINT *, "   ", pj(2,2)
	PRINT *, "   ", pj(2,3)
	PRINT *, "p2=", pj(3,4)
	PRINT *, "   ", pj(3,1)
	PRINT *, "   ", pj(3,2)
	PRINT *, "   ", pj(3,3)
	PRINT *, "p3=", pj(4,4)
	PRINT *, "   ", pj(4,1)
	PRINT *, "   ", pj(4,2)
	PRINT *, "   ", pj(4,3)
	PRINT *, "m1=", m1
	PRINT *, "m2=", m2
	PRINT *, "m3=", m3
	PRINT *, "m4=", m4

	call samurai(num, tot, totr, pj, mj, 4, rank, 1, mudim, ok)
	SADcut = tot(epscoeff)
	PRINT *, "SADcut=", SADcut
	end

