* SAEcut.F
* the FormCalc wrapper for the Samurai five-point function
* this file is part of FormCalc
* last modified 31 Mar 11 th


	double complex function SAEcut(rank, num, numeps,
     &    p1, p2, p3, p4, m1, m2, m3, m4, m5)
	implicit none
	integer rank
	double complex num, numeps
	external num, numeps
	integer*8 p1, p2, p3, p4
	double precision m1, m2, m3, m4, m5

#include "SA.h"

	double precision pj(5,4), mj(5)
	double complex tot(-2:0), totr
	logical ok

	pj(1,1) = 0
	pj(1,2) = 0
	pj(1,3) = 0
	pj(1,4) = 0
	call SAMomDecode(p1, pj(2,1), 5)
	call SAMomDecode(p2, pj(3,1), 5)
	call SAMomDecode(p3, pj(4,1), 5)
	call SAMomDecode(p4, pj(5,1), 5)

	mj(1) = m1
	mj(2) = m2
	mj(3) = m3
	mj(4) = m4
	mj(5) = m5

	call samurai(num, tot, totr, pj, mj, 5, rank, 1, mudim, ok)
	SAEcut = tot(epscoeff)
	end

