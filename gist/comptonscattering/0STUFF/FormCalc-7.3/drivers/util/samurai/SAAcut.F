* SAAcut.F
* the FormCalc wrapper for the Samurai one-point function
* this file is part of FormCalc
* last modified 31 Mar 11 th


	double complex function SAAcut(rank, num, numeps, m)
	implicit none
	integer rank
	double complex num, numeps
	external num, numeps
	double precision m

#include "SA.h"

	double precision pj(4,1)
	double complex tot(-2:0)
	integer ncut

	pj(1,1) = 0
	pj(2,1) = 0
	pj(3,1) = 0
	pj(4,1) = 0

	PRINT *, "SAAcut:", rank
	PRINT *, "m1=", m

	tot(-2) = 0
	tot(-1) = m
	tot(0) = m - log(m/mudim)
	SAAcut = tot(epscoeff)*num(ncut, pj, 0D0)

	PRINT *, "SAAcut=", SAAcut
	end

