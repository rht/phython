* CTAcut.F
* the FormCalc wrapper for the CutTools one-point function
* this file is part of FormCalc
* last modified 3 Aug 11 th


	double complex function CTAcut(rank, num, numeps, m)
	implicit none
	integer rank
	external num, numeps
	double precision m

	double precision pj(0:3,1)
	logical stable, forcemp
	double complex res(0:2), res1
	external CTmpnum

	double precision rootsvalue, muscale
	common /cuttools_para/ rootsvalue, muscale

	pj(0,1) = 0
	pj(1,1) = 0
	pj(2,1) = 0
	pj(3,1) = 0

	PRINT *, "CTAcut:", rank
	PRINT *, "m1=", m

	call ctsxcut(rootsvalue, muscale, 1, num, CTmpnum,
     &    rank, pj, m, res, res1, stable, forcemp)
	CTAcut = res(0)
	PRINT *, "CTAcut=", CTAcut
	end

