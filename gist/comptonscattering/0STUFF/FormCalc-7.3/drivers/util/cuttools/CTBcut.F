* CTBcut.F
* the FormCalc wrapper for the CutTools two-point function
* this file is part of FormCalc
* last modified 3 Aug 11 th


	double complex function CTBcut(rank, num, numeps,
     &    p, m1, m2)
	implicit none
	integer rank
	external num, numeps
	integer*8 p
	double precision m1, m2

	double precision pj(0:3,2), mj(2)
	logical stable, forcemp
	double complex res(0:2), res1
	external CTmpnum

	double precision rootsvalue, muscale
	common /cuttools_para/ rootsvalue, muscale

	pj(0,1) = 0
	pj(1,1) = 0
	pj(2,1) = 0
	pj(3,1) = 0
	call CTMomDecode(p, pj(0,2))

	mj(1) = m1
	mj(2) = m2

	PRINT *, "CTBcut:", rank
	PRINT *, "p1=", pj(0,2)
	PRINT *, "   ", pj(1,2)
	PRINT *, "   ", pj(2,2)
	PRINT *, "   ", pj(3,2)
	PRINT *, "m1=", m1
	PRINT *, "m2=", m2

	call ctsxcut(rootsvalue, muscale, 2, num, CTmpnum,
     &    rank, pj, mj, res, res1, stable, forcemp)
	CTBcut = res(0)
	PRINT *, "CTBcut=", CTBcut
	end

