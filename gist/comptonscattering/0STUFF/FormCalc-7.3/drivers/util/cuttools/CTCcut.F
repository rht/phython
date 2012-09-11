* CTCcut.F
* the FormCalc wrapper for the CutTools three-point function
* this file is part of FormCalc
* last modified 3 Aug 11 th


	double complex function CTCcut(rank, num, numeps,
     &    p1, p2, m1, m2, m3)
	implicit none
	integer rank
	external num, numeps
	integer*8 p1, p2
	double precision m1, m2, m3

	double precision pj(0:3,3), mj(3)
	logical stable, forcemp
	double complex res(0:2), res1
	external CTmpnum

	double precision rootsvalue, muscale
	common /cuttools_para/ rootsvalue, muscale

	pj(0,1) = 0
	pj(1,1) = 0
	pj(2,1) = 0
	pj(3,1) = 0
	call CTMomDecode(p1, pj(0,2))
	call CTMomDecode(p2, pj(0,3))

	mj(1) = m1
	mj(2) = m2
	mj(3) = m3

	PRINT *, "Ccut:", rank
	PRINT *, "p1=", pj(0,2)
	PRINT *, "   ", pj(1,2)
	PRINT *, "   ", pj(2,2)
	PRINT *, "   ", pj(3,2)
	PRINT *, "p2=", pj(0,3)
	PRINT *, "   ", pj(1,3)
	PRINT *, "   ", pj(2,3)
	PRINT *, "   ", pj(3,3)
	PRINT *, "m1=", m1
	PRINT *, "m2=", m2
	PRINT *, "m3=", m3

	call ctsxcut(rootsvalue, muscale, 3, num, CTmpnum,
     &    rank, pj, mj, res, res1, stable, forcemp)
	CTCcut = res(0)
	PRINT *, "Ccut=", CTCcut
	end

