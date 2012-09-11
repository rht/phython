* CTDcut.F
* the FormCalc wrapper for the CutTools four-point function
* this file is part of FormCalc
* last modified 3 Aug 11 th


	double complex function CTDcut(rank, num, numeps,
     &    p1, p2, p3, m1, m2, m3, m4)
	implicit none
	integer rank
	external num, numeps
	integer*8 p1, p2, p3
	double precision m1, m2, m3, m4

	double precision pj(0:3,4), mj(4)
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
	call CTMomDecode(p3, pj(0,4))

	mj(1) = m1
	mj(2) = m2
	mj(3) = m3
	mj(4) = m4

	PRINT *, "CTDcut:", rank
	PRINT *, "p1=", pj(0,2)
	PRINT *, "   ", pj(1,2)
	PRINT *, "   ", pj(2,2)
	PRINT *, "   ", pj(3,2)
	PRINT *, "p2=", pj(0,3)
	PRINT *, "   ", pj(1,3)
	PRINT *, "   ", pj(2,3)
	PRINT *, "   ", pj(3,3)
	PRINT *, "p3=", pj(0,4)
	PRINT *, "   ", pj(1,4)
	PRINT *, "   ", pj(2,4)
	PRINT *, "   ", pj(3,4)
	PRINT *, "m1=", m1
	PRINT *, "m2=", m2
	PRINT *, "m3=", m3
	PRINT *, "m4=", m4

	call ctsxcut(rootsvalue, muscale, 4, num, CTmpnum,
     &    rank, pj, mj, res, res1, stable, forcemp)
	CTDcut = res(0)
	PRINT *, "CTDcut=", CTDcut
	end

