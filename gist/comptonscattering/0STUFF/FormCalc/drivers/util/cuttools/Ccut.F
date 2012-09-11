* Ccut.F
* the FormCalc wrapper for the CutTools three-point function
* this file is part of FormCalc
* last modified 18 Nov 09 th


	double complex function Ccut(rank, num,
     &    p1, p2, m1, m2, m3)
	implicit none
	integer rank
	external num
	integer*8 p1, p2
	double precision m1, m2, m3

	double precision pi(0:3,0:2), mi(3)
	external mpnum

	pi(0,0) = 0
	pi(1,0) = 0
	pi(2,0) = 0
	pi(3,0) = 0
	call MomDecode(p1, pi(0,1))
	call MomDecode(p2, pi(0,2))

	mi(1) = m1
	mi(2) = m2
	mi(3) = m3

	PRINT *, "Ccut:", rank
	PRINT *, "p1=", pi(0,1)
	PRINT *, "   ", pi(1,1)
	PRINT *, "   ", pi(2,1)
	PRINT *, "   ", pi(3,1)
	PRINT *, "p2=", pi(0,2)
	PRINT *, "   ", pi(1,2)
	PRINT *, "   ", pi(2,2)
	PRINT *, "   ", pi(3,2)
	PRINT *, "m1=", m1
	PRINT *, "m2=", m2
	PRINT *, "m3=", m3

c	call ctsxcut(3, num, mpnum, rank, pi, mi, Ccut)
	call ctsxcut(3, num, mpnum, 3, pi, mi, Ccut)
	PRINT *, "Ccut=", Ccut
	end

