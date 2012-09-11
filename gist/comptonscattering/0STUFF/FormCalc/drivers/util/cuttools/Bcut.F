* Bcut.F
* the FormCalc wrapper for the CutTools two-point function
* this file is part of FormCalc
* last modified 21 Oct 09 th


	double complex function Bcut(rank, num, p, m1, m2)
	implicit none
	integer rank
	external num
	integer*8 p
	double precision m1, m2

	double precision pi(0:3,0:1), mi(2)
	external mpnum

	pi(0,0) = 0
	pi(1,0) = 0
	pi(2,0) = 0
	pi(3,0) = 0
	call MomDecode(p, pi(0,1))

	mi(1) = m1
	mi(2) = m2

	PRINT *, "Bcut:", rank
	PRINT *, "p1=", pi(0,1)
	PRINT *, "   ", pi(1,1)
	PRINT *, "   ", pi(2,1)
	PRINT *, "   ", pi(3,1)
	PRINT *, "m1=", m1
	PRINT *, "m2=", m2

c	call ctsxcut(2, num, mpnum, rank, pi, mi, Bcut)
	call ctsxcut(2, num, mpnum, 2, pi, mi, Bcut)
	PRINT *, "Bcut=", Bcut
	end

