* Dcut.F
* the FormCalc wrapper for the CutTools four-point function
* this file is part of FormCalc
* last modified 21 Oct 09 th


	double complex function Dcut(rank, num,
     &    p1, p2, p3, m1, m2, m3, m4)
	implicit none
	integer rank
	external num
	integer*8 p1, p2, p3
	double precision m1, m2, m3, m4

	double precision pi(0:3,0:3), mi(4)
	external mpnum

	pi(0,0) = 0
	pi(1,0) = 0
	pi(2,0) = 0
	pi(3,0) = 0
	call MomDecode(p1, pi(0,1))
	call MomDecode(p2, pi(0,2))
	call MomDecode(p3, pi(0,3))

	mi(1) = m1
	mi(2) = m2
	mi(3) = m3
	mi(4) = m4

	PRINT *, "Dcut:", rank
	PRINT *, "p1=", pi(0,1)
	PRINT *, "   ", pi(1,1)
	PRINT *, "   ", pi(2,1)
	PRINT *, "   ", pi(3,1)
	PRINT *, "p2=", pi(0,2)
	PRINT *, "   ", pi(1,2)
	PRINT *, "   ", pi(2,2)
	PRINT *, "   ", pi(3,2)
	PRINT *, "p3=", pi(0,3)
	PRINT *, "   ", pi(1,3)
	PRINT *, "   ", pi(2,3)
	PRINT *, "   ", pi(3,3)
	PRINT *, "m1=", m1
	PRINT *, "m2=", m2
	PRINT *, "m3=", m3
	PRINT *, "m4=", m4

c	call ctsxcut(4, num, mpnum, rank, pi, mi, Dcut)
	call ctsxcut(4, num, mpnum, 4, pi, mi, Dcut)
	PRINT *, "Dcut=", Dcut
	end

