* Ecut.F
* the FormCalc wrapper for the CutTools five-point function
* this file is part of FormCalc
* last modified 21 Oct 09 th


	double complex function Ecut(rank, num,
     &    p1, p2, p3, p4, m1, m2, m3, m4, m5)
	implicit none
	integer rank
	external num
	integer*8 p1, p2, p3, p4
	double precision m1, m2, m3, m4, m5

	double precision pi(0:3,0:4), mi(5)
	external mpnum

	pi(0,0) = 0
	pi(1,0) = 0
	pi(2,0) = 0
	pi(3,0) = 0
	call MomDecode(p1, pi(0,1))
	call MomDecode(p2, pi(0,2))
	call MomDecode(p3, pi(0,3))
	call MomDecode(p4, pi(0,4))

	mi(1) = m1
	mi(2) = m2
	mi(3) = m3
	mi(4) = m4
	mi(5) = m5

c	call ctsxcut(5, num, mpnum, rank, pi, mi, Ecut)
	call ctsxcut(5, num, mpnum, 5, pi, mi, Ecut)
	end

