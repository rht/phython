* Fcut.F
* the FormCalc wrapper for the CutTools six-point function
* this file is part of FormCalc
* last modified 21 Oct 09 th


	double complex function Fcut(rank, num,
     &    p1, p2, p3, p4, p5, m1, m2, m3, m4, m5, m6)
	implicit none
	integer rank
	external num
	integer*8 p1, p2, p3, p4, p5
	double precision m1, m2, m3, m4, m5, m6

	double precision pi(0:3,0:5), mi(6)
	external mpnum

	pi(0,0) = 0
	pi(1,0) = 0
	pi(2,0) = 0
	pi(3,0) = 0
	call MomDecode(p1, pi(0,1))
	call MomDecode(p2, pi(0,2))
	call MomDecode(p3, pi(0,3))
	call MomDecode(p4, pi(0,4))
	call MomDecode(p5, pi(0,5))

	mi(1) = m1
	mi(2) = m2
	mi(3) = m3
	mi(4) = m4
	mi(5) = m5
	mi(6) = m6

c	call ctsxcut(6, num, mpnum, rank, pi, mi, Fcut)
	call ctsxcut(6, num, mpnum, 6, pi, mi, Fcut)
	end

