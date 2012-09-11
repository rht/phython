* Acut.F
* the FormCalc wrapper for the CutTools one-point function
* this file is part of FormCalc
* last modified 21 Oct 09 th


	double complex function Acut(rank, num, m)
	implicit none
	integer rank
	external num
	double precision m

	double precision pi(0:3,0:0)
	external mpnum

	pi(0,0) = 0
	pi(1,0) = 0
	pi(2,0) = 0
	pi(3,0) = 0

	PRINT *, "Acut:", rank
	PRINT *, "m1=", m

c	call ctsxcut(1, num, mpnum, rank, pi, m, Acut)
	call ctsxcut(1, num, mpnum, 1, pi, m, Acut)
	PRINT *, "Acut=", Acut
	end

