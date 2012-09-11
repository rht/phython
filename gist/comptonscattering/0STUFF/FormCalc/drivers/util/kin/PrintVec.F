* PrintVec.F
* print out the Minkowskian representation of a vector
* this file is part of FormCalc
* last modified 19 Oct 09 th

* Note: the arguments a, b specify the vector's position in
* the array vec, NOT the vector itself.  The macros k, e, ec, s
* (defined in util.h) convert from 1...LEGS to a, b.

	subroutine PrintVec(a, s)
	implicit none
	integer a
	character*(*) s

	double complex vec(2,2, 1)
	common /vectors/ vec

#define A(i,j) vec(i,j,a)

	print *, s, " =", .5D0*(DBLE(A(1,1)) + DBLE(A(2,2)))
	print *, "    ", DBLE(A(2,1))
	print *, "    ", DIMAG(A(2,1))
	print *, "    ", .5D0*(DBLE(A(1,1)) - DBLE(A(2,2)))
	end

