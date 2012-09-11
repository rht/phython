* Pair.F
* calculates the scalar product of two four-vectors
* this file is part of FormCalc
* last modified 20 Oct 05 th


* Note: the arguments a, b specify the vector's position in
* the array vec, NOT the vector itself.  The macros k, e, ec, s
* (defined in util.h) convert from 1...LEGS to a, b.

	double complex function Pair(a, b)
	implicit none
	integer a, b

	double complex vec(2,2, 1)
	common /vectors/ vec

#define A(x,y) vec(x,y, a)
#define B(x,y) vec(x,y, b)

	Pair = .5D0*(
     &    DBLE(A(1,1))*DBLE(B(2,2)) +
     &    DBLE(A(2,2))*DBLE(B(1,1)) -
     &    A(1,2)*B(2,1) - A(2,1)*B(1,2) )
	end

