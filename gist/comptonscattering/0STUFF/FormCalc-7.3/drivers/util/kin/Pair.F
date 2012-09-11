* Pair.F
* calculates the scalar product of two four-vectors
* this file is part of FormCalc
* last modified 30 Nov 11 th

#include "types.h"


* Note: the arguments a, b specify the vector's position in
* the array vec, NOT the vector itself.  The macros k, e, ec, s
* (defined in util.h) convert from 1...LEGS to a, b.

	ComplexType function Pair(a, b)
	implicit none
	integer a, b

	ComplexType vec(2,2, 1)
	common /vectors/ vec

#define A(x,y) vec(x,y, a)
#define B(x,y) vec(x,y, b)

	Pair = .5D0*(
     &    A(1,1)*B(2,2) + A(2,2)*B(1,1) -
     &    A(1,2)*B(2,1) - A(2,1)*B(1,2) )
	end

