* ThreeMom.F
* this file is part of FormCalc
* last modified 30 Nov 11 th

#include "types.h"


* ThreeMom computes the length of \vec p_b in the frame in which
* \vec p_a + \vec p_b vanishes.  With sqrtS = sqrt((p_a + p_b)^2),
*
* ThreeMom(sqrtS, ma, mb) = sqrt(lambda(sqrtS, ma**2, mb**2))/(2*sqrtS),
*
* where lambda is the Kallen function.


	RealType function ThreeMom(sqrtS, ma, mb)
	implicit none
	RealType sqrtS, ma, mb

	RealType kb0

c	kb0 = ((sqrtS - ma)*(sqrtS + ma) + m2)/(2*sqrtS)
	kb0 = .5D0*(sqrtS - (ma - mb)*(ma + mb)/sqrtS)
	ThreeMom = sqrt((kb0 - mb)*(kb0 + mb))
	end

