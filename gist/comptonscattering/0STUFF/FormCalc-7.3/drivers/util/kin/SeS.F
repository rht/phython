* SeS.F
* computes Spinor.eps.Spinor
* this file is part of FormCalc
* last modified 30 Nov 11 th

#include "types.h"


	ComplexType function SeS(a, b)
	implicit none
	integer a, b

	ComplexType spi(2, 1)
	common /vectors/ spi

#define A(x) spi(x, a)
#define B(x) spi(x, b)

	SeS = A(1)*B(2) - A(2)*B(1)
	end

