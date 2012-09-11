* TInvariant.F
* calculates an t-type invariant, i.e. t_{ij} = (k_i - k_j)^2
* this file is part of FormCalc
* last modified 30 Nov 11 th

#include "types.h"


	RealType function TInvariant(i, j)
	implicit none
	integer i, j

	ComplexType vec(2,2, 8, 0:1)
	common /vectors/ vec

#define Ki(x,y) vec(x,y, 1, i)
#define Kj(x,y) vec(x,y, 1, j)

	RealType T11, T22
	ComplexType T12

	T11 = Re(Ki(1,1)) - Re(Kj(1,1))
	T22 = Re(Ki(2,2)) - Re(Kj(2,2))
	T12 = Ki(1,2) - Kj(1,2)

        TInvariant = T11*T22 - Sq(T12)
	end

