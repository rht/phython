* SInvariant.F
* calculates an s-type invariant, i.e. s_{ij} = (k_i + k_j)^2
* this file is part of FormCalc
* last modified 30 Nov 11 th

#include "types.h"


	RealType function SInvariant(i, j)
	implicit none
	integer i, j

	ComplexType vec(2,2, 8, 0:1)
	common /vectors/ vec

#define Ki(x,y) vec(x,y, 1, i)
#define Kj(x,y) vec(x,y, 1, j)

	RealType S11, S22
	ComplexType S12

	S11 = Re(Ki(1,1)) + Re(Kj(1,1))
	S22 = Re(Ki(2,2)) + Re(Kj(2,2))
	S12 = Ki(1,2) + Kj(1,2)

	SInvariant = S11*S22 - Sq(S12)
	end

