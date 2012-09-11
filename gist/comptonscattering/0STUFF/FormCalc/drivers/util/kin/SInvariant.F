* SInvariant.F
* calculates an s-type invariant, i.e. s_{ij} = (k_i + k_j)^2
* this file is part of FormCalc
* last modified 20 Oct 05 th


	double precision function SInvariant(i, j)
	implicit none
	integer i, j

	double complex vec(2,2, 8, 0:1)
	common /vectors/ vec

#define Ki(x,y) vec(x,y, 1, i)
#define Kj(x,y) vec(x,y, 1, j)

	double precision S11, S22
	double complex S12

	S11 = DBLE(Ki(1,1)) + DBLE(Kj(1,1))
	S22 = DBLE(Ki(2,2)) + DBLE(Kj(2,2))
	S12 = Ki(1,2) + Kj(1,2)

	SInvariant = S11*S22 - DBLE(S12)**2 - DIMAG(S12)**2
	end

