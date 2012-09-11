* TInvariant.F
* calculates an t-type invariant, i.e. t_{ij} = (k_i - k_j)^2
* this file is part of FormCalc
* last modified 20 Oct 05 th


	double precision function TInvariant(i, j)
	implicit none
	integer i, j

	double complex vec(2,2, 8, 0:1)
	common /vectors/ vec

#define Ki(x,y) vec(x,y, 1, i)
#define Kj(x,y) vec(x,y, 1, j)

	double precision T11, T22
	double complex T12

	T11 = DBLE(Ki(1,1)) - DBLE(Kj(1,1))
	T22 = DBLE(Ki(2,2)) - DBLE(Kj(2,2))
	T12 = Ki(1,2) - Kj(1,2)

        TInvariant = T11*T22 - DBLE(T12)**2 - DIMAG(T12)**2
	end

