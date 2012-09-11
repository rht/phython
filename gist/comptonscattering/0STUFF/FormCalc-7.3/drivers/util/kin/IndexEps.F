* IndexEps.F
* this file is part of FormCalc
* last modified 26 Jan 07 th


	integer function IndexEps(i, j, k)
	implicit none
	integer i, j, k

#define signbit(i) ibits(i,31,1)

	integer cmp, x
	cmp(x) = 2*signbit(x) - signbit(ior(x, -x))

	IndexEps = cmp(j - i)*cmp(j - k)*cmp(k - i)
	end

