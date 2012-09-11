* MomDecode.F
* form a linear combination of momenta
* and return it as a Minkowski four-vector
* the momenta's (integer*1) prefactors are encoded as the
* bytes of an integer*8
* this file is part of FormCalc
* last modified 20 Oct 05 th


	subroutine MomDecode(code, p)
	implicit none
	integer*8 code
	double precision p(0:3)

	integer i
	double precision a(4)

	integer*8 c
	integer*1 f
	equivalence (f, c)

	double complex vec(2,2, 1)
	common /vectors/ vec

	a(1) = 0
	a(2) = 0
	a(3) = 0
	a(4) = 0

	c = code
	i = 1

	do while( c .ne. 0 )
	  i = i + 8
	  a(1) = a(1) + f*DBLE(vec(1,1,i))
	  a(2) = a(2) + f*DBLE(vec(2,1,i))
	  a(3) = a(3) + f*DIMAG(vec(2,1,i))
	  a(4) = a(4) + f*DBLE(vec(2,2,i))
	  c = ishft(c, -8)
	enddo

	p(0) = .5D0*(a(1) + a(4))
	p(3) = .5D0*(a(1) - a(4))
	p(1) = a(2)
	p(2) = a(3)
	end

