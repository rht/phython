* SAMomDecode.F
* form a linear combination of momenta in Samurai conventions
* and return it as a Minkowski four-vector
* the momenta's (integer*1) prefactors are encoded as the
* bytes of an integer*8
* this file is part of FormCalc
* last modified 7 Dec 10 th


	subroutine SAMomDecode(code, p,ldp)
	implicit none
	integer*8 code
	integer ldp
	double precision p(ldp,4)

	integer i
	double complex a(2,2)

	integer*8 c
	integer*1 f
	equivalence (f, c)

	double complex vec(2,2, 1)
	common /vectors/ vec

	a(1,1) = 0
	a(2,1) = 0
	a(1,2) = 0
	a(2,2) = 0

	c = code
	i = 1

	do while( c .ne. 0 )
	  i = i + 8
	  a(1,1) = a(1,1) + f*vec(1,1,i)
	  a(2,1) = a(2,1) + f*vec(2,1,i)
	  a(1,2) = a(1,2) + f*vec(1,2,i)
	  a(2,2) = a(2,2) + f*vec(2,2,i)
	  c = ishft(c, -8)
	enddo

	p(1,4) = .5D0*DBLE(a(1,1) + a(2,2))
	p(1,3) = .5D0*DBLE(a(1,1) - a(2,2))
	p(1,1) = .5D0*DBLE(a(2,1) + a(1,2))
	p(1,2) = .5D0*DIMAG(a(2,1) - a(1,2))
	end

