* Eigen2x2.F
* diagonalization of a Hermitian 2-by-2 matrix
* this file is part of FormCalc
* last modified 30 Nov 11 th

#include "diag-f.h"


	subroutine Eigen2x2(fail, m, msq, U, M11, M22, M12)
	implicit none
	integer fail
	RealType m(2), msq(2)
	ComplexType U(2,2)
	RealType M11, M22
	ComplexType M12

	RealType M12sq, delta, t, c, m1, m2
	ComplexType s

	M12sq = DBLE(M12)**2 + DIMAG(M12)**2
	t = .5D0*(M11 - M22)
	t = 1/(t + sign(sqrt(t**2 + M12sq), t))

	delta = t*M12sq
	m1 = M11 + delta
	m2 = M22 - delta

	c = 1/sqrt(delta*t + 1)
	s = c*t*M12

	if( m1 .gt. m2 ) then
	  msq(1) = m2
	  msq(2) = m1
	  U(1,1) = -DCONJG(s)
	  U(2,2) = s
	  U(1,2) = c
	  U(2,1) = c
	else
	  msq(1) = m1
	  msq(2) = m2
	  U(1,1) = c
	  U(2,2) = c
	  U(1,2) = s
	  U(2,1) = -DCONJG(s)
	endif

	if( msq(1) .lt. 0 ) then
	  print *, "Negative mass-squares."
	  fail = 1
	  return
	endif

	m(1) = sqrt(msq(1))
	m(2) = sqrt(msq(2))
	fail = 0
	end

