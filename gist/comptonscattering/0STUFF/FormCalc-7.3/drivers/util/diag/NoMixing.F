* NoMixing.f
* same invocation as Eigen2x2, but ignores M12 (i.e., no mixing)
* used for testing
* this file is part of FormCalc
* last modified 30 Nov 11 th

#include "diag-f.h"


	subroutine NoMixing(fail, m, msq, U, M11, M22, M12)
	implicit none
	integer fail
	RealType m(2), msq(2)
	ComplexType U(2,2)
	RealType M11, M22
	ComplexType M12

	if( M11 .gt. M22 ) then
	  msq(1) = M22
	  msq(2) = M11
	  U(1,1) = 0
	  U(2,2) = 0
	  U(1,2) = 1
	  U(2,1) = 1
	else
	  msq(1) = M11
	  msq(2) = M22
	  U(1,1) = 1
	  U(2,2) = 1
	  U(1,2) = 0
	  U(2,1) = 0
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

