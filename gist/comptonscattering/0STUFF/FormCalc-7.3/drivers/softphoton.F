* softphoton.F
* routines for the treatment of soft-photon radiation
* this file is part of FormCalc
* last modified 30 Nov 11 th


************************************************************************
** SoftPhotonFactor adds up the contributions from soft-photon emission
** off the different legs.  The result then has to be multiplied with
** the Born cross-section.

	RealType function SoftPhotonFactor()
	implicit none

#include "xsection.h"

	RealType PhaseSpaceIntegral
	external PhaseSpaceIntegral

	integer i, j

	RealType loglambda
	common /infrared/ loglambda

	RealType charge(LEGS), r

* not a very nice solution, the following, but necessary
* because there may be no "/" in a constant in the data
* statement (such as in 2/3D0)
#if LEGS == 3
	RealType ch1, ch2, ch3
	parameter (ch1 = CHARGE1)
	parameter (ch2 = CHARGE2)
	parameter (ch3 = CHARGE3)
	data charge / ch1, ch2, ch3 /
#elif LEGS == 4
	RealType ch1, ch2, ch3, ch4
	parameter (ch1 = CHARGE1)
	parameter (ch2 = CHARGE2)
	parameter (ch3 = CHARGE3)
	parameter (ch4 = CHARGE4)
	data charge / ch1, ch2, ch3, ch4 /
#elif LEGS == 5
	RealType ch1, ch2, ch3, ch4, ch5
	parameter (ch1 = CHARGE1)
	parameter (ch2 = CHARGE2)
	parameter (ch3 = CHARGE3)
	parameter (ch4 = CHARGE4)
	parameter (ch5 = CHARGE5)
	data charge / ch1, ch2, ch3, ch4, ch5 /
#else
#error No definition for LEGS external legs
#endif

	loglambda = log(4*(ESOFTMAX)**2/getlambda())

	SoftPhotonFactor = 0
	do i = 1, LEGS
	  do j = i, LEGS
	    r = charge(i)*charge(j)
	    if( r .ne. 0 ) SoftPhotonFactor = SoftPhotonFactor +
     &        sign(1, LEGS_IN - i)*sign(1, LEGS_IN - j)*
     &        r*PhaseSpaceIntegral(i, j)
	  enddo
	enddo
	SoftPhotonFactor = -4*pi*Alfa/(2*pi)**3*SoftPhotonFactor
	end


************************************************************************
** PhaseSpaceIntegral calculates the integral over the soft photon's
** phase space.  See A. Denner, Fortschr. d. Physik, 41 (1993) 4.

	RealType function PhaseSpaceIntegral(i, j)
	implicit none
	integer i, j

#include "xsection.h"

	RealType mi, ki0, ki, kidelta, kilog
	RealType mj, kj0, kj, kjdelta, kjlog
	RealType kikj, alpha, d, vi

	RealType loglambda
	common /infrared/ loglambda

	mi = momspec(SPEC_M, i)
	ki0 = momspec(SPEC_E, i)
	ki = momspec(SPEC_K, i)
	kidelta = momspec(SPEC_DELTAK, i)
	kilog = log(kidelta/(ki0 + ki))

	if( i .eq. j ) then
	  if( ki .lt. 1D-14 ) then
	    PhaseSpaceIntegral = pi*(loglambda - 2)
	  else
	    PhaseSpaceIntegral = pi*(loglambda + ki0/ki*kilog)
	  endif
	  return
	endif

	mj = momspec(SPEC_M, j)
	kj0 = momspec(SPEC_E, j)
	kj = momspec(SPEC_K, j)
	kjdelta = momspec(SPEC_DELTAK, j)
	kjlog = log(kjdelta/(kj0 + kj))

	kikj = Re(Pair(k(i), k(j)))

	d = sqrt((kikj - mi*mj)*(kikj + mi*mj))
	alpha = (kikj + d)/mi**2
	if( alpha*ki0 .lt. kj0 ) alpha = (kikj - d)/mi**2

	d = 1/((alpha*mi - mj)*(alpha*mi + mj))
	vi = 2*(alpha*ki0 - kj0)*d

	PhaseSpaceIntegral = 4*pi*alpha*kikj*d * Re(
     &    log(mi/mj*abs(alpha))*loglambda +
     &    .25D0*(kilog - kjlog)*(kilog + kjlog) +
     &    Li2omx((ki0 + ki)*alpha*vi) + Li2omx(kidelta*alpha*vi) -
     &    Li2omx((kj0 + kj)*vi) - Li2omx(kjdelta*vi) )
	end

