* lumi_photon.F
* parton luminosity for a photon collider
* this file is part of FormCalc
* last modified 13 Jul 07 th


#define PARTON_T -1
#define PARTON_X -2
#define MINVAR PARTON_X


************************************************************************

	subroutine LumiDefaults
	implicit none

#include "xsection.h"

	if( Lower(SQRTS) .ne. Upper(SQRTS) )
     &    stop "Cannot handle a range of beam energies"
	end


************************************************************************

	subroutine LumiIni(fail)
	implicit none
	integer fail

#include "xsection.h"

	call SetEnergy(fail, Var(SQRTS))
	if( fail .ne. 0 ) return

	Lower(PARTON_T) = threshold/sqrtS
	Upper(PARTON_T) = 1
	Step(PARTON_T) = 0

	Lower(PARTON_X) = 0
	Upper(PARTON_X) = 1
	Step(PARTON_X) = 0
	end


************************************************************************

	subroutine Luminosity(fac)
	implicit none
	double precision fac

#include "xsection.h"

	integer fail
	double precision tau, x, xrange, beamenergy

	double precision CompAZ
	external CompAZ

	call SetEnergy(fail, Var(PARTON_T)*Var(SQRTS))
	if( fail .ne. 0 ) then
	  fac = 0
	  return
	endif

	tau = Var(PARTON_T)**2
	xrange = 1 - tau
	x = tau + xrange*Var(PARTON_X)

	beamenergy = .5D0*Var(SQRTS)

	fac = 2*fac*Var(PARTON_T)*xrange/x*
     &    CompAZ(x, beamenergy, 0)*
     &    CompAZ(tau/x, beamenergy, 0)
	end

