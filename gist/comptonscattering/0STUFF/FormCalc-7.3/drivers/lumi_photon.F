* lumi_photon.F
* parton luminosity for a photon collider
* this file is part of FormCalc
* last modified 30 Nov 11 th


#define PARTON_T -1
#define PARTON_X -2
#define MINVAR PARTON_X


************************************************************************

	subroutine LumiDefaults(argc, argv)
	implicit none
	integer argc
	character*128 argv(*)

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
	RealType fac

#include "xsection.h"

	integer fail
	RealType tau, x, xrange, beamenergy

	RealType CompAZ
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

