* lumi_hadron.F
* parton luminosity for a hadron collider
* this file is part of FormCalc
* last modified 13 Feb 07 th

* This code links to the LHAPDF library.
* The following preprocessor variables are used:
*
* - PARTON1, PARTON2 are the PDG parton codes,
* - SCALE determines the factorization scale,
* - PDFSET specifies the file name of the PDF-set to be used.


#define PARTON_T -1
#define PARTON_X -2
#define MINVAR PARTON_X

#define HAVE_ALPHASPDF


************************************************************************

	subroutine LumiDefaults
	implicit none

#include "xsection.h"

	integer partonsign
	common /hadron/ partonsign

	if( Lower(SQRTS) .ne. Upper(SQRTS) )
     &    stop "Cannot handle a range of beam energies"

	if( Lower(SQRTS) .eq. 2000 ) then
	  partonsign = -1
	else if( Lower(SQRTS) .eq. 14000 ) then
	  partonsign = 1
	else
	  stop "Unknown collider"
	endif

	call InitPDFsetByName(PDFSET)
	call InitPDF(PDFMEM)
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

	double precision xPDF(-6:6), tPDF(-6:6)
	double precision tau, x, xrange
	integer fail

	integer partonsign
	common /hadron/ partonsign

	call SetEnergy(fail, Var(PARTON_T)*Var(SQRTS))
	if( fail .ne. 0 ) then
	  fac = 0
	  return
	endif

	tau = Var(PARTON_T)**2
	xrange = 1 - tau
	x = tau + xrange*Var(PARTON_X)

	if( x .lt. 1D-10 .or. x .gt. (1 - 1D-10) ) then
	  fac = 0
	  return
	endif

	call evolvePDF(x, scale, xPDF)
	call evolvePDF(tau/x, scale, tPDF)

	fac = fac*xrange/(Var(PARTON_T)*x)*
     &    ( tPDF(partonsign*PARTON1)*xPDF(PARTON2) +
     &      xPDF(partonsign*PARTON2)*tPDF(PARTON1) )
#if PARTON1 != PARTON2
     &    * 2
#endif
	end

