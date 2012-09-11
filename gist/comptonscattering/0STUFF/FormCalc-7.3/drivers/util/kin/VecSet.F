* VecSet.F
* construct the CMS representation of a momentum and the
* three corresponding polarization vectors
* this file is part of FormCalc
* last modified 30 Nov 11 th

#include "types.h"


* i: the index of the momentum
* m, p: mass and three-momentum of the particle
* ex,ey,ez: the unit three-vector of the momentum

	subroutine VecSet(i, m, p, ex,ey,ez)
	implicit none
	integer i
	RealType m, p, ex,ey,ez

	ComplexType vec(2,2, 8, 0:1)
	common /vectors/ vec

#define K(x,y)     vec(x,y, 1, i)
#define EpsTR(x,y) vec(x,y, 2, i)
#define EpsL(x,y)  vec(x,y, 3, i)
#define EpsTL(x,y) vec(x,y, 4, i)

	ComplexType spi(2, 16, 0:1)
	equivalence (spi, vec)

#define Spinor(hel,om,x) spi(x, hel*2+om+5, i)
#define DottedSpinor(hel,om,x) spi(x, hel*2+om+7, i)

* i.e., spi(x,  9, i) = left-handed,  undotted, upper
*       spi(x, 10, i) = left-handed,  undotted, lower
*       spi(x, 11, i) = left-handed,  dotted,   upper
*       spi(x, 12, i) = left-handed,  dotted,   lower
*       spi(x, 13, i) = right-handed, undotted, upper
*       spi(x, 14, i) = right-handed, undotted, lower
*       spi(x, 15, i) = right-handed, dotted,   upper
*       spi(x, 16, i) = right-handed, dotted,   lower

	RealType momspec(16,1)
	common /momenta/ momspec

* the following must match the definitions in util.h:
#define SPEC_M 1
#define SPEC_K 2
#define SPEC_E 3
#define SPEC_KT 4
#define SPEC_ET 5
#define SPEC_PRAP 6
#define SPEC_RAP 7
#define SPEC_DELTAK 8
#define SPEC_PHI 9
#define SPEC_EX 10
#define SPEC_EY 11
#define SPEC_EZ 12

	RealType m2, p2, p0, deltap, sump
	RealType sinth, onePez, oneMez
	ComplexType expIphi

	RealType sqrt2
	parameter (sqrt2 = 1.4142135623730950488016887242096981D0)

	momspec(SPEC_M,i) = m
	momspec(SPEC_K,i) = p

	m2 = m**2
	if( m2 .lt. 1D-14 ) then
	  p0 = p
	  deltap = 0
	else
	  p2 = p**2
	  p0 = sqrt(p2 + m2)
	  deltap = m2/(p0 + p)
	endif

	momspec(SPEC_E,i) = p0
	momspec(SPEC_DELTAK,i) = deltap
	momspec(SPEC_EX,i) = ex
	momspec(SPEC_EY,i) = ey
	momspec(SPEC_EZ,i) = ez

	sinth = ex**2 + ey**2
	onePez = 1 + ez
	if( onePez .lt. .5D0 ) onePez = sinth/(1 - ez)
	oneMez = 1 - ez
	if( oneMez .lt. .5D0 ) oneMez = sinth/(1 + ez)

	if( sinth .lt. 1D-14 ) then
* phi is irrelevant when theta = 0
	  expIphi = 1
	  momspec(SPEC_PHI,i) = 0
	else
	  sinth = sqrt(sinth)
	  expIphi = ToComplex(ex, ey)/sinth
	  momspec(SPEC_PHI,i) = atan2(ey, ex)
	endif

	K(1,1) = p0*onePez - deltap*ez
	K(2,2) = p0*oneMez + deltap*ez
	K(2,1) = p*ToComplex(ex, ey)
	K(1,2) = Conjugate(K(2,1))

	if( m .ne. 0 ) then
	  EpsL(1,1) = (p*onePez + deltap*ez)/m
	  EpsL(2,2) = (p*oneMez - deltap*ez)/m
	  EpsL(2,1) = p0/m*ToComplex(ex, ey)
	  EpsL(1,2) = Conjugate(EpsL(2,1))
	endif

	EpsTR(2,2) = 1/sqrt2*sinth
	EpsTL(2,2) = EpsTR(2,2)
	EpsTR(1,1) = -EpsTR(2,2)
	EpsTL(1,1) = EpsTR(1,1)
	EpsTR(2,1) = -1/sqrt2*oneMez*expIphi
	EpsTL(1,2) = Conjugate(EpsTR(2,1))
	EpsTL(2,1) = 1/sqrt2*onePez*expIphi
	EpsTR(1,2) = Conjugate(EpsTL(2,1))

#define safelog(x) log(max(x, 1D-300))

	momspec(SPEC_KT,i) = p*sinth
	momspec(SPEC_ET,i) = sqrt(momspec(SPEC_KT,i)**2 + m2)
	momspec(SPEC_PRAP,i) = safelog(onePez) -
     &    safelog(sinth)
	momspec(SPEC_RAP,i) = safelog(Re(K(1,1))) -
     &    safelog(momspec(SPEC_ET,i))

* this is E^(I phi) cos(th/2) = 1/sqrt2*sqrt(1 + ez)*expIphi:
	expIphi = 1/sqrt2*sqrt(onePez)*expIphi
* this is sin(th/2):
	sinth = 1/sqrt2*sqrt(oneMez)
	sump = sqrt(p0 + p)
	deltap = sqrt(deltap)

	Spinor(-1, 6, 1) = deltap*sinth
	DottedSpinor(-1, 6, 1) = Spinor(-1, 6, 1)
	Spinor(-1, 6, 2) = -deltap*expIphi
	DottedSpinor(-1, 6, 2) = Conjugate(Spinor(-1, 6, 2))

	Spinor(-1, 7, 1) = sump*sinth
	DottedSpinor(-1, 7, 1) = Spinor(-1, 7, 1)
	Spinor(-1, 7, 2) = -sump*expIphi
	DottedSpinor(-1, 7, 2) = Conjugate(Spinor(-1, 7, 2))

	DottedSpinor(1, 6, 1) = sump*expIphi
	Spinor(1, 6, 1) = Conjugate(DottedSpinor(1, 6, 1))
	DottedSpinor(1, 6, 2) = sump*sinth
	Spinor(1, 6, 2) = DottedSpinor(1, 6, 2)

	DottedSpinor(1, 7, 1) = deltap*expIphi
	Spinor(1, 7, 1) = Conjugate(DottedSpinor(1, 7, 1))
	DottedSpinor(1, 7, 2) = deltap*sinth
	Spinor(1, 7, 2) = DottedSpinor(1, 7, 2)
	end

