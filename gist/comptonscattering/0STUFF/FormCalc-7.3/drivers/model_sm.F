* model_sm.F
* all parameters for the Standard Model are initialized here
* this file is part of FormCalc
* last modified 30 Nov 11 th


	subroutine ModelDefaults(argc, argv)
	implicit none
	integer argc
	character*128 argv(*)

#include "decl.h"

#ifdef NO_CKM

	CKM(1,1) = 1
	CKM(2,2) = 1
	CKM(3,3) = 1
	CKM(2,1) = 0
	CKM(3,1) = 0
	CKM(1,2) = 0
	CKM(3,2) = 0
	CKM(1,3) = 0
	CKM(2,3) = 0

#else

	RealType s12, s23, s13
	RealType c12, c23, c13
	ComplexType s13c

	s12 = CKMlambda
	c12 = sqrt((1 - s12)*(1 + s12))

	s23 = s12**2*CKMA
	c23 = sqrt((1 - s23)*(1 + s23))

	s13c = s12*s23*ToComplex(CKMrhobar, CKMetabar)
	s13c = c23/c12*s13c/(1 - s23/s12*s13c)
	s13 = abs(s13c)
	c13 = sqrt((1 - s13)*(1 + s13))

	CKM(1,1) = c12*c13
	CKM(2,1) = -s12*c23 - c12*s23*s13c
	CKM(3,1) = s12*s23 - c12*c23*s13c
	CKM(1,2) = s12*c13
	CKM(2,2) = c12*c23 - s12*s23*s13c
	CKM(3,2) = -c12*s23 - s12*c23*s13c
	CKM(1,3) = Conjugate(s13c)
	CKM(2,3) = s23*c13
	CKM(3,3) = c23*c13

#endif

	Mf(1,1) = 0
	Mf(2,1) = ME
	Mf(3,1) = MU
	Mf(4,1) = MD
	Mf(1,2) = 0
	Mf(2,2) = MM
	Mf(3,2) = MC
	Mf(4,2) = MS
	Mf(1,3) = 0
	Mf(2,3) = ML
	Mf(3,3) = MT
	Mf(4,3) = MB

	Mf2(1,1) = 0
	Mf2(2,1) = ME2
	Mf2(3,1) = MU2
	Mf2(4,1) = MD2
	Mf2(1,2) = 0
	Mf2(2,2) = MM2
	Mf2(3,2) = MC2
	Mf2(4,2) = MS2
	Mf2(1,3) = 0
	Mf2(2,3) = ML2
	Mf2(3,3) = MT2
	Mf2(4,3) = MB2

	EL = sqrt(4*pi*Alfa)
	SW = sqrt(SW2)
	end


************************************************************************

	subroutine ModelConstIni(fail)
	implicit none
	integer fail

#include "decl.h"

	fail = 0

	MH2 = MH**2
	end


************************************************************************

	subroutine ModelVarIni(fail, scale)
	implicit none
	integer fail
	RealType scale

#include "decl.h"

	integer nf

	fail = 0

	call AlphaS(Alfas, nf, scale)
	Alfas2 = Alfas**2
	GS = sqrt(4*pi*Alfas)
	end

#include "AlphaS.F"


************************************************************************

	subroutine ModelDigest
	implicit none

#include "decl.h"

	INFO "============================"
	INFO "SM parameters"
	INFO "----------------------------"
	INFO "1/Alfa =", 1/Alfa
	INFO "Alfas =", Alfas
	INFO "GF =", GF
	INFO "SW2 =", SW2
	INFO "----------------------------"
	INFO "Gauge-boson masses:"
	INFO "MW =", MW
	INFO "MZ =", MZ
	INFO "----------------------------"
	INFO "Fermion masses:"
	INFO "ME =", ME
	INFO "MM =", MM
	INFO "ML =", ML
	INFO "MU =", MU
	INFO "MC =", MC
	INFO "MT =", MT
	INFO "MD =", MD
	INFO "MS =", MS
	INFO "MB =", MB
	INFO "----------------------------"
	INFO "Higgs mass:"
	INFO "MH =", MH
	INFO "============================"
	end

