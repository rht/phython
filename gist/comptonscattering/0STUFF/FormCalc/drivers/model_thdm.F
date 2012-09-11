* model_thdm.F
* all parameters for the THDM are initialized here
* this file is part of FormCalc
* Abdesslam Arhrib, 2 Jun 03
* last modified 1 Aug 08 th


* The input parameters for model_thdm are:

* SA = sin alpha
* TB = tan beta
* Mh0 = Mass(H1)
* MHH = Mass(H2)
* MA0 = Mass(H3)
* MHp = Mass(Hp)
* Lambda5

* The following preprocessor variables can be defined
* (they are usually set in run.F):

* MODEL_TYPE_I, MODEL_TYPE_II
*   Which type of two-higgs-doublet model to select.

* NO_EXCLUSION_LIMITS
*   If defined, the experimental exclusion limits are ignored.


* Exclusion limits (only checked if NO_EXCLUSION_LIMITS is not set):

* PDG
#define DeltaRho_UPPERBOUND .003D0


#define ModelConstIni SMConstIni
#define ModelDigest SMDigest
#include "model_sm.F"
#undef ModelConstIni
#undef ModelDigest


************************************************************************

	subroutine ModelConstIni(fail)
	implicit none
	integer fail

#include "decl.h"

	double precision Lambda1, Lambda2, Lambda3, Lambda4, Lambda6
	double precision DeltaRho

#if !defined(MODEL_TYPE_I) && !defined(MODEL_TYPE_II)
	stop "Neither MODEL_TYPE_I nor MODEL_TYPE_II defined"
#elif defined(MODEL_TYPE_I) && defined(MODEL_TYPE_II)
	stop "Both MODEL_TYPE_I and MODEL_TYPE_II defined"
#endif

	call SMConstIni(fail)
	if( fail .ne. 0 ) return

	TB2 = TB**2
	CB2 = 1/(1 + TB2)
	SB2 = CB2*TB2
	CB = sqrt(CB2)
	SB = CB*TB
	S2B = 2*TB*CB2
	C2B = (CB - SB)*(CB + SB)

	SA2 = SA**2
	CA2 = (1 - SA)*(1 + SA)
	CA = sqrt(CA2)
	S2A = 2*SA*CA
	C2A = (CA - SA)*(CA + SA)
	SAB = SA*CB + CA*SB
	CAB = CA*CB - SA*SB
	SBA = SB*CA - CB*SA
	CBA = CB*CA + SB*SA

#ifdef MODEL_TYPE_I
	Yuk1 = CA/SB
	Yuk2 = SA/SB
	Yuk3 = -1/TB
#else
	Yuk1 = -SA/CB
	Yuk2 = CA/CB
	Yuk3 = TB
#endif

	Mh02 = Mh0**2
	MH = Mh0
	MH2 = Mh02
	MHH2 = MHH**2
	MA02 = MA0**2
	MHp2 = MHp**2

#ifndef NO_EXCLUSION_LIMITS
	Lambda1 = .25D0*Lambda5*(1 - TB2) +
     &    Alfa*Pi/(2*MW2*SW2)/CB2*
     &      (MHH2*CA2 + Mh02*SA2 - .5D0*(MHH2 - Mh02)*S2A/TB)
	Lambda2 = .25D0*Lambda5*(1 - 1/TB2) +
     &    Alfa*Pi/(2*MW2*SW2)/SB2*
     &      (Mh02*CA2 + MHH2*SA2 - .5D0*(MHH2 - Mh02)*S2A*TB)
        Lambda3 = -.25D0*Lambda5 +
     &    Alfa*Pi/(2*MW2*SW2)*(MHH2 - Mh02)*S2A/S2B
        Lambda4 = (2*Alfa*Pi)/(MW2*SW2)*MHp2
        Lambda6 = (2*Alfa*Pi)/(MW2*SW2)*MA02

	if( max(abs(Lambda1), abs(Lambda2), abs(Lambda3),
     &          Lambda4, abs(Lambda5), Lambda6) .gt. 8*pi ) then
	  Error("Perturbativity violated")
	  return
	endif

	if( abs(MA02 - MHH2) .gt. 1D-13 .and.
     &      abs(MHp2 - MHH2) .gt. 1D-13 ) then
	  DeltaRho = GF/(8*sqrt2*pi**2)*(
     &      MHp2*(1 - MA02/(MHp2 - MA02)*log(MHp2/MA02)) +
     &      CBA**2*Mh02*(MA02/(MA02 - Mh02)*log(MA02/Mh02) -
     &                   MHp2/(MHp2 - Mh02)*log(MHp2/Mh02)) +
     &      SBA**2*MHH2*(MA02/(MA02 - MHH2)*log(MA02/MHH2) -
     &                   MHp2/(MHp2 - MHH2)*log(MHp2/MHH2)) )

	  if( DeltaRho .gt. DeltaRho_UPPERBOUND ) then
	    Error("DeltaRho > DeltaRho_UPPERBOUND")
	    return
	  endif
	endif
#endif
	end


************************************************************************

	subroutine ModelDigest
	implicit none

#include "model.h"

	call SMDigest

	INFO "THDM parameters"
	INFO "----------------------------"
#ifdef MODEL_TYPE_I
	INFO "using type-I model"
#else
	INFO "using type-II model"
#endif
	INFO "TB =", TB
	INFO "SA =", SA
	INFO "Mh0 =", Mh0
	INFO "MHH =", MHH
	INFO "MA0 =", MA0
	INFO "MHp =", MHp
	INFO "Lambda5 =", Lambda5
	INFO "============================"
	end

