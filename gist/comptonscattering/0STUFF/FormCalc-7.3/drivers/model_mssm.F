* model_mssm.F
* all parameters for the MSSM are initialized here
* this file is part of FormCalc
* Christian Schappacher, 17 Feb 99
* last modified 30 Nov 11 th


* The input parameters for model_mssm are:

* Higgs sector:
*   TB
*   MA0 or MHp

* Chargino/Neutralino sector:
*   Mino2
*   MUE

* Sfermion sector:
*   MSusy
*   Atau, At, Ab
* If one wants to specify the breaking and mixing in the
* sfermion sector more precisely, the following preprocessor
* variables can be used:
*   MSQ(g), MSU(g), MSD(g), MSL(g), MSE(g) (= MSusy by default)
*   Ae(g) (= Atau by default)
*   Au(g) (= At by default)
*   Ad(g) (= Ab by default)


* The following preprocessor variables can be defined
* (they are usually set in run.F):

* NO_SQUARK_MIXING
*   If defined, sets Atau, At, and Ab such that the off-diagonal entries
*   of the sfermion mixing matrices vanish, i.e. gauge eigenstates =
*   mass eigenstates.

* HIGGS_MASSES can take the values
* TREE -- uses only tree-level Higgs masses,
* SIMPLE -- uses a simple one-loop approximation of the
*    Higgs masses which is valid for all parameters.
* FEYNHIGGS -- invokes FeynHiggs to compute the Higgs masses,
* (undefined) -- uses a fairly precise two-loop approximation,
*    which is valid though only for real parameters.

* SM_HIGGS
*   If defined, reverse-engineers the mixing in the Higgs sector
*   (alpha and beta) such that the MSSM Higgs sector looks like
*   a SM Higgs sector.

* NMFV
*   If defined, initialize also the fully 6x6 mixed sfermion mass
*   matrices with the flavour-violating parameters deltaSf.

* NO_EXCLUSION_LIMITS
*   If defined, the experimental exclusion limits are ignored.

* NO_GUT_RELATION
*   If defined, does not use the GUT relation which relates
*   the gaugino mass paramters.  In this case, Mino1 and Mino3 are
*   also input parameters.


#ifndef Ae
#define Ae(g) Atau
#endif

#ifndef Au
#define Au(g) At
#endif

#ifndef Ad
#define Ad(g) Ab
#endif

#ifndef MSQ
#define MSQ(g) MSusy
#endif

#ifndef MSU
#define MSU(g) MSusy
#endif

#ifndef MSL
#define MSL(g) MSusy
#endif

#ifndef MSD
#define MSD(g) MSusy
#endif

#ifndef MSE
#define MSE(g) MSusy
#endif


* Exclusion limits (only checked if NO_EXCLUSION_LIMITS is not set):

* hep-ex/9910020 (L3)
#define MStop_LOWERBOUND 80

* hep-ex/9910020 (L3)
#define MSbottom_LOWERBOUND 70

* hep-ex/9910020 (L3)
#define MSquark_LOWERBOUND 150

* hep-ex/9910006 (L3)
#define MSlepton_LOWERBOUND 70

* PDG
#define DeltaRho_UPPERBOUND .003D0

#ifdef COMPLEX_PARAMETERS
* Kane & Wang, Phys Lett B488 (2000) 383
#define Mh0_LOWERBOUND 85
#else
* http://lephiggs.web.cern.ch/LEPHIGGS/papers/osaka_note.ps
#define Mh0_LOWERBOUND 91
#endif

* hep-ex/9910007 (L3) and hep-ex/9909051 (OPAL)
#define MNeu_LOWERBOUND 30

* hep-ex/9910007 (L3) and hep-ex/9909051 (OPAL)
#define MCha_LOWERBOUND 90

* hep-ex/9910020 (L3) and PDG
#define MGl_LOWERBOUND 175


#define TREE 100
#define SIMPLE 200
#define FEYNHIGGS 300

#define ModelDefaults SMDefaults
#define ModelConstIni SMConstIni
#define ModelVarIni SMVarIni
#define ModelDigest SMDigest
#include "model_sm.F"
#undef ModelDefaults
#undef ModelConstIni
#undef ModelVarIni
#undef ModelDigest


************************************************************************

	subroutine ModelDefaults(argc, argv)
	implicit none
	integer argc
	character*128 argv(*)

#include "decl.h"

	integer i

	do i = 1, 3*3*3
	  Af_flat(i) = 0
	  Xf_flat(i) = 0
	enddo

#ifdef NMFV
	do i = 1, 6*6*2
	  deltaSf_flat(i) = 0
	enddo
#endif

	TB = -1
	MA0 = -1
	MHp = -1

	call SMDefaults(argc, argv)
	end


************************************************************************

	subroutine ModelConstIni(fail)
	implicit none
	integer fail

#include "decl.h"

	integer nf

	call SMConstIni(fail)
	if( fail .ne. 0 ) return

	call AlphaS(AlfasMT, nf, MT)

	TB2 = TB**2
	CB2 = 1/(1 + TB2)
	SB2 = CB2*TB2
	CB = sqrt(CB2)
	SB = CB*TB
	S2B = 2*TB*CB2
	C2B = (CB - SB)*(CB + SB)

	call SfermionMasses(fail)
	if( fail .ne. 0 ) return

	call CharginoMasses(fail)
	if( fail .ne. 0 ) return

	call NeutralinoMasses(fail)
	if( fail .ne. 0 ) return

	call HiggsMassesConst(fail)
	end


************************************************************************

	subroutine ModelVarIni(fail, scale)
	implicit none
	integer fail
	RealType scale

#include "decl.h"

	call SMVarIni(fail, scale)
	if( fail .ne. 0 ) return

#ifndef NO_GUT_RELATION
* Gluino mass: GUT relation
* we use 1/Alfa(s) approx 1/Alfa(MZ2) - 20/(9*pi)*log(s/MZ2)
	Mino3 = Alfas*(1/AlfaMZ - 40/(9*pi)*log(scale/MZ))*SW2*Mino2
#endif

	MGl = abs(Mino3)
#ifndef NO_EXCLUSION_LIMITS
	if( MGl .lt. MGl_LOWERBOUND ) then
	  Error(fail, "MGl < MGl_LOWERBOUND")
	  return
	endif
#endif
	MGl2 = MGl**2

	SqrtEGl = Mino3
	if( MGl .gt. 1D-14 ) SqrtEGl = sqrt(SqrtEGl/MGl)

	call HiggsMassesVar(fail)
	end


************************************************************************

	subroutine ModelDigest
	implicit none

#include "decl.h"

	call SMDigest

	INFO "MSSM parameters"
	INFO "----------------------------"
	INFO "TB =", TB
	INFO "Mino1 =", Mino1
	INFO "Mino2 =", Mino2
	INFO "Mino3 =", Mino3
	INFO "MUE =", MUE
	INFO "MSusy =", MSusy
	INFO "Atau =", Af(2,3,3)
	INFO "At =", Af(3,3,3)
	INFO "Ab =", Af(4,3,3)
	INFO "----------------------------"
	INFO "Higgs masses:"
#if HIGGS_MASSES == FEYNHIGGS
	INFO "(using FeynHiggs)"
#elif HIGGS_MASSES == TREE
	INFO "(using tree-level Higgs masses)"
#elif HIGGS_MASSES == SIMPLE
	INFO "(using simple approximation)"
#else
	INFO "(using 2-loop approximation)"
#endif
	INFO "Mh0 =", Mh0
	INFO "MHH =", MHH
	INFO "MA0 =", MA0
	INFO "MHp =", MHp
	INFO "----------------------------"
	INFO "Sfermion masses:"
#ifdef NO_SQUARK_MIXING
	INFO "(diagonal sfermion mass matrix: L-R basis)"
#else
	INFO "(non-diagonal sfermion mass matrix: 1-2 basis)"
#endif
	INFO "Sneutrino  e:  ", MSf(1,1,1)
	INFO "           mu: ", MSf(1,1,2)
	INFO "           tau:", MSf(1,1,3)
	INFO "Slepton    e:  ", MSf(1,2,1), MSf(2,2,1)
	INFO "           mu: ", MSf(1,2,2), MSf(2,2,2)
	INFO "           tau:", MSf(1,2,3), MSf(2,2,3)
	INFO "u-Squark   u:  ", MSf(1,3,1), MSf(2,3,1)
	INFO "           c:  ", MSf(1,3,2), MSf(2,3,2)
	INFO "           t:  ", MSf(1,3,3), MSf(2,3,3)
	INFO "d-Squark   d:  ", MSf(1,4,1), MSf(2,4,1)
	INFO "           s:  ", MSf(1,4,2), MSf(2,4,2)
	INFO "           b:  ", MSf(1,4,3), MSf(2,4,3)
#ifdef NMFV
	INFO "u-Squark fv    ", MASf(1,3), MASf(2,3)
	INFO "               ", MASf(3,3), MASf(4,3)
	INFO "               ", MASf(5,3), MASf(6,3)
	INFO "d-Squark fv    ", MASf(1,4), MASf(2,4)
	INFO "               ", MASf(3,4), MASf(4,4)
	INFO "               ", MASf(5,4), MASf(6,4)
#endif
	INFO "----------------------------"
	INFO "Chargino masses:"
	INFO "MCha(1) =", MCha(1)
	INFO "MCha(2) =", MCha(2)
	INFO "----------------------------"
	INFO "Neutralino masses:"
	INFO "MNeu(1) =", MNeu(1)
	INFO "MNeu(2) =", MNeu(2)
	INFO "MNeu(3) =", MNeu(3)
	INFO "MNeu(4) =", MNeu(4)
	INFO "----------------------------"
	INFO "Gluino mass:"
	INFO "MGl =", MGl
	INFO "============================"
	end


*-----------------------------------------------------------------------
* Sfermion sector:

	subroutine SfermionMasses(fail)
	implicit none
	integer fail

#include "decl.h"

	integer t, g
	RealType i3part, qpart

	RealType Qe, Qu, Qd
	parameter (Qe = -1, Qu = 2/3D0, Qd = -1/3D0)

#ifndef NO_EXCLUSION_LIMITS
	RealType DeltaRho, f1f2
	external f1f2
#endif

#ifdef NMFV
	integer as, g2
	ComplexType MX2(6,6), m1, m2, m3

#define MX2_LL(i,j) MX2(i,j)
#define MX2_LR(i,j) MX2(i,j+3)
#define MX2_RL(i,j) Conjugate(MX2(j,i+3))
#define MX2_RR(i,j) MX2(i+3,j+3)

	ComplexType MX2_flat(6*6)
	equivalence (MX2, MX2_flat)
#endif

	i3part = .5D0*MZ2*C2B
	qpart = SW2*MZ2*C2B
	DSf(1,1) = +i3part
	DSf(2,1) = 0
	DSf(2,2) = Qe*qpart
	DSf(1,2) = -i3part - DSf(2,2)
	DSf(2,3) = Qu*qpart
	DSf(1,3) = +i3part - DSf(2,3)
	DSf(2,4) = Qd*qpart
	DSf(1,4) = -i3part - DSf(2,4)

	MUETB(3) = -Conjugate(MUE)/TB
	MUETB(2) = -Conjugate(MUE)*TB
	MUETB(4) = MUETB(2)

	do g = 1, 3
	  MSS(1,2,g) = MSL(g)
	  MSS2(1,2,g) = MSS(1,2,g)**2
	  MSS(2,2,g) = MSE(g)
	  MSS2(2,2,g) = MSS(2,2,g)**2
	  MSS(1,3,g) = MSQ(g)
	  MSS2(1,3,g) = MSS(1,3,g)**2
	  MSS(2,3,g) = MSU(g)
	  MSS2(2,3,g) = MSS(2,3,g)**2
	  MSS(1,4,g) = MSQ(g)
	  MSS2(1,4,g) = MSS(1,4,g)**2
	  MSS(2,4,g) = MSD(g)
	  MSS2(2,4,g) = MSS(2,4,g)**2

#ifdef NO_SQUARK_MIXING
#define Eigen2x2 NoMixing
	  Af(2,g,g) = -MUETB(2)
	  Af(3,g,g) = -MUETB(3)
	  Af(4,g,g) = -MUETB(4)
	  Xf(2,g) = 0
	  Xf(3,g) = 0
	  Xf(4,g) = 0
#else
	  Af(2,g,g) = Ae(g)
	  Xf(2,g) = Af(2,g,g) + MUETB(2)
	  Af(3,g,g) = Au(g)
	  Xf(3,g) = Af(3,g,g) + MUETB(3)
	  Af(4,g,g) = Ad(g)
	  Xf(4,g) = Af(4,g,g) + MUETB(4)
#endif

	  call NoMixing(fail, MSf(1,1,g), MSf2(1,1,g), USf(1,1,1,g),
     &      MSS2(1,2,g) + DSf(1,1), 1D300, DCMPLX(0D0))

	  do t = 2, 4
	    call Eigen2x2(fail, MSf(1,t,g), MSf2(1,t,g), USf(1,1,t,g),
     &        MSS2(1,t,g) + DSf(1,t) + Mf2(t,g),
     &        MSS2(2,t,g) + DSf(2,t) + Mf2(t,g),
     &        Mf(t,g)*Conjugate(Xf(t,g)))
	    if( fail .ne. 0 ) return

	    UCSf(1,1,t,g) = Sq(USf(1,1,t,g))
	    UCSf(2,1,t,g) = Sq(USf(2,1,t,g))
	    UCSf(1,2,t,g) = Sq(USf(1,2,t,g))
	    UCSf(2,2,t,g) = Sq(USf(2,2,t,g))

	    UCSf(3,1,t,g) = USf(1,1,t,g)*USfC(2,1,t,g)
	    UCSf(3,2,t,g) = USf(1,2,t,g)*USfC(2,2,t,g)
	    UCSf(1,3,t,g) = USf(1,1,t,g)*USfC(1,2,t,g)
	    UCSf(2,3,t,g) = USf(2,1,t,g)*USfC(2,2,t,g)
	    UCSf(3,3,t,g) = USf(1,1,t,g)*USfC(2,2,t,g)
	    UCSf(3,4,t,g) = USf(1,2,t,g)*USfC(2,1,t,g)

	    UUSf(1,1,t,g) = USf(1,1,t,g)**2
	    UUSf(2,1,t,g) = USf(2,1,t,g)**2
	    UUSf(1,2,t,g) = USf(1,2,t,g)**2
	    UUSf(2,2,t,g) = USf(2,2,t,g)**2

	    UUSf(3,1,t,g) = USf(1,1,t,g)*USf(2,1,t,g)
	    UUSf(3,2,t,g) = USf(1,2,t,g)*USf(2,2,t,g)
	    UUSf(1,3,t,g) = USf(1,1,t,g)*USf(1,2,t,g)
	    UUSf(2,3,t,g) = USf(2,1,t,g)*USf(2,2,t,g)
	    UUSf(3,3,t,g) = USf(1,1,t,g)*USf(2,2,t,g)
	    UUSf(3,4,t,g) = USf(1,2,t,g)*USf(2,1,t,g)
	  enddo
	enddo

#ifndef NO_EXCLUSION_LIMITS
	if( MSf(1,4,3) .lt. MSbottom_LOWERBOUND ) then
	  Error(fail, "MSf(1,4,3) < MSbottom_LOWERBOUND")
	  return
	endif

	if( MSf(1,3,3) .lt. MStop_LOWERBOUND ) then
	  Error(fail, "MSf(1,3,3) < MStop_LOWERBOUND")
	  return
	endif

	if( MSf(1,2,1) .lt. MSlepton_LOWERBOUND ) then
	  Error(fail, "MSf(1,2,1) < MSlepton_LOWERBOUND")
	  return
	endif

	if( MSf(1,3,1) .lt. MSquark_LOWERBOUND ) then
	  Error(fail, "MSf(1,3,1) < MSquark_LOWERBOUND")
	  return
	endif

	DeltaRho =
     &    USf2(1,1,3,3)*USf2(1,1,4,3)*f1f2(MSf2(1,3,3), MSf2(1,4,3)) +
     &    USf2(1,1,3,3)*USf2(2,1,4,3)*f1f2(MSf2(1,3,3), MSf2(2,4,3)) +
     &    USf2(2,1,3,3)*USf2(1,1,4,3)*f1f2(MSf2(2,3,3), MSf2(1,4,3)) +
     &    USf2(2,1,3,3)*USf2(2,1,4,3)*f1f2(MSf2(2,3,3), MSf2(2,4,3)) -
     &    USf2(1,2,3,3)*USf2(2,2,3,3)*f1f2(MSf2(1,3,3), MSf2(2,3,3)) -
     &    USf2(1,2,4,3)*USf2(2,2,4,3)*f1f2(MSf2(1,4,3), MSf2(2,4,3))

	if( DeltaRho .gt. DeltaRho_UPPERBOUND ) then
	  Error(fail, "DeltaRho > DeltaRho_UPPERBOUND")
	  return
	endif
#endif

#ifdef NMFV
* Note: we compute two sets of squark masses.
* 1) the "regular" ones with only left-right mixing,
* 2) the fully 6x6 mixed ones.
* The "regular" ones are needed for computing the Higgs mass.

	do t = 3, 4
	  MX2_LL(1,1) = MSS2(1,4,1)
	  MX2_LL(2,2) = MSS2(1,4,2)
	  MX2_LL(3,3) = MSS2(1,4,3)
	  MX2_LL(1,2) = MSS(1,4,1)*MSS(1,4,2)*deltaSf_LL(1,2,4)
	  MX2_LL(2,3) = MSS(1,4,2)*MSS(1,4,3)*deltaSf_LL(2,3,4)
	  MX2_LL(1,3) = MSS(1,4,1)*MSS(1,4,3)*deltaSf_LL(1,3,4)
	  if( t .eq. 3 ) then
	    do g2 = 1, 3
	      m1 = MX2_LL(1,1)*CKMC(g2,1) +
     &             MX2_LL(1,2)*CKMC(g2,2) +
     &             MX2_LL(1,3)*CKMC(g2,3)
	      m2 = Conjugate(MX2_LL(1,2))*CKMC(g2,1) +
     &             MX2_LL(2,2)*CKMC(g2,2) +
     &             MX2_LL(2,3)*CKMC(g2,3)
	      m3 = Conjugate(MX2_LL(1,3))*CKMC(g2,1) +
     &             Conjugate(MX2_LL(2,3))*CKMC(g2,2) +
     &             MX2_LL(3,3)*CKMC(g2,3)
	      do g = 1, g2
	        MX2_LL(g,g2) = CKM(g,1)*m1 + CKM(g,2)*m2 + CKM(g,3)*m3
	      enddo
	    enddo
	  endif
	  MX2_LL(1,1) = MX2_LL(1,1) + DSf(1,t) + Mf2(t,1)
	  MX2_LL(2,2) = MX2_LL(2,2) + DSf(1,t) + Mf2(t,2)
	  MX2_LL(3,3) = MX2_LL(3,3) + DSf(1,t) + Mf2(t,3)

	  MX2_RR(1,1) = MSS2(2,t,1) + DSf(2,t) + Mf2(t,1)
	  MX2_RR(2,2) = MSS2(2,t,2) + DSf(2,t) + Mf2(t,2)
	  MX2_RR(3,3) = MSS2(2,t,3) + DSf(2,t) + Mf2(t,3)
	  MX2_RR(1,2) = MSS(2,t,1)*MSS(2,t,2)*deltaSf_RR(1,2,t)
	  MX2_RR(2,3) = MSS(2,t,2)*MSS(2,t,3)*deltaSf_RR(2,3,t)
	  MX2_RR(1,3) = MSS(2,t,1)*MSS(2,t,3)*deltaSf_RR(1,3,t)

	  MX2_LR(1,1) = Mf(t,1)*Xf(t,1)
	  MX2_LR(2,2) = Mf(t,2)*Xf(t,2)
	  MX2_LR(3,3) = Mf(t,3)*Xf(t,3)
	  MX2_LR(1,2) = MSS(1,t,1)*MSS(2,t,2)*deltaSf_LR(1,2,t)
	  Af(t,1,2) = MX2_LR(1,2)/Mf(t,1)
	  MX2_LR(2,3) = MSS(1,t,2)*MSS(2,t,3)*deltaSf_LR(2,3,t)
	  Af(t,2,3) = MX2_LR(2,3)/Mf(t,2)
	  MX2_LR(1,3) = MSS(1,t,1)*MSS(2,t,3)*deltaSf_LR(1,3,t)
	  Af(t,1,3) = MX2_LR(1,3)/Mf(t,1)
	  MX2_LR(2,1) = Conjugate(MSS(2,t,1)*MSS(1,t,2)*deltaSf_RL(1,2,t))
	  Af(t,2,1) = MX2_LR(2,1)/Mf(t,2)
	  MX2_LR(3,2) = Conjugate(MSS(2,t,2)*MSS(1,t,3)*deltaSf_RL(2,3,t))
	  Af(t,3,2) = MX2_LR(3,2)/Mf(t,3)
	  MX2_LR(3,1) = Conjugate(MSS(2,t,1)*MSS(1,t,3)*deltaSf_RL(1,3,t))
	  Af(t,3,1) = MX2_LR(3,1)/Mf(t,3)

	  call HEigensystem(6, MX2,6, MASf2(1,t), UASf(1,1,t),6, 1)

	  if( MASf2(1,t) .lt. 0 ) then
	    Error(fail, "Negative MASf2")
	    return
	  endif

	  do as = 1, 6
	    MASf(as,t) = sqrt(MASf2(as,t))
	  enddo
	enddo
#endif
	end


*-----------------------------------------------------------------------
* this function is needed for the calculation of Delta rho

	RealType function f1f2(m1, m2)
	implicit none
	RealType m1, m2

#include "decl.h"

	RealType f1, f2, logm1m2, c1

	logm1m2 = log(m1/m2)
	c1 = 2*m1*m2/(m1 - m2)*logm1m2
	f1 = m1 + m2 - c1
	f2 = m1 + m2 - c1*(2 + m1/m2*logm1m2)
     &   + (m1 + m2)*m1**2/(m1 - m2)**2 * logm1m2**2
     &   - 2*(m1 - m2) * Li2(1 - m1/m2)

	f1f2 = 3*GF/(8*sqrt2*pi**2) * (f1 + 2/(3*pi)*AlfasMT*f2)
	end


*-----------------------------------------------------------------------
* Chargino sector:

	subroutine CharginoMasses(fail)
	implicit none
	integer fail

#include "decl.h"

	ComplexType X(2,2)

	X(1,1) = Mino2
	X(2,1) = sqrt2*MW*CB
	X(1,2) = sqrt2*MW*SB
	X(2,2) = MUE

	call SVD(2, 2, X,2, MCha, UCha,2, VCha,2, 1)

#ifndef NO_EXCLUSION_LIMITS
	if( MCha(1) .lt. MCha_LOWERBOUND ) then
	  Error(fail, "MCha(1) < MCha_LOWERBOUND")
	  return
	endif
#endif

	MCha2(1) = MCha(1)**2
	MCha2(2) = MCha(2)**2
	end


*-----------------------------------------------------------------------
* Neutralino sector:

	subroutine NeutralinoMasses(fail)
	implicit none
	integer fail

#include "decl.h"

	ComplexType Y(4,4)

#ifndef NO_GUT_RELATION
* GUT relation (Higgs Hunters Guide p. 287)
	Mino1 = 5/3D0*SW2/CW2*Mino2
#endif

	Y(1,1) = Mino1
	Y(2,2) = Mino2
	Y(3,1) = -MZ*SW*CB
	Y(1,3) = Y(3,1)
	Y(4,1) = MZ*SW*SB
	Y(1,4) = Y(4,1)
	Y(3,2) = MZ*CW*CB
	Y(2,3) = Y(3,2)
	Y(4,2) = -MZ*CW*SB
	Y(2,4) = Y(4,2)
	Y(4,3) = -MUE
	Y(3,4) = Y(4,3)
	Y(2,1) = 0
	Y(1,2) = 0
	Y(3,3) = 0
	Y(4,4) = 0

	call TakagiFactor(4, Y,4, MNeu, ZNeu,4, 1)

#ifndef NO_EXCLUSION_LIMITS
	if( MNeu(1) .lt. MNeu_LOWERBOUND ) then
	  Error(fail, "MNeu(1) < MNeu_LOWERBOUND")
	  return
	endif
#endif

	MNeu2(1) = MNeu(1)**2
	MNeu2(2) = MNeu(2)**2
	MNeu2(3) = MNeu(3)**2
	MNeu2(4) = MNeu(4)**2
	end


*-----------------------------------------------------------------------
* Higgs sector:

	subroutine HiggsMassesConst(fail)
	implicit none
	integer fail

#include "decl.h"

	RealType avgMH2, deltaMH2
	integer inputmass
	data inputmass /0/

	if( inputmass .eq. 0 ) then
	  if( MA0 .ge. 0 ) then
	    inputmass = 3
	  else if( MHp .gt. MW ) then
	    inputmass = 4
	  else
	    Error(fail, "Neither MA0 nor MHp given")
	    return
	  endif
	endif

	if( inputmass .eq. 3 ) then
	  MA0tree2 = MA0**2
	  MHptree2 = MA0tree2 + MW2
	else if( inputmass .eq. 4 ) then
	  MA0tree2 = (MHp - MW)*(MHp + MW)
	  MHptree2 = MHp**2
	endif

	avgMH2 = .5D0*(MA0tree2 + MZ2)
	deltaMH2 = sqrt(MA0tree2*MZ2)*C2B
	deltaMH2 = sqrt((avgMH2 - deltaMH2)*(avgMH2 + deltaMH2))
	Mh0tree2 = avgMH2 - deltaMH2
	MHHtree2 = avgMH2 + deltaMH2

	CA2 = .5D0*(SB2*MA0tree2 + CB2*MZ2 - Mh0tree2)/deltaMH2
	CA = sqrt(CA2)
	SA = -SB*CB*avgMH2/(CA*deltaMH2)

#if HIGGS_MASSES != FEYNHIGGS
	Mh02 = Mh0tree2
	MHH2 = MHHtree2
	MA02 = MA0tree2
	MHp2 = MHptree2
#if HIGGS_MASSES != TREE
	call HiggsCorr(fail)
	if( fail .ne. 0 ) return
#endif
	Mh0 = sqrt(Mh02)
	MHH = sqrt(MHH2)
	MA0 = sqrt(MA02)
	MHp = sqrt(MHp2)
#endif

	SA2 = SA**2
	S2A = 2*SA*CA
	C2A = (CA - SA)*(CA + SA)
	SAB = SA*CB + CA*SB
	CAB = CA*CB - SA*SB
	SBA = SB*CA - CB*SA
	SBA2 = SBA**2
	CBA = CB*CA + SB*SA
	CBA2 = CBA**2
	end


*-----------------------------------------------------------------------
* Higgs sector: some parts need to be evaluated in the `Var' part,
* e.g. FeynHiggs needs the gluino mass, which is known only after
* alpha_s(s) is set.

	subroutine HiggsMassesVar(fail)
	implicit none
	integer fail

#include "decl.h"

#if HIGGS_MASSES == FEYNHIGGS
	call FeynHiggs(fail)
	if( fail .ne. 0 ) return
#endif

#ifndef NO_EXCLUSION_LIMITS
	if( Mh0 .lt. Mh0_LOWERBOUND ) then
	  Error(fail, "Mh0 < Mh0_LOWERBOUND")
	  return
	endif
#endif

	MH = Mh0
	MH2 = Mh02

#ifdef SM_HIGGS
* With SM_HIGGS defined, the values of the MSSM parameters are set
* such that the SM is reproduced.  See Higgs Hunters Guide p. 356.
* Note one tricky thing: the Higgs mass Mh0 is first calculated as
* it would be in the MSSM (to get the same value when comparing
* MSSM <-> SM), but then beta is reverse-engineered to be consistent
* with the SM relations.  In other words, we start with the MSSM-
* value of beta, calculate Mh0, and then redefine beta to fulfill
* cos(2 beta) = Mh0/MZ.  In principle, beta could also be fixed such
* that cos(2 beta) = -Mh0/MZ, but this is insignificant here since
* we only need beta to be consistent with the values of Mh0 and MZ.

	C2B = Mh0/MZ
	if( C2B .gt. 1 ) then
	  fail = 4
	  return
	endif
	CB2 = .5D0*(1 + C2B)
	CB = sqrt(CB2)
	SB2 = (1 - CB)*(1 + CB)
	TB2 = SB2/CB2
	TB = sqrt(TB2)
	SB = CB*TB
	S2B = 2*TB*CB2

	SA = -CB
	CA = SB
	SA2 = CB2
	CA2 = SB2
	S2A = -S2B
	C2A = -C2B
	SAB = -C2B
	CAB = S2B
	SBA = 1
	SBA2 = 1
	CBA = 0
	CBA2 = 0
#endif
	end


*-----------------------------------------------------------------------

#if HIGGS_MASSES != FEYNHIGGS && HIGGS_MASSES != TREE

	subroutine HiggsCorr(fail)
	implicit none
	integer fail

#include "decl.h"

	RealType MTr2, Qsusy, Lam
	RealType Xt, topstop, toplog, mix
	RealType Pb, Pf, Pg, Pg_, P1h, P2h, P2h_
	RealType pref, pref2, tmp, tmplog
	RealType sig11, sig12, sig22
	RealType avgMH2, deltaMH2

	RealType MtLR, MbLR
	equivalence (Xf(3,3), MtLR)
	equivalence (Xf(4,3), MbLR)

	integer Nc, Ng
	RealType Qt, Qb
	parameter (Nc = 3, Ng = 3, Qt = 2/3D0, Qb = -1/3D0)

* this is the running top mass:
	MTr2 = MT2/(1 + 4/(3*pi)*AlfasMT)**2

#define SusyScale(mt2) sqrt((MSS2(1,3,3) + mt2)*(MSS2(2,3,3) + mt2))

#if HIGGS_MASSES == SIMPLE

* Simple approximation: (for complex and/or real case)
	Qsusy = SusyScale(MTr2)
	sig22 = -3*GF/(sqrt2*pi**2)*MTr2**2/SB2*log(Qsusy/MTr2)

#else

* FeynHiggsFast approximation: hep-ph/9903404 (only for real case)
	Qsusy = SusyScale(MT2)
	mix = MtLR**2/Qsusy

* The one-loop top/stop contributions:
	Lam = 1/8D0 - SW2*(1/3D0 - 4/9D0*SW2)
	topstop = MT2/Qsusy
	toplog = log(topstop)
	tmp = MZ2/MT2
	sig11 = GF*sqrt2/pi**2*MZ2**2*Lam*CB2 * toplog
	sig12 = -GF*sqrt2/pi**2*MZ2/TB *
     &    (-3/8D0*MT2 + MZ2*Lam*SB2) * toplog
	sig22 = GF*sqrt2/(8*pi**2)*MT2**2/SB2 *
     &    ( tmp*(-2 + 11/10D0*tmp)
     &       + (12 + SB2*tmp*(-6 + 8*Lam*SB2*tmp)) * toplog
     &       + mix*(-12 + 4*tmp + 6*topstop
     &         + mix*(1 + topstop*(-4 + 3*topstop)
     &           + mix*topstop
     &             *(3/5D0 + topstop*(-12/5D0 + 2*topstop)
     &             + mix*topstop
     &               *(3/7D0 + topstop*(-12/7D0 + 3/2D0*topstop))
     &       ))) )

* The dominant two-loop top/stop contributions:
	toplog = log(MTr2/Qsusy)
	tmp = (MSf2(2,3,3) - MSf2(1,3,3))/MTr2
     &    * (USf(1,1,3,3)*USf(2,1,3,3))**2
	tmplog = log(MSf2(2,3,3)/MSf2(1,3,3))
	Xt = tmp*(
     &    tmp*(2 - (MSf2(2,3,3) + MSf2(1,3,3))/
     &             (MSf2(2,3,3) - MSf2(1,3,3))*tmplog)
     &    + 2*tmplog)
	tmplog = log(MSf(2,3,3)*MSf(1,3,3)/MTr2)
	sig22 = sig22 + MTr2**2/SB2 * (
     &    + GF*sqrt2/pi**3*AlfasMT*(toplog*(3*toplog - 6)
     &      - 6*sign(sqrt(mix), MtLR) + mix*(-3*toplog + 3/4D0*mix))
     &    - 9*GF**2/(16*pi**4)*MTr2* tmplog*(Xt + tmplog) )

* The one-loop leading-log contributions from all other sectors:
* [NB: MSQ is used here instead of MSusy.]
	Pb = Nc*(1 + 4*Qb*SW2*(1 + 2*Qb*SW2))
	Pf = Nc*(Ng - 1)*(2 - SW2*(4 - 8*(Qt**2 + Qb**2)*SW2))
     &    + Ng*(2 - SW2*(4 - 8*SW2))
	Pg = -44 + SW2*(106 - 62*SW2)
	Pg_ = 10 + SW2*(34 - 26*SW2)
	P1h = (1 - 2*SW2*CW2 - 9*C2B**2)*C2B**2
	P2h = -10 + 2*SW2*CW2
	P2h_ = 8 - SW2*(22 - 10*SW2)
	tmplog = log(MSS2(1,3,3)/MZ2)
	tmp = MB2/(MZ2*CB2)
	pref = GF*MZ2**2/(12*sqrt2*pi**2)
	pref2 = GF*Nc/(4*sqrt2*pi**2)*MB2/MSS2(1,3,3)
	sig11 = sig11
     &    - pref*CB2
     &      *(6*Nc*tmp*(2*tmp - 1) + Pb + Pf + Pg + P2h)*tmplog
     &    - pref2*Af(4,3,3)
     &      *(4*MB2*MbLR/CB2*(1 - Af(4,3,3)*MbLR/(12*MSS2(1,3,3)))
     &        - MZ2*(MbLR + Af(4,3,3)/3D0))
	sig12 = sig12
     &    - pref*SB*CB*(3*Nc*tmp - Pb - Pf - Pg_ - P2h_)*tmplog
     &    + .5D0*pref2
     &      *(4*MB2*MUE*MbLR/CB2*(1 - Af(4,3,3)*MbLR/(6*MSS2(1,3,3)))
     &        - MZ2*TB*(MbLR*(Af(4,3,3) + MUE/TB)
     &                  + (MUE**2 + Af(4,3,3)**2)/3D0))
	sig22 = sig22
     &    - pref*SB2*(Pb + Pf + Pg + P2h)*tmplog
     &    + pref2*MUE
     &      *(MB2*MUE*MbLR**2/(3*CB2*MSS2(1,3,3))
     &        + MZ2*TB*(MbLR + MUE*TB/3D0))
	if( MA02 .gt. MZ2 ) then
	  tmplog = log(MA02/MZ2)
	  sig11 = sig11 - pref*CB2*(P1h - P2h)*tmplog
	  sig12 = sig12 - pref*SB*CB*(P1h + P2h_)*tmplog
	  sig22 = sig22 - pref*SB2*(P1h - P2h)*tmplog
	endif
#endif

	avgMH2 = .5D0*(MA02 + MZ2 - sig22 - sig11)
	deltaMH2 = .5D0*sqrt(
     &    (C2B*(MA02 - MZ2) + sig11 - sig22)**2 +
     &    (S2B*(MA02 + MZ2) + 2*sig12)**2 )
	Mh02 = avgMH2 - deltaMH2
	MHH2 = avgMH2 + deltaMH2

	if( Mh02 .lt. 0 ) then
	  Error(fail, "Negative Mh0^2")
	  return
	endif

	CA2 = .5D0*(SB2*MA02 + CB2*MZ2 - Mh02 - sig11)/deltaMH2
	CA = sqrt(CA2)
	SA = -.5D0*(SB*CB*(MA02 + MZ2) + sig12)/(CA*deltaMH2)

	if( MA0 .le. 250 .and. MSS(1,3,3) .ge. 250 ) then
	  MHp2 = MHp2 +
     &      Alfa*MW2/(12*pi*SW2)*
     &        (Nc*(Ng - 1) + Ng - 9 + 15*SW2/CW2)*
     &        log(MSS2(1,3,3)/MW2) +
     &      Nc*Alfa/(8*pi*SW2*MW2)*
     &        (2*MTr2*MB2/(SB2*CB2)
     &          - MW2*(MTr2/SB2 + MB2/CB2 - 2/3D0*MW2))*
     &        log(MSS2(1,3,3)/MTr2)
	endif
	end

#endif


*-----------------------------------------------------------------------

#if HIGGS_MASSES == FEYNHIGGS

	subroutine FeynHiggs(fail)
	implicit none
	integer fail

#include "decl.h"

	integer mssmpart, fieldren, tanbren, higgsmix, p2approx
	integer looplevel, runningMT, botResum, tlCplxApprox
	parameter (mssmpart = 4)
	parameter (fieldren = 0)
	parameter (tanbren = 0)
	parameter (higgsmix = 2)
	parameter (p2approx = 0)
	parameter (looplevel = 2)
	parameter (runningMT = 1)
	parameter (botResum = 1)
	parameter (tlCplxApprox = 3)

	RealType scalefactor
	RealType Qtau, Qt, Qb
	parameter (scalefactor = 1)
	parameter (Qtau = 0)
	parameter (Qt = 0)
	parameter (Qb = 0)

	ComplexType SAeff

	call FHSetFlags(fail, mssmpart, fieldren, tanbren, higgsmix,
     &    p2approx, looplevel, runningMT, botResum, tlCplxApprox)
	if( fail .ne. 0 ) return

	call FHSetSMPara(fail,
     &    1/Alfa, AlfasMZ, GF,
     &    ME, MU, MD, MM, MC, MS, ML, MBatMB,
     &    MW, MZ,
#ifdef NO_CKM
     &    0D0, 0D0, 0D0, 0D0)
#else
     &    CKMlambda, CKMA, CKMrhobar, CKMetabar)
#endif
	if( fail .ne. 0 ) return

	call FHSetPara(fail, scalefactor,
     &    MT, TB, MA0, MHp,
     &    MSS(1,2,3), MSS(2,2,3), MSS(1,3,3), MSS(2,3,3), MSS(2,4,3),
     &    MSS(1,2,2), MSS(2,2,2), MSS(1,3,2), MSS(2,3,2), MSS(2,4,2),
     &    MSS(1,2,1), MSS(2,2,1), MSS(1,3,1), MSS(2,3,1), MSS(2,4,1),
     &    MUE,
     &    Af(2,3,3), Af(3,3,3), Af(4,3,3),
     &    Af(2,2,2), Af(3,2,2), Af(4,2,2),
     &    Af(2,1,1), Af(3,1,1), Af(4,1,1),
     &    Mino1, Mino2, Mino3,
     &    Qtau, Qt, Qb)
        if( fail .ne. 0 ) return

#ifdef NMFV
	call FHSetNMFV(fail,
     &    deltaLL12, deltaLL23, deltaLL13,
     &    deltaLRuc, deltaLRct, deltaLRut,
     &    deltaRLuc, deltaRLct, deltaRLut,
     &    deltaRRuc, deltaRRct, deltaRRut,
     &    deltaLRds, deltaLRsb, deltaLRdb,
     &    deltaRLds, deltaRLsb, deltaRLdb,
     &    deltaRRds, deltaRRsb, deltaRRdb)
	if( fail .ne. 0 ) return
#endif

	call FHHiggsCorr(fail, MHiggs, SAeff, UHiggs, ZHiggs)

	MHiggs2(1) = MHiggs(1)**2
	MHiggs2(2) = MHiggs(2)**2
	MHiggs2(3) = MHiggs(3)**2
	MHiggs2(4) = MHiggs(4)**2
	end

#endif

