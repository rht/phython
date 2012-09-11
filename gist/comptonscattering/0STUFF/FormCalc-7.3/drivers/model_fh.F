* model_fh.F
* parameter initialization via FeynHiggs
* this file is part of FormCalc
* last modified 30 Nov 11 th


* NO_EXCLUSION_LIMITS
*   If defined, the experimental exclusion limits are ignored.


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


************************************************************************

	subroutine ModelDefaults(argc, argv)
	implicit none
	integer argc
	character*128 argv(*)

#include "decl.h"
#include "SLHA.h"

	ComplexType slhadata(nslhadata)
	integer error

	character*16 flags
	data flags /"400202113"/

	if( argc .lt. 1 .or. argc .gt. 2 )
     &    stop "model_fh: need arguments <file> [<flags>]"

	if( argc .ge. 2 ) flags = argv(2)
	call FHSetFlagsString(error, flags)
	if( error .ne. 0 ) stop

	call FHReadRecord(error, record, slhadata, argv(1))
* 0 = SLHA, 2 = FH file, all others are true errors
	if( ibclr(error, 1) .ne. 0 ) stop

* turn on verbose mode:   
	record(iAdmin,iStep) = 1
	end


************************************************************************

	subroutine ModelConstIni(fail)
	implicit none
	integer fail

#include "decl.h"

	RealType invAlfaMZ, scalefactor, asMB
	integer i, nf
#ifndef NO_EXCLUSION_LIMITS
	RealType gm2, DeltaRho
	RealType MWMSSM, MWSM, SW2MSSM, SW2SM
	RealType edmeTh, edmn, edmHg
#endif

	call FHLoopRecord(fail, record)
	if( fail .ne. 0 ) return

	call FHSetRecord(fail, record)
	if( fail .ne. 0 ) return

	call FHRetrieveSMPara(fail,
     &    invAlfaMZ, AlfasMZ, GF,
     &    ME, MU, MD, MM, MC, MS, ML, MBatMB,
     &    MW, MZ,
     &    CKMlambda, CKMA, CKMrhobar, CKMetabar)
	if( fail .ne. 0 ) return

	call FHGetSMPara(fail, CKM)
	if( fail .ne. 0 ) return

	call FHRetrieveOSPara(fail, scalefactor,
     &    MT, TB, MA0, MHp,
     &    MSS(1,3), MSS(2,3), MSS(3,3), MSS(4,3), MSS(5,3),
     &    MSS(1,2), MSS(2,2), MSS(3,2), MSS(4,2), MSS(5,2),
     &    MSS(1,1), MSS(2,1), MSS(3,1), MSS(4,1), MSS(5,1),
     &    MUE,
     &    Afd(2,3), Afd(3,3), Afd(4,3),
     &    Afd(2,2), Afd(3,2), Afd(4,2),
     &    Afd(2,1), Afd(3,1), Afd(4,1),
     &    Mino1, Mino2, Mino3)
	if( fail .ne. 0 ) return

	call FHGetPara(fail, nmfv, MASf, UASf,
     &    MCha, UCha, VCha, MNeu, ZNeu, Deltab, MGl,
     &    MHtree, SA)
	if( fail .ne. 0 ) return

	call FHGetNMFV(fail, MSS2, Kf)
	if( fail .ne. 0 ) return

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

	if( MCha(1) .lt. MCha_LOWERBOUND ) then
	  Error(fail, "MCha(1) < MCha_LOWERBOUND")
	  return
	endif

	if( MNeu(1) .lt. MNeu_LOWERBOUND ) then
	  Error(fail, "MNeu(1) < MNeu_LOWERBOUND")
	  return
	endif
#endif

	call FHHiggsCorr(fail, MHiggs, SAeff, UHiggs, ZHiggs)
	if( fail .ne. 0 ) return

#ifndef NO_EXCLUSION_LIMITS
	if( Mh0 .lt. Mh0_LOWERBOUND ) then
	  Error(fail, "Mh0 < Mh0_LOWERBOUND")
	  return
	endif

	call FHConstraints(fail, gm2, DeltaRho,
     &    MWMSSM, MWSM, SW2MSSM, SW2SM,
     &    edmeTh, edmn, edmHg)

	if( DeltaRho .gt. DeltaRho_UPPERBOUND ) then
	  Error(fail, "DeltaRho > DeltaRho_UPPERBOUND")
	  return
	endif
#endif

	ME2 = ME**2
	MM2 = MM**2
	ML2 = ML**2
	MU2 = MU**2
	MC2 = MC**2
	MT2 = MT**2
	MD2 = MD**2
	MS2 = MS**2

	call AlphaS(asMB, nf, MBatMB)
	MB = MBatMB*sqrt(1 + 4*2/(3*pi)*asMB)
	MB2 = MB**2

	PRINT *, "MBatMB=", MBatMB
	PRINT *, "MBosMS=", MB
	PRINT *, "MBosDR=", MBatMB*sqrt(1 + 5*2/(3*pi)*asMB)
	PRINT *, "MBosDR=", MBatMB*(1 + 4/(3*pi)*asMB)

	MW2 = MW**2
	MZ2 = MZ**2
	CW = MW/MZ
	CW2 = CW**2
	SW2 = (1 - CW)*(1 + CW)
	SW = sqrt(SW2)

	AlfaGF = sqrt2/pi*GF*MW2*SW2

	AlfaMZ = 1/invAlfaMZ
	ELMZ = sqrt(4*pi*AlfaMZ)

	Alfa = Alfa0
	Alfa2 = Alfa**2
	EL = sqrt(4*pi*Alfa)

	Mh0tree2 = Mh0tree**2
	MHHtree2 = MHHtree**2
	MA0tree2 = MA0tree**2
	MHptree2 = MHptree**2

	MH = Mh0
	MH2 = MH**2

	Mh02 = Mh0**2
	MHH2 = MHH**2
	MA02 = MA0**2
	MHp2 = MHp**2

	do i = 1, 6*4
	  MASf2_flat(i) = MASf_flat(i)**2
	enddo

	MCha2(1) = MCha(1)**2
	MCha2(2) = MCha(2)**2

	MNeu2(1) = MNeu(1)**2
	MNeu2(2) = MNeu(2)**2
	MNeu2(3) = MNeu(3)**2
	MNeu2(4) = MNeu(4)**2

	TB2 = TB**2
	CB2 = 1/(1 + TB2)
	SB2 = CB2*TB2
	CB = sqrt(CB2)
	SB = CB*TB
	S2B = 2*TB*CB2
	C2B = (CB - SB)*(CB + SB)

	SA2 = SA**2  
	CA2 = (1 - SA)*(1 + SA)
	CA  = sqrt(CA2)
	S2A = 2*SA*CA
	C2A = (CA - SA)*(CA + SA)
	SAB = SA*CB + CA*SB
	CAB = CA*CB - SA*SB
	SBA = SB*CA - CB*SA
	SBA2 = SBA**2 
	CBA = CB*CA + SB*SA
	CBA2 = CBA**2
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

	integer error

	call FHOutput(error, "stdout", 7, -1D0)
	end

