* model_sm.h
* declarations for model_sm.F
* this file is part of FormCalc
* last modified 30 Nov 11 th


	RealType CKMlambda, CKMA, CKMrhobar, CKMetabar
	parameter (CKMlambda = .2253D0)
	parameter (CKMA = .808D0)
	parameter (CKMrhobar = .132D0)
	parameter (CKMetabar = .341D0)

	RealType MZ, MZ2, MW, MW2, CW, CW2, SW2
	parameter (MZ = 91.1876D0, MZ2 = MZ**2)
	parameter (MW = 80.399D0, MW2 = MW**2)
	parameter (CW = MW/MZ, CW2 = CW**2)
	parameter (SW2 = 1 - CW2)

	RealType GF, Alfa, Alfa2, AlfaMZ, AlfasMZ
	parameter (GF = 1.16637D-5)
	parameter (Alfa = 1/137.035999679D0, Alfa2 = Alfa**2)
c	parameter (Alfa = sqrt2/pi*GF*MW2*SW2, Alfa2 = Alfa**2)
	parameter (AlfaMZ = 1/127.934D0)
	parameter (AlfasMZ = .1184D0)

	RealType ME, ME2, MM, MM2, ML, ML2
	parameter (ME = .510998910D-3, ME2 = ME**2)
	parameter (MM = 105.658367D-3, MM2 = MM**2)
	parameter (ML = 1776.82D-3, ML2 = ML**2)

	RealType MU, MU2, MC, MC2, MT, MT2
	parameter (MU = 53.8D-3, MU2 = MU**2)
	parameter (MC = 1.27D0, MC2 = MC**2)
	parameter (MT = 172D0, MT2 = MT**2)

	RealType MD, MD2, MS, MS2, MB, MB2, MBatMB
	parameter (MD = 53.8D-3, MD2 = MD**2)
	parameter (MS = 101D-3, MS2 = MS**2)
	parameter (MB = 4.67D0, MB2 = MB**2)
	parameter (MBatMB = 4.25D0)

	ComplexType CKM(3,3)
	RealType Mf(4,3), Mf2(4,3)
	RealType MH, MH2
	RealType EL, GS, Alfas, Alfas2, SW

	common /sm_para/ CKM
	common /sm_para/ Mf, Mf2
	common /sm_para/ MH, MH2
	common /sm_para/ EL, GS, Alfas, Alfas2, SW

#ifndef CKMC
#define CKMC(i,j) Conjugate(CKM(i,j))
#endif

