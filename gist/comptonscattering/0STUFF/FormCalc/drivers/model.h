* model.h
* common blocks for the model parameters
* this file is part of FormCalc
* last modified 18 Oct 09 th


	double precision pi, degree, sqrt2, hbar_c2

	parameter (pi = 3.1415926535897932384626433832795029D0)
	parameter (degree = pi/180D0)
	parameter (sqrt2 = 1.41421356237309504880168872421D0)

	parameter (hbar_c2 = 3.8937966D8)
*         = hbar c^2 in picobarn

	double complex bogus, cI

	parameter (bogus = (-1D123, -2D123))
*	  some weird number likely to noticeably distort the final result;
*	  used for initializing arrays to check that all components
*	  have been calculated

	parameter (cI = (0D0, 1D0))

	double precision Divergence
	common /renorm/ Divergence


* SM parameters

	double precision sin12, sin23, sin13, delta13
	parameter (sin12 = .2243D0)
	parameter (sin23 = .0413D0)
	parameter (sin13 = .0037D0)
	parameter (delta13 = 60*degree)

	double precision MZ, MZ2, MW, MW2, CW, CW2, SW2
	parameter (MZ = 91.1875D0, MZ2 = MZ**2)
	parameter (MW = 80.450D0, MW2 = MW**2)
	parameter (CW = MW/MZ, CW2 = CW**2)
	parameter (SW2 = 1 - CW2)

	double precision GF, Alfa, Alfa2, AlfaMZ, AlfasMZ
	parameter (GF = 1.16639D-5)
	parameter (Alfa = 1/137.0359895D0, Alfa2 = Alfa**2)
c	parameter (Alfa = sqrt2/pi*GF*MW2*SW2, Alfa2 = Alfa**2)
	parameter (AlfaMZ = 1/127.934D0)
	parameter (AlfasMZ = .118D0)

	double precision ME, ME2, MM, MM2, ML, ML2
	parameter (ME = .51099907D-3, ME2 = ME**2)
	parameter (MM = 105.658389D-3, MM2 = MM**2)
	parameter (ML = 1777D-3, ML2 = ML**2)

	double precision MU, MU2, MC, MC2, MT, MT2
	parameter (MU = 53.8D-3, MU2 = MU**2)
	parameter (MC = 1.50D0, MC2 = MC**2)
	parameter (MT = 170.9D0, MT2 = MT**2)

	double precision MD, MD2, MS, MS2, MB, MB2, MBatMB
	parameter (MD = 53.8D-3, MD2 = MD**2)
	parameter (MS = 150D-3, MS2 = MS**2)
	parameter (MB = 4.7D0, MB2 = MB**2)
	parameter (MBatMB = 4.25D0)

	double complex CKM(3,3), CKMrhoeta
	double precision CKMlambda, CKMA
	double precision Mf(4,3), Mf2(4,3)
	double precision MH, MH2
	double precision EL, GS, Alfas, Alfas2, AlfasMT, SW
	logical sm_digest

	common /sm_para/ CKM, CKMrhoeta, CKMlambda, CKMA
	common /sm_para/ Mf, Mf2
	common /sm_para/ MH, MH2
	common /sm_para/ EL, GS, Alfas, Alfas2, AlfasMT, SW
	common /sm_para/ sm_digest

#ifndef CKMC
#define CKMC(i,j) DCONJG(CKM(i,j))
#endif


* MSSM parameters

	double complex UCha(2,2), VCha(2,2), ZNeu(4,4)
	double complex USf(2,2,4,3), UCSf(3,4,2:4,3), UUSf(3,4,2:4,3)
	double complex XHiggs(3,3,2)
	double complex Af(2:4,3,3), Xf(2:4,3,3)
	double complex Atau, At, Ab, MUE
	double complex Mino1, Mino2, Mino3, SqrtEGl
	double precision MCha(2), MCha2(2), MNeu(4), MNeu2(4)
	double precision MSS(2,2:4,3), MSS2(2,2:4,3), DSf(2,4)
	double precision MSf(2,4,3), MSf2(2,4,3), MSusy, MGl, MGl2
	double precision MHiggs(4), MHiggs2(4), MHiggstree2(4)
	double precision CB, SB, TB, CB2, SB2, TB2, C2B, S2B
	double precision CA, SA, CA2, SA2, C2A, S2A
	double precision CAB, SAB, CBA, SBA, CBA2, SBA2
	logical mssm_digest

	common /mssm_para/ UCha, VCha, ZNeu
	common /mssm_para/ USf, UCSf, UUSf
	common /mssm_para/ XHiggs
	common /mssm_para/ Af, Xf
	common /mssm_para/ Atau, At, Ab, MUE
	common /mssm_para/ Mino1, Mino2, Mino3, SqrtEGl
	common /mssm_para/ MCha, MCha2, MNeu, MNeu2
	common /mssm_para/ MSS, MSS2, DSf
	common /mssm_para/ MSf, MSf2, MSusy, MGl, MGl2
	common /mssm_para/ MHiggs, MHiggs2, MHiggstree2
	common /mssm_para/ CB, SB, TB, CB2, SB2, TB2, C2B, S2B
	common /mssm_para/ CA, SA, CA2, SA2, C2A, S2A
	common /mssm_para/ CAB, SAB, CBA, SBA, CBA2, SBA2
	common /mssm_para/ mssm_digest

#ifndef USfC
#define USfC(i,j,t,g) DCONJG(USf(i,j,t,g))
#define UCSfC(i,j,t,g) DCONJG(UCSf(i,j,t,g))
#define UUSfC(i,j,t,g) DCONJG(UUSf(i,j,t,g))
#define USf2(i,j,t,g) DBLE(UCSf(i,j,t,g))
#define VChaC(i,j) DCONJG(VCha(i,j))
#define UChaC(i,j) DCONJG(UCha(i,j))
#define ZNeuC(i,j) DCONJG(ZNeu(i,j))
#define UHiggsC(i,j) DCONJG(UHiggs(i,j))
#define ZHiggsC(i,j) DCONJG(ZHiggs(i,j))
#define AfC(t,g1,g2) DCONJG(Af(t,g1,g2))
#define Mino3C DCONJG(Mino3)
#define MUEC DCONJG(MUE)
#define SqrtEGlC DCONJG(SqrtEGl)
#endif

	double precision Mh0, Mh02, MHH, MHH2, MA0, MA02, MHp, MHp2
	equivalence (MHiggs(1), Mh0), (MHiggs2(1), Mh02)
	equivalence (MHiggs(2), MHH), (MHiggs2(2), MHH2)
	equivalence (MHiggs(3), MA0), (MHiggs2(3), MA02)
	equivalence (MHiggs(4), MHp), (MHiggs2(4), MHp2)

	double precision Mh0tree2, MHHtree2, MA0tree2, MHptree2
	equivalence (MHiggstree2(1), Mh0tree2)
	equivalence (MHiggstree2(2), MHHtree2)
	equivalence (MHiggstree2(3), MA0tree2)
	equivalence (MHiggstree2(4), MHptree2)

	double complex UHiggs(3,3), ZHiggs(3,3)
	equivalence (XHiggs(1,1,1), UHiggs)
	equivalence (XHiggs(1,1,2), ZHiggs)

	double precision Af_flat(3*3*3), Xf_flat(3*3*3)
	equivalence (Af, Af_flat)
	equivalence (Xf, Xf_flat)

	double precision ReImAtau(2), ReAtau, ImAtau
	equivalence (Atau, ReImAtau)
	equivalence (ReImAtau(1), ReAtau), (ReImAtau(2), ImAtau)

	double precision ReImAt(2), ReAt, ImAt
	equivalence (At, ReImAt)
	equivalence (ReImAt(1), ReAt), (ReImAt(2), ImAt)

	double precision ReImAb(2), ReAb, ImAb
	equivalence (Ab, ReImAb)
	equivalence (ReImAb(1), ReAb), (ReImAb(2), ImAb)

	double precision ReImMUE(2), ReMUE, ImMUE
	equivalence (MUE, ReImMUE)
	equivalence (ReImMUE(1), ReMUE), (ReImMUE(2), ImMUE)

	double precision ReImMino1(2), ReMino1, ImMino1
	double precision M_1, ReM_1, ImM_1
	equivalence (Mino1, M_1, ReImMino1)
	equivalence (ReImMino1(1), ReMino1, ReM_1)
	equivalence (ReImMino1(2), ImMino1, ImM_1)

	double precision ReImMino2(2), ReMino2, ImMino2
	double precision M_2, ReM_2, ImM_2
	equivalence (Mino2, M_2, ReImMino2)
	equivalence (ReImMino2(1), ReMino2, ReM_2)
	equivalence (ReImMino2(2), ImMino2, ImM_2)

	double precision ReImMino3(2), ReMino3, ImMino3
	double precision M_3, ReM_3, ImM_3
	equivalence (Mino3, M_3, ReImMino3)
	equivalence (ReImMino3(1), ReMino3, ReM_3)
	equivalence (ReImMino3(2), ImMino3, ImM_3)


* flavour-violating parameters

	double complex deltaSf(6,6,3:4)
	double complex UASf(6,6,3:4)
	double precision MASf(6,3:4), MASf2(6,3:4)

	common /fv_para/ UASf, MASf, MASf2, deltaSf

#ifndef UASfC
#define UASfC(i,j,t) DCONJG(UASf(i,j,t))

#define deltaSf_LL(i,j,t) deltaSf(i,j,t)
#define deltaSf_LR(i,j,t) deltaSf(i,j+3,t)
#define deltaSf_RL(i,j,t) deltaSf(j,i+3,t)
#define deltaSf_RR(i,j,t) deltaSf(i+3,j+3,t)

#define ImdeltaSf_LL(i,j,t) ReImdeltaSf(2,i,j,t)
#define ImdeltaSf_LR(i,j,t) ReImdeltaSf(2,i,j+3,t)
#define ImdeltaSf_RL(i,j,t) ReImdeltaSf(2,j,i+3,t)
#define ImdeltaSf_RR(i,j,t) ReImdeltaSf(2,i+3,j+3,t)
#endif

	double complex deltaSf_flat(6*6*2)
	equivalence (deltaSf, deltaSf_flat)

	double precision ReImdeltaSf(2,6,6,3:4)
	equivalence (deltaSf, ReImdeltaSf)

	double complex deltaLLuc, deltaLRuc
	double complex deltaRLuc, deltaRRuc
	double precision RedeltaLLuc, ImdeltaLLuc
	double precision RedeltaLRuc, ImdeltaLRuc
	double precision RedeltaRLuc, ImdeltaRLuc
	double precision RedeltaRRuc, ImdeltaRRuc
	equivalence (  deltaSf_LL(1,2,3), deltaLLuc, RedeltaLLuc)
	equivalence (ImdeltaSf_LL(1,2,3), ImdeltaLLuc)
	equivalence (  deltaSf_LR(1,2,3), deltaLRuc, RedeltaLRuc)
	equivalence (ImdeltaSf_LR(1,2,3), ImdeltaLRuc)
	equivalence (  deltaSf_RL(1,2,3), deltaRLuc, RedeltaRLuc)
	equivalence (ImdeltaSf_RL(1,2,3), ImdeltaRLuc)
	equivalence (  deltaSf_RR(1,2,3), deltaRRuc, RedeltaRRuc)
	equivalence (ImdeltaSf_RR(1,2,3), ImdeltaRRuc)

	double complex deltaLLct, deltaLRct
	double complex deltaRLct, deltaRRct
	double precision RedeltaLLct, ImdeltaLLct
	double precision RedeltaLRct, ImdeltaLRct
	double precision RedeltaRLct, ImdeltaRLct
	double precision RedeltaRRct, ImdeltaRRct
	equivalence (  deltaSf_LL(2,3,3), deltaLLct, RedeltaLLct)
	equivalence (ImdeltaSf_LL(2,3,3), ImdeltaLLct)
	equivalence (  deltaSf_LR(2,3,3), deltaLRct, RedeltaLRct)
	equivalence (ImdeltaSf_LR(2,3,3), ImdeltaLRct)
	equivalence (  deltaSf_RL(2,3,3), deltaRLct, RedeltaRLct)
	equivalence (ImdeltaSf_RL(2,3,3), ImdeltaRLct)
	equivalence (  deltaSf_RR(2,3,3), deltaRRct, RedeltaRRct)
	equivalence (ImdeltaSf_RR(2,3,3), ImdeltaRRct)

	double complex deltaLLut, deltaLRut
	double complex deltaRLut, deltaRRut
	double precision RedeltaLLut, ImdeltaLLut
	double precision RedeltaLRut, ImdeltaLRut
	double precision RedeltaRLut, ImdeltaRLut
	double precision RedeltaRRut, ImdeltaRRut
	equivalence (  deltaSf_LL(1,3,3), deltaLLut, RedeltaLLut)
	equivalence (ImdeltaSf_LL(1,3,3), ImdeltaLLut)
	equivalence (  deltaSf_LR(1,3,3), deltaLRut, RedeltaLRut)
	equivalence (ImdeltaSf_LR(1,3,3), ImdeltaLRut)
	equivalence (  deltaSf_RL(1,3,3), deltaRLut, RedeltaRLut)
	equivalence (ImdeltaSf_RL(1,3,3), ImdeltaRLut)
	equivalence (  deltaSf_RR(1,3,3), deltaRRut, RedeltaRRut)
	equivalence (ImdeltaSf_RR(1,3,3), ImdeltaRRut)

	double complex deltaLLds, deltaLRds
	double complex deltaRLds, deltaRRds
	double precision RedeltaLLds, ImdeltaLLds
	double precision RedeltaLRds, ImdeltaLRds
	double precision RedeltaRLds, ImdeltaRLds
	double precision RedeltaRRds, ImdeltaRRds
	equivalence (  deltaSf_LL(1,2,4), deltaLLds, RedeltaLLds)
	equivalence (ImdeltaSf_LL(1,2,4), ImdeltaLLds)
	equivalence (  deltaSf_LR(1,2,4), deltaLRds, RedeltaLRds)
	equivalence (ImdeltaSf_LR(1,2,4), ImdeltaLRds)
	equivalence (  deltaSf_RL(1,2,4), deltaRLds, RedeltaRLds)
	equivalence (ImdeltaSf_RL(1,2,4), ImdeltaRLds)
	equivalence (  deltaSf_RR(1,2,4), deltaRRds, RedeltaRRds)
	equivalence (ImdeltaSf_RR(1,2,4), ImdeltaRRds)

	double complex deltaLLsb, deltaLRsb
	double complex deltaRLsb, deltaRRsb
	double precision RedeltaLLsb, ImdeltaLLsb
	double precision RedeltaLRsb, ImdeltaLRsb
	double precision RedeltaRLsb, ImdeltaRLsb
	double precision RedeltaRRsb, ImdeltaRRsb
	equivalence (  deltaSf_LL(2,3,4), deltaLLsb, RedeltaLLsb)
	equivalence (ImdeltaSf_LL(2,3,4), ImdeltaLLsb)
	equivalence (  deltaSf_LR(2,3,4), deltaLRsb, RedeltaLRsb)
	equivalence (ImdeltaSf_LR(2,3,4), ImdeltaLRsb)
	equivalence (  deltaSf_RL(2,3,4), deltaRLsb, RedeltaRLsb)
	equivalence (ImdeltaSf_RL(2,3,4), ImdeltaRLsb)
	equivalence (  deltaSf_RR(2,3,4), deltaRRsb, RedeltaRRsb)
	equivalence (ImdeltaSf_RR(2,3,4), ImdeltaRRsb)

	double complex deltaLLdb, deltaLRdb
	double complex deltaRLdb, deltaRRdb
	double precision RedeltaLLdb, ImdeltaLLdb
	double precision RedeltaLRdb, ImdeltaLRdb
	double precision RedeltaRLdb, ImdeltaRLdb
	double precision RedeltaRRdb, ImdeltaRRdb
	equivalence (  deltaSf_LL(1,3,4), deltaLLdb, RedeltaLLdb)
	equivalence (ImdeltaSf_LL(1,3,4), ImdeltaLLdb)
	equivalence (  deltaSf_LR(1,3,4), deltaLRdb, RedeltaLRdb)
	equivalence (ImdeltaSf_LR(1,3,4), ImdeltaLRdb)
	equivalence (  deltaSf_RL(1,3,4), deltaRLdb, RedeltaRLdb)
	equivalence (ImdeltaSf_RL(1,3,4), ImdeltaRLdb)
	equivalence (  deltaSf_RR(1,3,4), deltaRRdb, RedeltaRRdb)
	equivalence (ImdeltaSf_RR(1,3,4), ImdeltaRRdb)


* THDM parameters

	double precision Lambda5
	double precision Yuk1, Yuk2, Yuk3
	logical thdm_digest

	common /thdm_para/ Lambda5
	common /thdm_para/ Yuk1, Yuk2, Yuk3
	common /thdm_para/ thdm_digest

