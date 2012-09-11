* model_mssm.h
* declarations for model_mssm.F
* this file is part of FormCalc
* last modified 30 Nov 11 th


#include "model_sm.h"

	ComplexType UCha(2,2), VCha(2,2), ZNeu(4,4)
	ComplexType USf(2,2,4,3), UCSf(3,4,2:4,3), UUSf(3,4,2:4,3)
	ComplexType XHiggs(3,3,2)
	ComplexType Af(2:4,3,3), Xf(2:4,3), MUETB(2:4)
	ComplexType Atau, At, Ab, MUE
	ComplexType Mino1, Mino2, Mino3, SqrtEGl
	RealType MCha(2), MCha2(2), MNeu(4), MNeu2(4)
	RealType MSS(2,2:4,3), MSS2(2,2:4,3), DSf(2,4)
	RealType MSf(2,4,3), MSf2(2,4,3), MSusy, MGl, MGl2
	RealType MHiggs(4), MHiggs2(4), MHiggstree2(4)
	RealType CB, SB, TB, CB2, SB2, TB2, C2B, S2B
	RealType CA, SA, CA2, SA2, C2A, S2A
	RealType CAB, SAB, CBA, SBA, CBA2, SBA2
	RealType AlfasMT

	common /mssm_para/ UCha, VCha, ZNeu
	common /mssm_para/ USf, UCSf, UUSf
	common /mssm_para/ XHiggs
	common /mssm_para/ Af, Xf
	common /mssm_para/ Atau, At, Ab, MUE, MUETB
	common /mssm_para/ Mino1, Mino2, Mino3, SqrtEGl
	common /mssm_para/ MCha, MCha2, MNeu, MNeu2
	common /mssm_para/ MSS, MSS2, DSf
	common /mssm_para/ MSf, MSf2, MSusy, MGl, MGl2
	common /mssm_para/ MHiggs, MHiggs2, MHiggstree2
	common /mssm_para/ CB, SB, TB, CB2, SB2, TB2, C2B, S2B
	common /mssm_para/ CA, SA, CA2, SA2, C2A, S2A
	common /mssm_para/ CAB, SAB, CBA, SBA, CBA2, SBA2
	common /mssm_para/ AlfasMT

#ifndef USfC
#define USfC(i,j,t,g) Conjugate(USf(i,j,t,g))
#define UCSfC(i,j,t,g) Conjugate(UCSf(i,j,t,g))
#define UUSfC(i,j,t,g) Conjugate(UUSf(i,j,t,g))
#define USf2(i,j,t,g) Re(UCSf(i,j,t,g))
#define VChaC(i,j) Conjugate(VCha(i,j))
#define UChaC(i,j) Conjugate(UCha(i,j))
#define ZNeuC(i,j) Conjugate(ZNeu(i,j))
#define UHiggsC(i,j) Conjugate(UHiggs(i,j))
#define ZHiggsC(i,j) Conjugate(ZHiggs(i,j))
#define AfC(t,g1,g2) Conjugate(Af(t,g1,g2))
#define Mino3C Conjugate(Mino3)
#define MUEC Conjugate(MUE)
#define SqrtEGlC Conjugate(SqrtEGl)
#endif

	RealType Mh0, Mh02, MHH, MHH2, MA0, MA02, MHp, MHp2
	equivalence (MHiggs(1), Mh0), (MHiggs2(1), Mh02)
	equivalence (MHiggs(2), MHH), (MHiggs2(2), MHH2)
	equivalence (MHiggs(3), MA0), (MHiggs2(3), MA02)
	equivalence (MHiggs(4), MHp), (MHiggs2(4), MHp2)

	RealType Mh0tree2, MHHtree2, MA0tree2, MHptree2
	equivalence (MHiggstree2(1), Mh0tree2)
	equivalence (MHiggstree2(2), MHHtree2)
	equivalence (MHiggstree2(3), MA0tree2)
	equivalence (MHiggstree2(4), MHptree2)

	ComplexType UHiggs(3,3), ZHiggs(3,3)
	equivalence (XHiggs(1,1,1), UHiggs)
	equivalence (XHiggs(1,1,2), ZHiggs)

	RealType Af_flat(3*3*3), Xf_flat(3*3*3)
	equivalence (Af, Af_flat)
	equivalence (Xf, Xf_flat)

	RealType ReImAtau(2), ReAtau, ImAtau
	equivalence (Atau, ReImAtau)
	equivalence (ReImAtau(1), ReAtau), (ReImAtau(2), ImAtau)

	RealType ReImAt(2), ReAt, ImAt
	equivalence (At, ReImAt)
	equivalence (ReImAt(1), ReAt), (ReImAt(2), ImAt)

	RealType ReImAb(2), ReAb, ImAb
	equivalence (Ab, ReImAb)
	equivalence (ReImAb(1), ReAb), (ReImAb(2), ImAb)

	RealType ReImMUE(2), ReMUE, ImMUE
	equivalence (MUE, ReImMUE)
	equivalence (ReImMUE(1), ReMUE), (ReImMUE(2), ImMUE)

	RealType ReImMino1(2), ReMino1, ImMino1
	RealType M_1, ReM_1, ImM_1
	equivalence (Mino1, M_1, ReImMino1)
	equivalence (ReImMino1(1), ReMino1, ReM_1)
	equivalence (ReImMino1(2), ImMino1, ImM_1)

	RealType ReImMino2(2), ReMino2, ImMino2
	RealType M_2, ReM_2, ImM_2
	equivalence (Mino2, M_2, ReImMino2)
	equivalence (ReImMino2(1), ReMino2, ReM_2)
	equivalence (ReImMino2(2), ImMino2, ImM_2)

	RealType ReImMino3(2), ReMino3, ImMino3
	RealType M_3, ReM_3, ImM_3
	equivalence (Mino3, M_3, ReImMino3)
	equivalence (ReImMino3(1), ReMino3, ReM_3)
	equivalence (ReImMino3(2), ImMino3, ImM_3)


* flavour-violating parameters

	ComplexType deltaSf(6,6,3:4)
	ComplexType UASf(6,6,3:4)
	RealType MASf(6,3:4), MASf2(6,3:4)

	common /fv_para/ UASf, MASf, MASf2, deltaSf

#ifndef UASfC
#define UASfC(i,j,t) Conjugate(UASf(i,j,t))

#define deltaSf_LL(i,j,t) deltaSf(i,j,t)
#define deltaSf_LR(i,j,t) deltaSf(i,j+3,t)
#define deltaSf_RL(i,j,t) deltaSf(j,i+3,t)
#define deltaSf_RR(i,j,t) deltaSf(i+3,j+3,t)

#define ImdeltaSf_LL(i,j,t) ReImdeltaSf(2,i,j,t)
#define ImdeltaSf_LR(i,j,t) ReImdeltaSf(2,i,j+3,t)
#define ImdeltaSf_RL(i,j,t) ReImdeltaSf(2,j,i+3,t)
#define ImdeltaSf_RR(i,j,t) ReImdeltaSf(2,i+3,j+3,t)
#endif

	ComplexType deltaSf_flat(6*6*2)
	equivalence (deltaSf, deltaSf_flat)

	RealType ReImdeltaSf(2,6,6,3:4)
	equivalence (deltaSf, ReImdeltaSf)

	ComplexType deltaLRuc
	RealType RedeltaLRuc, ImdeltaLRuc
	equivalence (deltaSf_LR(1,2,3), deltaLRuc, RedeltaLRuc)
	equivalence (ImdeltaSf_LR(1,2,3), ImdeltaLRuc)
	ComplexType deltaLRct
	RealType RedeltaLRct, ImdeltaLRct
	equivalence (deltaSf_LR(2,3,3), deltaLRct, RedeltaLRct)
	equivalence (ImdeltaSf_LR(2,3,3), ImdeltaLRct)
	ComplexType deltaLRut
	RealType RedeltaLRut, ImdeltaLRut
	equivalence (deltaSf_LR(1,3,3), deltaLRut, RedeltaLRut)
	equivalence (ImdeltaSf_LR(1,3,3), ImdeltaLRut)

	ComplexType deltaRLuc
	RealType RedeltaRLuc, ImdeltaRLuc
	equivalence (deltaSf_RL(1,2,3), deltaRLuc, RedeltaRLuc)
	equivalence (ImdeltaSf_RL(1,2,3), ImdeltaRLuc)
	ComplexType deltaRLct
	RealType RedeltaRLct, ImdeltaRLct
	equivalence (deltaSf_RL(2,3,3), deltaRLct, RedeltaRLct)
	equivalence (ImdeltaSf_RL(2,3,3), ImdeltaRLct)
	ComplexType deltaRLut
	RealType RedeltaRLut, ImdeltaRLut
	equivalence (deltaSf_RL(1,3,3), deltaRLut, RedeltaRLut)
	equivalence (ImdeltaSf_RL(1,3,3), ImdeltaRLut)

	ComplexType deltaRRuc
	RealType RedeltaRRuc, ImdeltaRRuc
	equivalence (deltaSf_RR(1,2,3), deltaRRuc, RedeltaRRuc)
	equivalence (ImdeltaSf_RR(1,2,3), ImdeltaRRuc)
	ComplexType deltaRRct
	RealType RedeltaRRct, ImdeltaRRct
	equivalence (deltaSf_RR(2,3,3), deltaRRct, RedeltaRRct)
	equivalence (ImdeltaSf_RR(2,3,3), ImdeltaRRct)
	ComplexType deltaRRut
	RealType RedeltaRRut, ImdeltaRRut
	equivalence (deltaSf_RR(1,3,3), deltaRRut, RedeltaRRut)
	equivalence (ImdeltaSf_RR(1,3,3), ImdeltaRRut)

	ComplexType deltaLL12
	RealType RedeltaLL12, ImdeltaLL12
	equivalence (deltaSf_LL(1,2,4), deltaLL12, RedeltaLL12)
	equivalence (ImdeltaSf_LL(1,2,4), ImdeltaLL12)
	ComplexType deltaLL23
	RealType RedeltaLL23, ImdeltaLL23
	equivalence (deltaSf_LL(2,3,4), deltaLL23, RedeltaLL23)
	equivalence (ImdeltaSf_LL(2,3,4), ImdeltaLL23)
	ComplexType deltaLL13
	RealType RedeltaLL13, ImdeltaLL13
	equivalence (deltaSf_LL(1,3,4), deltaLL13, RedeltaLL13)
	equivalence (ImdeltaSf_LL(1,3,4), ImdeltaLL13)

	ComplexType deltaLRds
	RealType RedeltaLRds, ImdeltaLRds
	equivalence (deltaSf_LR(1,2,4), deltaLRds, RedeltaLRds)
	equivalence (ImdeltaSf_LR(1,2,4), ImdeltaLRds)
	ComplexType deltaLRsb
	RealType RedeltaLRsb, ImdeltaLRsb
	equivalence (deltaSf_LR(2,3,4), deltaLRsb, RedeltaLRsb)
	equivalence (ImdeltaSf_LR(2,3,4), ImdeltaLRsb)
	ComplexType deltaLRdb
	RealType RedeltaLRdb, ImdeltaLRdb
	equivalence (deltaSf_LR(1,3,4), deltaLRdb, RedeltaLRdb)
	equivalence (ImdeltaSf_LR(1,3,4), ImdeltaLRdb)

	ComplexType deltaRLds
	RealType RedeltaRLds, ImdeltaRLds
	equivalence (deltaSf_RL(1,2,4), deltaRLds, RedeltaRLds)
	equivalence (ImdeltaSf_RL(1,2,4), ImdeltaRLds)
	ComplexType deltaRLsb
	RealType RedeltaRLsb, ImdeltaRLsb
	equivalence (deltaSf_RL(2,3,4), deltaRLsb, RedeltaRLsb)
	equivalence (ImdeltaSf_RL(2,3,4), ImdeltaRLsb)
	ComplexType deltaRLdb
	RealType RedeltaRLdb, ImdeltaRLdb
	equivalence (deltaSf_RL(1,3,4), deltaRLdb, RedeltaRLdb)
	equivalence (ImdeltaSf_RL(1,3,4), ImdeltaRLdb)

	ComplexType deltaRRds
	RealType RedeltaRRds, ImdeltaRRds
	equivalence (deltaSf_RR(1,2,4), deltaRRds, RedeltaRRds)
	equivalence (ImdeltaSf_RR(1,2,4), ImdeltaRRds)
	ComplexType deltaRRsb
	RealType RedeltaRRsb, ImdeltaRRsb
	equivalence (deltaSf_RR(2,3,4), deltaRRsb, RedeltaRRsb)
	equivalence (ImdeltaSf_RR(2,3,4), ImdeltaRRsb)
	ComplexType deltaRRdb
	RealType RedeltaRRdb, ImdeltaRRdb
	equivalence (deltaSf_RR(1,3,4), deltaRRdb, RedeltaRRdb)
	equivalence (ImdeltaSf_RR(1,3,4), ImdeltaRRdb)
