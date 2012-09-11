* util.h
* prototypes for the functions in util.a
* this file is part of FormCalc
* last modified 21 Oct 09 th


	double precision ThreeMom
	double precision SInvariant, TInvariant
	double complex Pair, Eps
	double complex SxS, SeS
	integer VxS, VeS, BxS, BeS

	external ThreeMom
	external SInvariant, TInvariant
	external Pair, Eps
	external SxS, SeS
	external VxS, VeS, BxS, BeS

#ifndef LEGS
#define LEGS 1
#endif

	double complex vec(2,2,8,0:LEGS)
	common /vectors/ vec

	double complex kcomp(32,0:LEGS)
	equivalence (vec, kcomp)

	double precision momspec(16,LEGS)
	common /momenta/ momspec

	integer*8 two
	parameter (two = 2)


#ifndef SPEC_M

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

#define k(i) (8*i+1)
#define s(i) (8*i+3)
#define e(i) (8*i+3+Hel(i))
#define ec(i) (8*i+3-Hel(i))
#define Spinor(i,s,om) (s*2*Hel(i)+16*i+om+5)
#define DottedSpinor(i,s,om) (s*2*Hel(i)+16*i+om+7)

#define K(f,i) iand(f,255)*two**(8*(i-1))

#define signbit(i) ibits(i,31,1)
#define IndexDelta(i,j) signbit(ieor(i,j)-1)
#define Digit(i) char(i+48)
#define Polar(r,theta) r*exp(cI*degree*theta)
#define Sq(c) DBLE((c)*DCONJG(c))

#define Error(msg) call m_(fail, __LINE__, __FILE__, msg)
#define Warning(msg) call m_(0, 0, __FILE__, msg)
#define INFO print *,
#define DEB print *,

#define Cut(c,m) (m)*(c)

#define CUT_MIN 1
#define CUT_MAX 2

#define CUT_COSTH 4
#define CUT_COSTHCMS 16
#define CUT_COSTH_E 64
#define CUT_COSTH_K 65
#define CUT_MREM 256
#define CUT_MREM_E 1024
#define CUT_MREM_K 1025
#define CUT_MREM_ET 4096
#define CUT_MREM_KT 4097
#define CUT_MREM_RAP 16384
#define CUT_MREM_PRAP 16385

#endif

