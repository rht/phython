* opp.h
* declarations for the OPP routines
* this file is part of FormCalc
* last modified 30 Nov 11 th


#ifndef Acut
#ifdef SAMURAI
#define Acut SAAcut
#define Bcut SABcut
#define Ccut SACcut
#define Dcut SADcut
#define Ecut SAEcut
#define Fcut SAFcut
#elif defined CUTTOOLS
#define Acut CTAcut
#define Bcut CTBcut
#define Ccut CTCcut
#define Dcut CTDcut
#define Ecut CTEcut
#define Fcut CTFcut
#else
#error Neither SAMURAI nor CUTTOOLS defined
#endif
#endif

	ComplexType Acut, Bcut, Ccut, Dcut, Ecut, Fcut
	external Acut, Bcut, Ccut, Dcut, Ecut, Fcut

