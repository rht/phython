* types.h
* real-based type declarations
* this file is part of FormCalc
* last modified 30 Nov 11 th


#ifndef TYPES_H
#define TYPES_H

#ifdef QUAD
#define RealType real*16
#define ComplexType complex*32
#define Re QEXT
#define Im QIMAG
#define Conjugate QCONJG
#define ToComplex QCMPLX
#else
#define RealType double precision
#define ComplexType double complex
#define Re DBLE
#define Im DIMAG
#define Conjugate DCONJG
#define ToComplex DCMPLX
#endif

#define Sq(c) Re((c)*Conjugate(c))

#endif

