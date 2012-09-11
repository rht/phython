#ifndef FTYPES_H
#define FTYPES_H

#if NOUNDERSCORE
#define FORTRAN(s) s
#else
#define FORTRAN(s) s##_
#endif

#if QUAD

#pragma pack(push, 1)
typedef union {
  long double r10;
  struct {
    unsigned long long frac;
    unsigned short exp;
  } i10;
  struct {
    char zero[6];
    unsigned long long frac;
    unsigned short exp;
  } i16;
  unsigned long long i8[2];
} real16;
#pragma pack(pop)

static inline real16 ToREAL(long double r) {
  real16 new;
  new.i8[0] = 0;
  new.i16.frac = ((real16 *)&r)->i10.frac << 1;
  new.i16.exp = ((real16 *)&r)->i10.exp;
  return new;
}

static inline long double ToReal(real16 r) {
  real16 new;
  const long long z = r.i16.frac | (r.i16.exp & 0x7fff);
  new.i10.frac = (r.i16.frac >> 1) | ((z | -z) & 0x8000000000000000);
  new.i10.exp = r.i16.exp;
  return new.r10;
}

#define Real long double
typedef real16 REAL;

#else

#define ToReal(r) (r)
#define ToREAL(r) (r)
#define Real double
typedef double REAL;

#endif

typedef int INTEGER;
typedef const INTEGER CINTEGER;
typedef const REAL CREAL;
typedef struct { REAL re, im; } COMPLEX;
typedef const COMPLEX CCOMPLEX;
typedef char CHARACTER;
typedef const CHARACTER CCHARACTER;

#ifdef __cplusplus

#include <complex>
typedef std::complex<Real> Complex;
#define ToComplex(c) Complex(ToReal((c).re), ToReal((c).im))
#define ToComplex2(r,i) Complex(r, i)
#define Re(x) std::real(x)
#define Im(x) std::imag(x)

#else

#include <complex.h>
typedef Real complex Complex;
#define ToComplex(c) (ToReal((c).re) + I*ToReal((c).im))
#define ToComplex2(r,i) (r + I*(i))
#define Re(x) creal(x)
#define Im(x) cimag(x)

#endif

typedef const Real cReal;
typedef const Complex cComplex;

#endif

