* diag.h
* global declarations for the Diag routines
* this file is part of Diag
* last modified 23 Nov 07 th


* The maximum dimension of a matrix, needed for allocating internal
* memory, i.e. the routines handle at most MAXDIM-by-MAXDIM matrices.

#define MAXDIM 16


* A matrix is considered diagonal if the sum of the squares
* of the off-diagonal elements is less than EPS.  SYM_EPS is
* half of EPS since only the upper triangle is counted for
* symmetric matrices.
* 52 bits is the mantissa length for IEEE double precision.

#define EPS 2D0**(-102)

#define SYM_EPS 2D0**(-103)

#define DBL_EPS 2D0**(-52)


* The transposed versions are needed for C, which has row-major
* matrix access.

#ifdef TRANSPOSE

#define Element(A,i,j) A(j,i)
#define HEigensystem HEigensystemT
#define SEigensystem SEigensystemT
#define CEigensystem CEigensystemT
#define TakagiFactor TakagiFactorT
#define SVD SVDT

#else

#define Element(A,i,j) A(i,j)

#endif

