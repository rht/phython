* lt.F
* the LoopTools command-line interface to Bget, Cget, Dget, Eget
* this file is part of LoopTools
* last modified 21 Jan 11 th

#include "defs.h"


	program LoopTools
	implicit none

#include "lt.h"

#if U77EXT
	integer iargc
	external iargc
#endif

	integer argc
	character argv*100
	double precision para(Pee)
	integer i, npoint, fail
	memindex b

	double complex A0, A00
	memindex cachelookup
	external cachelookup, A0, A00
	external Bcoeff, Ccoeff, Dcoeff, Ecoeff

	integer npara(1:5), key(1:5)
	data npara /1, Pbb, Pcc, Pdd, Pee/
	data key /KeyA0, KeyBget, KeyC0, KeyD0, KeyEget/

	argc = iargc() 
	do npoint = 1, 5
	  if( argc .eq. npara(npoint) .or.
     &        argc .eq. npara(npoint) + 1 ) goto 1
	enddo

999	print *, "Usage:  lt `parameters' [versionkey]"
	print *, "computes the n-point one-loop integrals"
	print *, "n depends on `parameters':"
	print *, "  n = 1: m"
	print *, "  n = 2: p m1 m2"
	print *, "  n = 3: p1 p2 p1p2 m1 m2 m3"
	print *, "  n = 4: p1 p2 p3 p4 p1p2 p2p3 m1 m2 m3 m4"
	print *, "  n = 5: p1 p2 p3 p4 p5 p1p2 p2p3 p3p4 p4p5 p5p1"//
     &                   " m1 m2 m3 m4 m5"
	print *, "versionkey can be one of:"
	print *, "  0 = compute version a (same as no versionkey)"
	print *, "  1 = compute version b"
	print *, "  2 = compute a and b, compare, return a"
	print *, "  3 = compute a and b, compare, return b"
	call exit(1)

1	do i = 1, npara(npoint)
	  call getarg(i, argv)
	  read(argv, *, iostat=fail, err=999) para(i)
	enddo

	call ltini

	if( argc .eq. i ) then
	  call getarg(i, argv)
	  read(argv, *, iostat=fail, err=999) i
	  if( i .lt. 0 .or. i .gt. 3 ) goto 999
	  print *, "using versionkey =", i
	  versionkey = ishft(i, key(npoint))
	endif

	if( npoint .eq. 1 ) then
	  print *, "m       =", para(1)
	  print *, "A0      =", A0(para(1))
	  print *, "A00     =", A00(para(1))
	else
	  call DumpPara(npoint, para, 1, " ")
	  if( npoint .eq. 2 ) then
	    b = cachelookup(para, Bval(1,0), Bcoeff, Pbb, Nbb)
	  else if( npoint .eq. 3 ) then
	    b = cachelookup(para, Cval(1,0), Ccoeff, Pcc, Ncc)
	  else if( npoint .eq. 4 ) then
	    b = cachelookup(para, Dval(1,0), Dcoeff, Pdd, Ndd)
	  else
	    b = cachelookup(para, Eval(1,0), Ecoeff, Pee, Nee)
	  endif
	  call DumpCoeff(npoint, Nval(npoint,1,b))
	endif

	call ltexi
	end

