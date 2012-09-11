* main.F
* this is the main program for the generated Fortran code,
* it scans the command line and then invokes ParameterScan
* this file is part of FormCalc
* last modified 30 Nov 11 th


#include "xsection.F"

#ifndef MMA

	program main
	implicit none

#if U77EXT
	integer iargc, lnblnk
	external iargc, lnblnk
#endif

	RealType sqrtSfrom, sqrtSto, sqrtSstep
	integer serialfrom, serialto, serialstep
	integer argc, next, fail
	character*LEGS pol
	character*128 argv(5), tmp
	character*512 dir

#ifdef DEBUG
c	call setfpu
#endif

	call getarg(0, dir)
	next = len_trim(dir) + 1

	argc = 0
1	call getarg(argc + 1, dir(next:))
	if( dir(next:next) .eq. ':' ) then
	  argc = argc + 1
	  argv(argc) = dir(next+1:)
	  next = next + len_trim(dir(next:))
	  goto 1
	endif

	pol = dir(next:)

	sqrtSfrom = -1
	sqrtSto = 0
	sqrtSstep = SQRTSSTEP
	call getarg(argc + 2, tmp)
	read(tmp, *, iostat=fail, end=991, err=999)
     &    sqrtSfrom, sqrtSto, sqrtSstep
991	if( sqrtSfrom .eq. -1 ) goto 999

	call ProcessIni(fail, pol,
     &    sqrtSfrom, sqrtSto, sqrtSstep,
     &    argc, argv)
	if( fail .ne. 0 ) goto 999

	if( sqrtSto .eq. 0 ) then
100	  format(".", A, ".", I5.5)
	  write(dir(next:), 100) pol, int(sqrtSfrom)
	else
101	  format(".", A, ".", I5.5, ",", I5.5, ",", I5.5)
	  write(dir(next:), 101) pol,
     &      int(sqrtSfrom), int(sqrtSto), int(sqrtSstep)
	endif

	serialfrom = 0
	serialto = 2**30
	serialstep = 1
	if( argc + 3 .le. iargc() ) then
	  call getarg(argc + 3, tmp)
	  serialfrom = -1
	  read(tmp, *, iostat=fail, end=992, err=999)
     &      serialfrom, serialto, serialstep
992	  if( serialfrom .eq. -1 ) goto 999
	endif

	call ParameterScan(dir, serialfrom, serialto, serialstep)

	call ProcessExi
	goto 2

99	format (/"USAGE:  ", A, " [:args] p1...p", I1,
     &    " sqrtSfrom[,to[,step]] [serialfrom[,to[,step]]]"
     &    //"p1..", I1, " = +  for right-circular polarization"
     &    /"        -  for left-circular       ''"
     &    /"        l  for longitudinal        ''"
     &    /"        t  for transverse          ''"
     &    /"        u  for unpolarized particle"
     &    //"sqrtSfrom = compute the differential cross-section"
     &    /"sqrtSfrom,to,step = compute the total cross-section"
     &    //"serialfrom,to,step = range of serial numbers"/)

999	call getarg(0, tmp)
	print 99, tmp(1:len_trim(tmp)), LEGS, LEGS
	stop

2	continue
	end

#endif

