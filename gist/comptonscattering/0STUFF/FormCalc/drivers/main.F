* main.F
* this is the main program for the generated Fortran code,
* it scans the command line and then invokes ParameterScan
* this file is part of FormCalc
* last modified 12 Jun 09 th


#include "xsection.F"

#ifndef MMA

	program main
	implicit none

#if U77EXT
	integer iargc, lnblnk
	external iargc, lnblnk
#endif

	double precision sqrtSfrom, sqrtSto, sqrtSstep
	integer serialfrom, serialto, serialstep
	integer argc, fail, next
	character*200 argv
	character*512 dir

#ifdef DEBUG
	call setfpu
#endif

	argc = iargc()
	if( argc .lt. 2 .or. argc .gt. 3 ) goto 999

	sqrtSfrom = -1
	sqrtSto = 0
	sqrtSstep = SQRTSSTEP
	call getarg(2, argv)
	read(argv, *, iostat=fail, end=991, err=999)
     &    sqrtSfrom, sqrtSto, sqrtSstep
991	if( sqrtSfrom .eq. -1 ) goto 999

	call getarg(1, argv)

	call ProcessIni(fail, argv, sqrtSfrom, sqrtSto, sqrtSstep)
	if( fail .ne. 0 ) goto 999

	call getarg(0, dir)
	next = lnblnk(dir) + 1
	if( sqrtSto .eq. 0 ) then
100	  format(".", A, ".", I5.5)
	  write(dir(next:), 100) argv(1:LEGS), int(sqrtSfrom)
	else
101	  format(".", A, ".", I5.5, ",", I5.5, ",", I5.5)
	  write(dir(next:), 101) argv(1:LEGS),
     &      int(sqrtSfrom), int(sqrtSto), int(sqrtSstep)
	endif

	serialfrom = 0
	serialto = 2**30
	serialstep = 1
	if( argc .gt. 2 ) then
	  call getarg(3, argv)
	  serialfrom = -1
	  read(argv, *, iostat=fail, end=992, err=999)
     &      serialfrom, serialto, serialstep
992	  if( serialfrom .eq. -1 ) goto 999
	endif

	call ParameterScan(dir, serialfrom, serialto, serialstep)

	call ffexi
	goto 1

99	format (/"USAGE:  ", A, " p1...p", I1,
     &    " sqrtSfrom[,to[,step]] [serialfrom[,to[,step]]]"
     &    //"p1..", I1, " = +  for right-circular polarization"
     &    /"        -  for left-circular       ''"
     &    /"        l  for longitudinal        ''"
     &    /"        t  for transverse          ''"
     &    /"        u  for unpolarized particle"
     &    //"sqrtSfrom = compute the differential cross-section"
     &    /"sqrtSfrom,to,step = compute the total cross-section"
     &    //"serialfrom,to,step = range of serial numbers"/)

999	call getarg(0, argv)
	print 99, argv(1:lnblnk(argv)), LEGS, LEGS
	stop

1	continue
	end

#endif

