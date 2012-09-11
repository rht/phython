* spectrum.f
* compute the CompAZ spectrum
* use e.g. in gnuplot as plot "< ./spectrum 250" with lines
* last modified 2 Feb 05 th


	program spectrum
	implicit none

	double precision energy, fraction
	character*100 argv

	double precision CompAZ
	integer iargc
	external iargc, CompAZ

	if( iargc() .ne. 1 ) then
	  call getarg(0, argv)
	  print *, "Usage: spectrum E"
	  print *, "Compute the CompAZ spectrum for energy E"
	  return
	endif

	call getarg(1, argv)
	read(argv, *) energy

	do fraction = 0, 1, .01
	  print *, fraction, CompAZ(fraction, energy, 0)
	enddo
	end

