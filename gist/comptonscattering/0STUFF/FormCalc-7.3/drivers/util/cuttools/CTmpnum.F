* CTmpnum.F
* a dummy multiprecision numerator function
* this file is part of FormCalc
* last modified 2 Dec 10 th


	subroutine CTmpnum(q1in, res)
	implicit none
	double precision q1in(0:3)
	double complex res

	integer ini
	data ini /0/
	save ini

	if( ini .eq. 0 ) then
	  print *, "MP version not yet available"
	  ini = 1
	endif

	res = 0
	end

