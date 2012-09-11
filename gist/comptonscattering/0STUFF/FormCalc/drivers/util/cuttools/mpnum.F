* mpnum.F
* a dummy multiprecision numerator function
* this file is part of FormCalc
* last modified 21 Oct 09 th


	subroutine mpnum(q1in, qt2, res)
	implicit none
	double precision q1in(0:3), qt2
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

