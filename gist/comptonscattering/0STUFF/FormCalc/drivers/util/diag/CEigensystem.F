* CEigensystem.F
* diagonalization of a complex n-by-n matrix using the Jacobi algorithm
* code adapted from the "Handbook" routines for complex A
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of the Diag library
* last modified 25 Nov 09 th

#include "diag.h"


************************************************************************
** CEigensystem diagonalizes a general complex n-by-n matrix.
** Input: n, A = n-by-n matrix
** Output: d = vector of eigenvalues, U = transformation matrix
** these fulfill diag(d) = U A U^-1.

	subroutine CEigensystem(n, A,ldA, d, U,ldU, sort)
	implicit none
	integer n, ldA, ldU, sort
	double complex A(ldA,*), U(ldU,*), d(*)

	integer p, q, j
	double precision red, off, thresh
	double complex delta, t, invc, sx, sy, tx, ty
	double complex x, y
	double complex ev(2,MAXDIM)

	integer sweep
	common /nsweeps/ sweep

	double precision sq
	double complex c
	sq(c) = DBLE(c*DCONJG(c))

	if( n .gt. MAXDIM ) then
	  print *, "Dimension too large"
	  d(1) = -999
	  return
	endif

	do p = 1, n
	  ev(1,p) = 0
	  ev(2,p) = A(p,p)
	  d(p) = ev(2,p)
	enddo

	do p = 1, n
	  do q = 1, n
	    U(q,p) = 0
	  enddo
	  U(p,p) = 1
	enddo

	red = .01D0/n**4

	do sweep = 1, 50
	  off = 0
	  do q = 2, n
	    do p = 1, q - 1
	      off = off + sq(A(p,q)) + sq(A(q,p))
	    enddo
	  enddo
	  if( .not. off .gt. EPS ) goto 1

	  thresh = 0
	  if( sweep .lt. 4 ) thresh = off*red

	  do q = 2, n
	    do p = 1, q - 1
	      off = sq(A(p,q)) + sq(A(q,p))
	      if( sweep .gt. 4 .and. off .lt.
     &              EPS*max(sq(ev(2,p)), sq(ev(2,q))) ) then
	        A(p,q) = 0
	        A(q,p) = 0
	      else
	        if( off .gt. thresh ) then
	          delta = A(p,q)*A(q,p)
	          x = .5D0*(ev(2,p) - ev(2,q))
	          y = sqrt(x**2 + delta)
	          t = x - y
	          x = x + y
	          if( abs(t) .lt. abs(x) ) t = x

	          t = 1/t
	          delta = delta*t
	          ev(1,p) = ev(1,p) + delta
	          ev(2,p) = d(p) + ev(1,p)
	          ev(1,q) = ev(1,q) - delta
	          ev(2,q) = d(q) + ev(1,q)

	          invc = sqrt(delta*t + 1)
	          x = t/invc
	          t = t/(invc + 1)
	          sx = x*Element(A,p,q)
	          ty = t*Element(A,p,q)
	          sy = x*Element(A,q,p)
	          tx = t*Element(A,q,p)

	          do j = 1, n
	            x = Element(A,j,p)
	            y = Element(A,j,q)
	            Element(A,j,p) = x + sy*(y - ty*x)
	            Element(A,j,q) = y - sx*(x + tx*y)
	            x = Element(A,p,j)
	            y = Element(A,q,j)
	            Element(A,p,j) = x + sx*(y - tx*x)
	            Element(A,q,j) = y - sy*(x + ty*y)
	          enddo

	          A(p,q) = 0
	          A(q,p) = 0

	          do j = 1, n
	            x = Element(U,p,j)
	            y = Element(U,q,j)
	            Element(U,p,j) = x + sx*(y - tx*x)
	            Element(U,q,j) = y - sy*(x + ty*y)
	          enddo
	        endif
	      endif
	    enddo
	  enddo

	  do p = 1, n
	    ev(1,p) = 0
	    d(p) = ev(2,p)
	  enddo
	enddo

	print *, "Bad convergence in CEigensystem"

1	if( sort .eq. 0 ) return

* sort the eigenvalues by their real part

	do p = 1, n - 1
	  j = p
	  t = d(p)
	  do q = p + 1, n
	    if( sort*(DBLE(t) - DBLE(d(q))) .gt. 0 ) then
	      j = q
	      t = d(q)
	    endif
	  enddo

	  if( j .ne. p ) then
	    d(j) = d(p)
	    d(p) = t
	    do q = 1, n
	      x = Element(U,p,q)
	      Element(U,p,q) = Element(U,j,q)
	      Element(U,j,q) = x
	    enddo
	  endif
	enddo
	end

