* TakagiFactor.F
* computes the Takagi factorization of a complex symmetric matrix
* code adapted from the "Handbook" routines
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of the Diag library
* last modified 17 Feb 10 th

#include "diag.h"


************************************************************************
** TakagiFactor factorizes a complex symmetric n-by-n matrix
** Input: n, A = n-by-n matrix, complex symmetric
** (only the upper triangle of A needs to be filled).
** Output: d = vector of diagonal values, U = transformation matrix
** these fulfill diag(d) = U^* A U^+ with U unitary.

	subroutine TakagiFactor(n, A,ldA, d, U,ldU, sort)
	implicit none
	integer n, ldA, ldU, sort
	double complex A(ldA,*), U(ldU,*)
	double precision d(*)

	integer p, q, j
	double precision red, off, thresh
	double precision sqp, sqq, t, invc
	double complex f, x, y
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
	enddo

	do p = 1, n
	  do q = 1, n
	    U(q,p) = 0
	  enddo
	  U(p,p) = 1
	enddo

	red = .04D0/n**4

	do sweep = 1, 50
	  off = 0
	  do q = 2, n
	    do p = 1, q - 1
	      off = off + sq(Element(A,p,q))
	    enddo
	  enddo
	  if( .not. off .gt. SYM_EPS ) goto 1

	  thresh = 0
	  if( sweep .lt. 4 ) thresh = off*red

	  do q = 2, n
	    do p = 1, q - 1
	      off = sq(Element(A,p,q))
	      sqp = sq(ev(2,p))
	      sqq = sq(ev(2,q))
	      if( sweep .gt. 4 .and.
     &            off .lt. SYM_EPS*max(sqp, sqq) ) then
	        Element(A,p,q) = 0
	      else
	        if( off .gt. thresh ) then
	          t = .5D0*abs(sqp - sqq)
	          if( t .eq. 0 ) then
	            f = 1
	            if( sqp .ne. 0 ) f = sqrt(ev(2,q)/ev(2,p))
	          else
	            f = sign(1D0, sqp - sqq)*
     &                (ev(2,q)*DCONJG(Element(A,p,q)) +
     &                  DCONJG(ev(2,p))*Element(A,p,q))
	          endif
	          t = t + sqrt(t**2 + DBLE(f)**2 + DIMAG(f)**2)
	          f = f/t

	          ev(1,p) = ev(1,p) + Element(A,p,q)*DCONJG(f)
	          ev(2,p) = A(p,p) + ev(1,p)
	          ev(1,q) = ev(1,q) - Element(A,p,q)*f
	          ev(2,q) = A(q,q) + ev(1,q)

	          t = DBLE(f)**2 + DIMAG(f)**2
	          invc = sqrt(t + 1)
	          f = f/invc
	          t = t/(invc*(invc + 1))

	          do j = 1, p - 1
	            x = Element(A,j,p)
	            y = Element(A,j,q)
	            Element(A,j,p) = x + (DCONJG(f)*y - t*x)
	            Element(A,j,q) = y - (f*x + t*y)
	          enddo

	          do j = p + 1, q - 1
	            x = Element(A,p,j)
	            y = Element(A,j,q)
	            Element(A,p,j) = x + (DCONJG(f)*y - t*x)
	            Element(A,j,q) = y - (f*x + t*y)
	          enddo

	          do j = q + 1, n
	            x = Element(A,p,j)
	            y = Element(A,q,j)
	            Element(A,p,j) = x + (DCONJG(f)*y - t*x)
	            Element(A,q,j) = y - (f*x + t*y)
	          enddo

	          Element(A,p,q) = 0

	          do j = 1, n
	            x = Element(U,p,j)
	            y = Element(U,q,j)
	            Element(U,p,j) = x + (f*y - t*x)
	            Element(U,q,j) = y - (DCONJG(f)*x + t*y)
	          enddo
	        endif
	      endif
	    enddo
	  enddo

	  do p = 1, n
	    ev(1,p) = 0
	    A(p,p) = ev(2,p)
	  enddo
	enddo

	print *, "Bad convergence in TakagiFactor"

1	continue

* make the diagonal elements nonnegative

	do p = 1, n
	  d(p) = abs(A(p,p))
	  if( d(p) .gt. DBL_EPS .and. d(p) .ne. DBLE(A(p,p)) ) then
	    f = sqrt(A(p,p)/d(p))
	    do q = 1, n
	      Element(U,p,q) = Element(U,p,q)*f
	    enddo
	  endif
	enddo

	if( sort .eq. 0 ) return

* sort the eigenvalues

	do p = 1, n - 1
	  j = p
	  t = d(p)
	  do q = p + 1, n
	    if( sort*(t - d(q)) .gt. 0 ) then
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

