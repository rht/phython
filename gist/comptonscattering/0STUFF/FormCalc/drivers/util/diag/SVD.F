* SVD.F
* singular value decomposition of an m-by-n matrix
* this file is part of the Diag library
* last modified 25 Nov 09 th

#include "diag.h"


************************************************************************
** SVD performs a singular value decomposition.
** Input: m, n, A = m-by-n matrix.
** Output: d = nm-vector of singular values,
** V = nm-by-m left transformation matrix,
** W = nm-by-n right transformation matrix, nm = min(m, n),
** these fulfill diag(d) = V^* A W^+.

	subroutine SVD(m, n, Ao,ldA, d, Vo,ldV, Wo,ldW, sort)
	implicit none
	integer m, n, ldA, ldV, ldW, sort
	double complex Ao(ldA,*), Vo(ldV,*), Wo(ldW,*)
	double precision d(*)

	integer nx, nm, p, q, j, rev, pi(MAXDIM)
	double precision red, off, thresh
	double precision t, dv, dw, xv, xw, invc
	double complex x, y, sv, sw, tv, tw, evp, evq, f
	double complex A(MAXDIM,MAXDIM)
	double complex VW(MAXDIM,MAXDIM,0:1)

* note: for better cache efficiency, the Vx, Wx, VWx arrays
* contain the *transpose* of the transformation matrices
	double complex V(MAXDIM,MAXDIM)
	double complex W(MAXDIM,MAXDIM)
	equivalence (VW(1,1,0), V)
	equivalence (VW(1,1,1), W)

	integer sweep
	common /nsweeps/ sweep

	double precision sq
	double complex c
	sq(c) = DBLE(c*DCONJG(c))

	nx = max(m, n)

	if( nx .gt. MAXDIM ) then
	  print *, "Dimension too large"
	  d(1) = -999
	  return
	endif

	do p = 1, nx
	  do q = 1, nx
	    A(q,p) = 0
	    V(q,p) = 0
	    W(q,p) = 0
	  enddo
	  V(p,p) = 1
	  W(p,p) = 1
	enddo

	if( m .ge. n ) then
	  do p = 1, n
	    do q = 1, m
	      A(q,p) = Element(Ao,q,p)
	    enddo
	  enddo
	else
	  do p = 1, n
	    do q = 1, m
	      A(p,q) = Element(Ao,q,p)
	    enddo
	  enddo
	endif

	red = .01D0/nx**4

	do sweep = 1, 50
	  off = 0
	  do q = 1, nx
	    do p = 1, q - 1
	      off = off + sq(A(p,q)) + sq(A(q,p))
	    enddo
	  enddo
	  if( .not. off .gt. EPS ) goto 1

	  thresh = 0
	  if( sweep .lt. 4 ) thresh = off*red

	  do q = 2, nx
	    do p = 1, q - 1
	      off = sq(A(p,q)) + sq(A(q,p))
	      if( sweep .gt. 4 .and. off .lt.
     &              EPS*max(sq(A(p,p)), sq(A(q,q))) ) then
	        A(p,q) = 0
	        A(q,p) = 0
	      else
	        if( off .gt. thresh ) then
	          xv = DBLE((A(p,p) - A(q,q))*DCONJG(A(p,p) + A(q,q)))
	          xw = DBLE((A(p,q) - A(q,p))*DCONJG(A(p,q) + A(q,p)))
	          dv = .5D0*(xv + xw)
	          dw = .5D0*(xv - xw)

	          tv = DCONJG(A(p,p))*A(q,p) + A(q,q)*DCONJG(A(p,q))
	          xv = sqrt(dv**2 + sq(tv))
	          tw = DCONJG(A(p,p))*A(p,q) + A(q,q)*DCONJG(A(q,p))
	          xw = sqrt(dw**2 + sq(tw))

	          t = sign(1D0, min(abs(dv + xv), abs(dw + xw)) -
     &                          min(abs(dv - xv), abs(dw - xw)))

	          tv = tv/(dv + t*xv)
	          invc = sqrt(1 + sq(tv))
	          sv = tv/invc
	          tv = tv/(invc + 1)

	          tw = tw/(dw + t*xw)
	          invc = sqrt(1 + sq(tw))
	          sw = tw/invc
	          tw = tw/(invc + 1)

	          x = A(p,p)
	          y = A(q,p)
	          evp = invc*(x + DCONJG(sv)*(y - tv*x))
	          x = A(p,q)
	          y = A(q,q)
	          evq = invc*(y - sv*(x + DCONJG(tv)*y))

	          do j = 1, nx
	            x = A(j,p)
	            y = A(j,q)
	            A(j,p) = x + DCONJG(sw)*(y - tw*x)
	            A(j,q) = y - sw*(x + DCONJG(tw)*y)
	            x = A(p,j)
	            y = A(q,j)
	            A(p,j) = x + DCONJG(sv)*(y - tv*x)
	            A(q,j) = y - sv*(x + DCONJG(tv)*y)
	          enddo

	          A(p,p) = evp
	          A(q,p) = 0
	          A(p,q) = 0
	          A(q,q) = evq

	          do j = 1, nx
	            x = V(j,p)
	            y = V(j,q)
	            V(j,p) = x + sv*(y - DCONJG(tv)*x)
	            V(j,q) = y - DCONJG(sv)*(x + tv*y)
	          enddo

	          do j = 1, nx
	            x = W(j,p)
	            y = W(j,q)
	            W(j,p) = x + sw*(y - DCONJG(tw)*x)
	            W(j,q) = y - DCONJG(sw)*(x + tw*y)
	          enddo
	        endif
	      endif
	    enddo
	  enddo
	enddo

	print *, "Bad convergence in SVD"

1	continue

	nm = min(m, n)
	rev = ibits(m - n, 15, 1)

* make the diagonal elements nonnegative

	do p = 1, nm
	  d(p) = abs(A(p,p))
	  if( d(p) .gt. DBL_EPS .and. d(p) .ne. DBLE(A(p,p)) ) then
	    f = A(p,p)/d(p)
	    do q = 1, nm
	      W(q,p) = W(q,p)*f
	    enddo
	  endif
	enddo

* sort the singular values

	do p = 1, nm
	  pi(p) = p
	enddo

	do p = 1, nm
	  j = p
	  t = d(p)
	  if( sort .ne. 0 ) then
	    do q = p + 1, nm
	      if( sort*(t - d(q)) .gt. 0 ) then
	        j = q
	        t = d(q)
	      endif
	    enddo
	  endif

	  d(j) = d(p)
	  d(p) = t

	  q = pi(j)
	  pi(j) = pi(p)

	  do j = 1, m
	    Element(Vo,p,j) = VW(j,q,rev)
	  enddo
	  do j = 1, n
	    Element(Wo,p,j) = VW(j,q,1-rev)
	  enddo
	enddo
	end

