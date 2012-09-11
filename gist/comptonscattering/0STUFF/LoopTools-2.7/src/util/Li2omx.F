* Li2omx.F
* the dilogarithm function of 1 - x
* this file is part of LoopTools
* last modified 25 Feb 11 th

#include "defs.h"


	double complex function XLi2omx(x)
	implicit none
	DVAR x

	double complex spence1

	XLi2omx = spence1(DCMPLX(1 - x), DCMPLX(x), 0D0)
	end

************************************************************************
* adapter code for C++

	subroutine XLi2omxsub(res, x)
	implicit none
	double complex res
	DVAR x

	double complex spence1
	external spence1

	res = spence1(DCMPLX(1 - x), DCMPLX(x), 0D0)
        end

