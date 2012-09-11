* xsection.F
* routines to compute the cross-section
* this file is part of FormCalc
* last modified 30 Nov 11 th

#include "decl.h"

#include "process.h"

#include MODEL

#ifdef BREMSSTRAHLUNG
#include "softphoton.F"
#endif

************************************************************************
** ProcessIni translates the polarization string into bit-encoded
** helicities, determines the averaging factor, and initializes the
** model defaults.

	subroutine ProcessIni(fail, pol,
     &    sqrtSfrom, sqrtSto, sqrtSstep,
     &    argc, argv)
	implicit none
	integer fail, argc
	character*(*) pol
	character*128 argv(*)
	RealType sqrtSfrom, sqrtSto, sqrtSstep

#include "xsection.h"

#if U77EXT
#define len_trim lnblnk
	integer lnblnk
	external lnblnk
#endif

	integer i, c, bits, df

#define SCALAR Z'001'
#define FERMION Z'10A'
#define PHOTON Z'20A'
#define GLUON PHOTON
#define VECTOR Z'20E'
#define GRAVITINO Z'31B'
#define GRAVITON Z'41B'
#define TENSOR Z'41F'

	integer type(LEGS)
	data type /TYPES/

	if( len_trim(pol) .ne. LEGS ) then
	  fail = 1
	  return
	endif

	df = 2**(LEGS)

	do i = 1, LEGS
	  bits = iand(type(i), 255)
	  c = ior(ichar(pol(i:i)), 32)
	  if( c .eq. ichar('t') ) then
	    c = ior(c, ibits(bits, 3, 1))
	    bits = iand(bits, 16+8+2+1)
	  else if( c .eq. ichar('p') ) then
	    bits = iand(bits, 16)
	  else if( c .eq. ichar('+') ) then
	    bits = iand(bits, 8)
	  else if( c .eq. ichar('l') ) then
	    bits = iand(bits, 4)
	  else if( c .eq. ichar('-') ) then
	    bits = iand(bits, 2)
	  else if( c .eq. ichar('m') ) then
	    bits = iand(bits, 1)
	  else if( c .ne. ichar('u') ) then
	    bits = 0
	  endif

	  if( bits .eq. 0 ) then
	    Error(fail, "Invalid polarization for leg "/ /Digit(i))
	    return
	  endif

	  if( i .le. LEGS_IN ) df = df*(ibits(bits, 4, 1) +
     &      ibits(bits, 3, 1) + ibits(bits, 2, 1) +
     &      ibits(bits, 1, 1) + ibits(bits, 0, 1))

#ifdef DIRACFERMIONS
	  if( type(i) .eq. FERMION ) then
	    df = df/2
	    bits = 1
	  endif
#endif

	  pol(i:i) = char(c)
	  helicities = helicities*32 + bits
	enddo

	Lower(SQRTS) = sqrtSfrom
	Upper(SQRTS) = sqrtSto
	Step(SQRTS) = sqrtSstep

	Var(FIXED) = 0
	Step(FIXED) = 1

	Var(TRIVIAL) = 0
	Step(TRIVIAL) = 0

	avgfac = 2**(LEGS)*Re(COLOURFACTOR)*Re(IDENTICALFACTOR)/df
	sqrtS = -1
	threshold = -1
	scale = -1
	sqrtSinvalid = 1

	call ltini

#ifdef CUTTOOLS
* args are:
* 1. rootsvalue = arbitrary OPP scale of the order of the scale of the process
* 2. limitvalue = upper limit of precision above which the mp routines activate
* 3. idig       = max number of digits of the mp routines
* 4. scaloop    = scalar integrals (1 is LoopTools)
* 5. muscale    = the renormalization scale
* 6. thrs       = internal collinear IR threshold for avh_olo
	rootsvalue = sqrtSfrom
	muscale = getmudim()
	call ctsinit(rootsvalue, 1D-2, 0, 1, muscale, 1D-6)
#endif

#ifdef SAMURAI
* args are:
* 1. imeth = "diag" (numerators) or "tree" (products of tree amps)
* 2. isca = 1 (QCDloop) or 2 (OneLOop)
* 3. verbosity = 0, 1, 2, 3
* 4. itest = 0 (none), 1 (powertest), 2 (nntest), 3 (lnntest)
	call initsamurai("diag", 2, 1, 1)
#endif

	call ModelDefaults(argc, argv)
	call LumiDefaults(argc, argv)

	fail = 0
	end


************************************************************************
** ProcessExi wraps up the calculation, e.g. prints a summary of
** messages, deallocates arrays etc.

	subroutine ProcessExi
	implicit none

	call ltexi
	end


************************************************************************
** SetEnergy sets the energy for the partonic scattering process.
** All scale-dependent quantities are initialized at SCALE.

	subroutine SetEnergy(fail, newsqrtS)
	implicit none
	integer fail
	RealType newsqrtS

#include "xsection.h"

	RealType oldscale, oldmass_in, oldmass_out
	logical reset
	integer i

	fail = 0
	reset = abs(newsqrtS - sqrtS) .gt. 1D-9
	sqrtS = newsqrtS
	oldscale = max(Re(SCALE), 1D0)

	if( sqrtSinvalid .eq. 1 .or.
     &      btest(flags, BIT_RESET) .or.
     &      abs(oldscale - scale) .gt. 1D-9 ) then
	  oldmass_in = MASS_IN
	  oldmass_out = MASS_OUT

	  do i = 1, 10
	    call ModelVarIni(fail, oldscale)

	    scale = max(Re(SCALE), 1D0)
	    mass_in = MASS_IN
	    mass_out = MASS_OUT

	    if( abs(scale - oldscale) +
     &          abs(mass_in - oldmass_in) +
     &          abs(mass_out - oldmass_out) .lt. 1D-9 ) goto 1

	    oldscale = scale
	    oldmass_in = mass_in
	    oldmass_out = mass_out
	  enddo

1	  threshold = max(mass_in, mass_out)
	  sqrtSinvalid = fail
	  reset = .TRUE.
	endif

	if( reset ) flags = ibset(flags, BIT_RESET)

	if( fail .ne. 0 .or. sqrtS .lt. threshold ) then
	  fail = 1
	else
	  call InitialState
	endif
	end


************************************************************************
** ParameterScan performs the actual calculation.
** It scans over the parameter loops LOOP* declared in run.F.

	subroutine ParameterScan(dir,
     &    serialfrom, serialto, serialstep)
	implicit none
	character*(*) dir
	integer serialfrom, serialto, serialstep

#include "xsection.h"

	integer openlog
	external openlog

	integer serial, next

	next = serialfrom
	serial = 0

#ifdef LOOP1
	LOOP1
#endif
#ifdef LOOP2
	LOOP2
#endif
#ifdef LOOP3
	LOOP3
#endif
#ifdef LOOP4
	LOOP4
#endif
#ifdef LOOP5
	LOOP5
#endif
#ifdef LOOP6
	LOOP6
#endif
#ifdef LOOP7
	LOOP7
#endif
#ifdef LOOP8
	LOOP8
#endif
#ifdef LOOP9
	LOOP9
#endif
#ifdef LOOP10
	LOOP10
#endif
#ifdef LOOP11
	LOOP11
#endif
#ifdef LOOP12
	LOOP12
#endif
#ifdef LOOP13
	LOOP13
#endif
#ifdef LOOP14
	LOOP14
#endif
#ifdef LOOP15
	LOOP15
#endif
#ifdef LOOP16
	LOOP16
#endif
#ifdef LOOP17
	LOOP17
#endif
#ifdef LOOP18
	LOOP18
#endif
#ifdef LOOP19
	LOOP19
#endif
#ifdef LOOP20
	LOOP20
#endif

	serial = serial + 1
	if( serial .lt. next ) goto 1

	call flush(6)
	if( openlog(dir, serial) .eq. 0 ) then
	  call IntegratedCS
	  call flush(6)
	  call closelog
	endif

	next = next + serialstep
	if( next .gt. serialto ) return

1	continue
	end


************************************************************************
** IntegratedCS computes the integrated cross-section at one point in
** parameter space.  The cross-section is differential in all variables
** with non-zero step, i.e. integration variables have step = 0.

	subroutine IntegratedCS
	implicit none

#include "xsection.h"

	integer nvars
	parameter (nvars = MAXVAR - (MINVAR) + 1)

	external DifferentialCS

	RealType result(NCOMP), error(NCOMP), show(nvars)
	integer fail, f, v, nfix, ndim, fix(nvars)

	call ModelConstIni(fail)
	if( fail .ne. 0 ) goto 999

	flags = 2**BIT_RESET + 2**BIT_LOOP

#ifdef MUDIM
	call setmudim(Re(MUDIM))
#endif
#ifdef DELTA
	call setdelta(Re(DELTA))
#endif
#ifdef LAMBDA
	call setlambda(Re(LAMBDA))
#endif
	Divergence = getdelta()
	mudim = getmudim()
	lambda = getlambda()
	if( lambda .eq. -1 ) then
	  epscoeff = -1
	  Finite = 0
	else if( lambda .eq. -2 ) then
	  epscoeff = -2
	  Finite = 0
	else
	  epscoeff = 0
	  Finite = 1
	endif

	call KinIni(fail)
	if( fail .ne. 0 ) goto 999

	nfix = 0
	do v = MINVAR, MAXVAR
	  if( Step(v) .ne. 0 ) then
	    nfix = nfix + 1
	    fix(nfix) = v
	    Var(v) = Lower(v)
	  endif
	enddo

	ndim = nvars - nfix

	call LumiIni(fail)
	if( fail .ne. 0 ) goto 999

	call ModelDigest

#define SHOW print 100,
100	format("|# ", A, "=", F10.4, SP, F10.4, " I")

#ifdef MMA
	call MmaSetPara
#endif

#ifdef PRINT1
	PRINT1
#endif
#ifdef PRINT2
	PRINT2
#endif
#ifdef PRINT3
	PRINT3
#endif
#ifdef PRINT4
	PRINT4
#endif
#ifdef PRINT5
	PRINT5
#endif
#ifdef PRINT6
	PRINT6
#endif
#ifdef PRINT7
	PRINT7
#endif
#ifdef PRINT8
	PRINT8
#endif
#ifdef PRINT9
	PRINT9
#endif
#ifdef PRINT10
	PRINT10
#endif
#ifdef PRINT11
	PRINT11
#endif
#ifdef PRINT12
	PRINT12
#endif
#ifdef PRINT13
	PRINT13
#endif
#ifdef PRINT14
	PRINT14
#endif
#ifdef PRINT15
	PRINT15
#endif
#ifdef PRINT16
	PRINT16
#endif
#ifdef PRINT17
	PRINT17
#endif
#ifdef PRINT18
	PRINT18
#endif
#ifdef PRINT19
	PRINT19
#endif
#ifdef PRINT20
	PRINT20
#endif

#ifdef MMA
	call MmaEndSet
	call MmaSetData
#endif

1	call Cuba(ndim, DifferentialCS, result, error)

	do f = 1, nfix
	  show(f) = Show(fix(f))
	enddo

#ifdef MMA
	call MmaData(show, nfix, result, error, NCOMP)
#else
* Note: "real" data lines are tagged with "|" in the output.
101	format("| ", 10(4G19.10, :, /"|+"))
	print 101, (show(f), f = 1, nfix)

102	format("|+   ", NCOMP G24.15)
	print 102, result
	print 102, error

	call flush(6)
#endif

	do f = nfix, 1, -1
	  v = fix(f)
	  Var(v) = Var(v) + Step(v)
	  if( (Var(v) - Upper(v))/Step(v) .lt. 1D-10 ) goto 1
	  Var(v) = Lower(v)
	enddo

#ifdef MMA
	call MmaEndSet
#else
103	format("|"/"|"/)
	print 103
#endif

999	continue
	end


************************************************************************
** DifferentialCS computes the differential cross-section at x.
** For all integration variables (those with zero step) it factors in
** the Jacobian, too.

	subroutine DifferentialCS(ndim, x, ncomp, result)
	implicit none
	integer ndim, ncomp
	RealType x(ndim), result(ncomp)

#include "xsection.h"

#ifdef BREMSSTRAHLUNG
	RealType SoftPhotonFactor
	external SoftPhotonFactor
#endif

	RealType fac, range, flux
	integer v, d, c

	fac = avgfac
	d = 0
	do v = MINVAR, MAXVAR
	  if( Step(v) .eq. 0 ) then
	    range = Upper(v) - Lower(v)
	    d = d + 1
	    Var(v) = Lower(v) + x(d)*range
	    fac = fac*range
	  endif
	  Show(v) = Var(v)
	enddo

	do c = 1, ncomp
	  result(c) = 0
	enddo

	call Luminosity(fac)
	if( fac .eq. 0 ) return

	call FinalState(fac)
	if( fac .eq. 0 ) return

	if( btest(flags, BIT_RESET) ) then
	  call clearcache
#ifndef NO_RENCONST
	  call CalcRenConst
#endif
	endif

	call SquaredME(result, helicities, flags)

	flags = ibclr(flags, BIT_RESET)

#if LEGS_IN == 1
	flux = 2*sqrtS
#else
	flux = 4/hbar_c2*sqrtS*momspec(SPEC_K, 1)
#endif
	fac = fac/((2*pi)**(3*LEGS_OUT - 4)*2*sqrtS*flux)

	do c = 1, ncomp
	  if( .not. abs(result(c)) .lt. 1D16 ) then
	    Warning("Got strange values from SquaredME:")
	    INFO result
	    INFO "(Did you compute the colour matrix elements?)"
	    stop
	  endif
	  result(c) = result(c)*fac
	enddo

#ifdef BREMSSTRAHLUNG
	result(2) = result(2) + SoftPhotonFactor()*result(1)
#endif

#ifdef WF_RENORMALIZATION
	result(2) = result(2) + (WF_RENORMALIZATION)*result(1)
#endif
	end


************************************************************************
** Cuba is a chooser for the Cuba routines, with special cases
** for ndim = 0 (integrand evaluation) and ndim = 1 (Patterson
** integration).

	subroutine Cuba(ndim, integrand, result, error)
	implicit none
	integer ndim
	external integrand
	RealType result(NCOMP), error(NCOMP)

	integer nregions, neval, fail, c
	RealType prob(NCOMP)

#include "xsection.h"

#define GAUSS 1
#define PATTERSON 2
#define VEGAS 3
#define SUAVE 4
#define DIVONNE 5
#define CUHRE 6

#define FLAGS_LAST 4
#define FLAGS_PSEUDO 8

	if( ndim .eq. 0 ) then

	  call integrand(0, 0D0, NCOMP, result)
	  do c = 1, NCOMP
	    error(c) = 0
	  enddo
	  return

	else if( ndim .eq. 1 ) then

#if METHOD == GAUSS

	  neval = 32
	  call Gauss(NCOMP, 0D0, 1D0, integrand,
     &      neval, result)
	  do c = 1, NCOMP
	    error(c) = -1
	    prob(c) = -1
	  enddo
	  nregions = 1
	  fail = 0
	  INFO "Gauss integration results:"

#else

	  call Patterson(NCOMP, 0D0, 1D0, integrand,
     &      Re(RELACCURACY), Re(ABSACCURACY),
     &      neval, fail, result, error)
	  do c = 1, NCOMP
	    prob(c) = -1
	  enddo
	  nregions = 1
	  INFO "Patterson integration results:"

#endif

	else

#if METHOD == VEGAS

	  call vegas(ndim, NCOMP, integrand, USERDATA,
     &      Re(RELACCURACY), Re(ABSACCURACY),
     &      VERBOSE, SEED, MINEVAL, MAXEVAL,
     &      NSTART, NINCREASE, NBATCH,
     &      GRIDNO, STATEFILE,
     &      neval, fail, result, error, prob)
	  nregions = 1
	  INFO "Vegas integration results:"

#elif METHOD == SUAVE

	  call suave(ndim, NCOMP, integrand, USERDATA,
     &      Re(RELACCURACY), Re(ABSACCURACY),
     &      VERBOSE + FLAGS_LAST, SEED, MINEVAL, MAXEVAL,
     &      NNEW, Re(FLATNESS),
     &      nregions, neval, fail, result, error, prob)
	  INFO "Suave integration results:"

#elif METHOD == DIVONNE

	  call divonne(ndim, NCOMP, integrand, USERDATA,
     &      Re(RELACCURACY), Re(ABSACCURACY),
     &      VERBOSE, SEED, MINEVAL, MAXEVAL,
     &      KEY1, KEY2, KEY3, MAXPASS,
     &      Re(BORDER), Re(MAXCHISQ), Re(MINDEVIATION),
     &      0, NDIM, 0, 0, 0,
     &      nregions, neval, fail, result, error, prob)
	  INFO "Divonne integration results:"

#else

	  call cuhre(ndim, NCOMP, integrand, USERDATA,
     &      Re(RELACCURACY), Re(ABSACCURACY),
     &      VERBOSE + FLAGS_LAST, MINEVAL, MAXEVAL,
     &      KEY,
     &      nregions, neval, fail, result, error, prob)
	  INFO "Cuhre integration results:"

#endif

	endif

	INFO "nregions =", nregions
	INFO "neval    =", neval
	INFO "fail     =", fail

	if( fail .gt. 0 ) then
	  Warning("Failed to reach the desired accuracy")
	else if( fail .lt. 0 ) then
	  Error(fail, "Integration error")
	endif

100	format(I2, G24.15, " +- ", G24.15, "    p = ", F6.3)
	print 100, (c, result(c), error(c), prob(c), c = 1, NCOMP)
	end

