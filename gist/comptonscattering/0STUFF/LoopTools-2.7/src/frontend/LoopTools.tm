:Evaluate: BeginPackage["LoopTools`"]

:Evaluate: A0::usage = "A0[m] is the one-point one-loop scalar integral.
	m is the mass squared."

:Evaluate: A00::usage = "A00[m] is the one-point one-loop tensor coefficient of g_{mu nu}.
	m is the mass squared."

:Evaluate: B0i::usage = "B0i[id, p, m1, m2] is the generic two-point loop integral which includes both scalar and tensor coefficients, as well as certain derivatives.
	For example, B0i[bb0, ...] is the scalar function B_0, B0i[bb11, ...] the tensor coefficient function B_11 etc.
	p is the external momentum squared and m1 and m2 are the masses squared."

:Evaluate: Bget::usage = "Bget[p, m1, m2] returns a list of all two-point coefficients."

:Evaluate: C0::usage = "C0[p1, p2, p1p2, m1, m2, m3] is the three-point one-loop scalar integral.
	p1, p2, and p1p2 are the external momenta squared and m1, m2, m3 are the masses squared."

:Evaluate: C0i::usage = "C0i[id, p1, p2, p1p2, m1, m2, m3] is the generic three-point loop integral which includes both scalar and tensor coefficients, specified by id.
	For example, C0i[cc0, ...] is the scalar function C_0, C0i[cc112, ...] the tensor coefficient function C_112 etc.
	p1, p2, and p1p2 are the external momenta squared and m1, m2, m3 are the masses squared."

:Evaluate: Cget::usage = "Cget[p1, p2, p1p2, m1, m2, m3] returns a list of all three-point coefficients."

:Evaluate: D0::usage = "D0[p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4] is the four-point scalar one-loop integral.
	p1...p4 are the external momenta squared, p1p2 and p2p3 are the squares of external momenta 1 + 2 and 2 + 3, respectively, and m1...m4 are the masses squared."

:Evaluate: D0i::usage = "D0i[id, p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4] is the generic four-point loop integral which includes both scalar and tensor coefficients, specified by id.
	For example, D0i[dd0, ...] is the scalar function D_0, D0i[dd1233, ...] the tensor function D_{1233} etc.
	p1...p4 are the external momenta squared, p1p2 and p2p3 are the squares of external momenta 1 + 2 and 2 + 3, respectively, and m1...m4 are the masses squared."

:Evaluate: Dget::usage = "Dget[p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4] returns a list of all four-point coefficients."

:Evaluate: E0::usage = "E0[p1, p2, p3, p4, p1p2, p2p3, p3p4, p4p5, p5p1, m1, m2, m3, m4, m5] is the five-point scalar one-loop integral.
	p1...p5 are the external momenta squared, pipj are the squares of (pi + pj), and m1...m5 are the masses squared."

:Evaluate: E0i::usage = "E0i[id, p1, p2, p3, p4, p5, p1p2, p2p3, p3p4, p4p5, p5p1, m1, m2, m3, m4, m5] is the generic five-point loop integral which includes both scalar and tensor coefficients, specified by id.
	For example, E0i[ee0, ...] is the scalar function E_0, E0i[ee3444, ...] the tensor function E_{3444} etc.
	p1...p5 are the external momenta squared, pipj are the squares of (pi + pj), and m1...m5 are the masses squared."

:Evaluate: Eget::usage = "Eget[p1, p2, p3, p4, p5, p1p2, p2p3, p3p4, p4p5, p5p1, m1, m2, m3, m4, m5] returns a list of all five-point coefficients."

:Evaluate: PaVe::usage = "PaVe[ind, {pi}, {mi}] is the generalized Passarino-Veltman function used by FeynCalc.
	It is converted to B0i, C0i, D0i, or E0i in LoopTools."

:Evaluate: Li2::usage = "Li2[x] returns the dilogarithm of x."

:Evaluate: Li2omx::usage = "Li2omx[x] returns the dilogarithm of 1 - x."

:Evaluate: SetMudim::usage = "SetMudim[m^2] sets the renormalization scale squared."

:Evaluate: GetMudim::usage = "GetMudim[] returns the current value for the renormalization scale squared."

:Evaluate: SetDelta::usage = "SetDelta[d] sets the numerical value of Delta which replaces the divergence 2/(4 - D) - EulerGamma + Log[4 Pi] in LoopTools."

:Evaluate: GetDelta::usage = "GetDelta[] returns the current numerical value of Delta which replaces the divergence 2/(4 - D) - EulerGamma + Log[4 Pi] in LoopTools."

:Evaluate: SetLambda::usage = "SetLambda[l^2] sets the infrared regulator mass squared."

:Evaluate: GetLambda::usage = "GetLambda[] returns the current value for the infrared regulator mass squared."

:Evaluate: SetMinMass::usage = "SetMinMass[m^2] sets the collinear cutoff mass."

:Evaluate: GetMinMass::usage = "GetMinMass[] returns the current value for the collinear cutoff mass squared."

:Evaluate: ClearCache::usage = "ClearCache[] clears the internal LoopTools caches."

:Evaluate: MarkCache::usage = "MarkCache[] marks the current positions of the internal LoopTools caches."

:Evaluate: RestoreCache::usage = "RestoreCache[] restores the internal LoopTools caches to the position when the last MarkCache was issued."

:Evaluate: SetMaxDev::usage = "SetMaxDev[d] sets the maximum relative deviation a result and its alternate derivation may have before a warning is issued."

:Evaluate: GetMaxDev::usage = "GetMaxDev[d] returns the maximum relative deviation a result and its alternate derivation may have before a warning is issued."

:Evaluate: SetWarnDigits::usage = "SetWarnDigits[n] sets the number of LoopTools' warning digits.
	If the number of digits presumed lost by FF is larger than the warning digits, either an alternate version is tried (if available) or a warning is issued."

:Evaluate: GetWarnDigits::usage = "GetWarnDigits[] returns the number of LoopTools' warning digits.
	If the number of digits presumed lost by FF is larger than the warning digits, either an alternate version is tried (if available) or a warning is issued."

:Evaluate: SetErrDigits::usage = "SetErrDigits[n] sets the number of LoopTools' error digits.
	If the number of digits presumed lost by FF is larger than the error digits, the alternate result is used instead of the FF result."

:Evaluate: GetErrDigits::usage = "GetErrDigits[] returns the number of LoopTools' error digits.
	If the number of digits presumed lost by FF is larger than the error digits, the alternate result is used instead of the FF result."

:Evaluate: SetVersionKey::usage = "SetVersionKey[key] sets the LoopTools version key.
	It determines which version of a loop integral is returned, and whether checks are performed."

:Evaluate: GetVersionKey::usage = "GetVersionKey[] returns the LoopTools version key.
	It determines which version of a loop integral is returned, and whether checks are performed."

:Evaluate: SetDebugKey::usage = "SetDebugKey[key] sets the LoopTools debug key.
	It determines how much debug information is printed for a loop integral."

:Evaluate: GetDebugKey::usage = "GetDebugKey[] returns the LoopTools debug key.
	It determines how much debug information is printed for a loop integral."

:Evaluate: SetDebugRange::usage = "SetDebugRange[from, to] sets the LoopTools debug range.
	The integrals printed out on screen as determined by the debug key are numbered consecutively.
	Setting a debug range restricts printing to the given range."

:Evaluate: SetCmpBits::usage = "SetCmpBits[bits] sets the number of bits compared in cache lookups.
	Setting it to less than 64 (double precision) makes the comparison more robust against numerical noise."

:Evaluate: GetCmpBits::usage = "GetCmpBits[] returns the number of bits compared of each real number in cache lookups."

:Evaluate: PaVe::usage = "PaVe[ind, {pi}, {mi}] is the generalized Passarino-Veltman function used by FeynCalc.
	It is converted to B0i, C0i, D0i, or E0i in LoopTools."

:Evaluate: Bids = {bb0, bb1, bb00, bb11, bb001, bb111,
	  dbb0, dbb1, dbb00, dbb11};
	Cids = {cc0, cc1, cc2, cc00, cc11, cc12, cc22,
	  cc001, cc002,	cc111, cc112, cc122, cc222,
	  cc0000, cc0011, cc0012, cc0022, cc1111, cc1112,
	  cc1122, cc1222, cc2222};
	Dids = {dd0, dd1, dd2, dd3, dd00, dd11, dd12, dd13, dd22, dd23,
	  dd33, dd001, dd002, dd003, dd111, dd112, dd113, dd122, dd123,
	  dd133, dd222, dd223, dd233, dd333, dd0000, dd0011, dd0012,
	  dd0013, dd0022, dd0023, dd0033, dd1111, dd1112, dd1113,
	  dd1122, dd1123, dd1133, dd1222, dd1223, dd1233, dd1333,
	  dd2222, dd2223, dd2233, dd2333, dd3333, dd00001, dd00002,
	  dd00003, dd00111, dd00112, dd00113, dd00122, dd00123, dd00133,
	  dd00222, dd00223, dd00233, dd00333, dd11111, dd11112, dd11113,
	  dd11122, dd11123, dd11133, dd11222, dd11223, dd11233, dd11333,
	  dd12222, dd12223, dd12233, dd12333, dd13333, dd22222, dd22223,
	  dd22233, dd22333, dd23333, dd33333};
	Eids = {ee0, ee1, ee2, ee3, ee4, ee00, ee11, ee12, ee13, ee14,
	  ee22, ee23, ee24, ee33, ee34, ee44, ee001, ee002, ee003,
	  ee004, ee111, ee112, ee113, ee114, ee122, ee123, ee124,
	  ee133, ee134, ee144, ee222, ee223, ee224, ee233, ee234,
	  ee244, ee333, ee334, ee344, ee444, ee0000, ee0011, ee0012,
	  ee0013, ee0014, ee0022, ee0023, ee0024, ee0033, ee0034,
	  ee0044, ee1111, ee1112, ee1113, ee1114, ee1122, ee1123,
	  ee1124, ee1133, ee1134, ee1144, ee1222, ee1223, ee1224,
	  ee1233, ee1234, ee1244, ee1333, ee1334, ee1344, ee1444,
	  ee2222, ee2223, ee2224, ee2233, ee2234, ee2244, ee2333, 
	  ee2334, ee2344, ee2444, ee3333, ee3334, ee3344, ee3444,
	  ee4444}

:Evaluate: KeyAll = Plus@@ ({KeyA0, KeyBget, KeyC0, KeyD0, KeyE0,
	  KeyEget, KeyCEget} = 4^Range[0, 6]);
	DebugAll = Plus@@ ({DebugB, DebugC, DebugD, DebugE} =
	  2^Range[0, 3])

:Evaluate: B0 = B0i[bb0, ##]&;
	B1 = B0i[bb1, ##]&;
	B00 = B0i[bb00, ##]&;
	B11 = B0i[bb11, ##]&;
	B001 = B0i[bb001, ##]&;
	B111 = B0i[bb111, ##]&;
	DB0 = B0i[dbb0, ##]&;
	DB1 = B0i[dbb1, ##]&;
	DB00 = B0i[dbb00, ##]&;
	DB11 = B0i[dbb11, ##]&

:Evaluate: Begin["`Private`"]

:Begin:
:Function:	mA0
:Pattern:	A0[m_?r]
:Arguments:	{N[m]}
:ArgumentTypes:	{Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mA0c
:Pattern:	A0[m_?c]
:Arguments:	{N[Re[m]], N[Im[m]]}
:ArgumentTypes:	{Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mA00
:Pattern:	A00[m_?r]
:Arguments:	{N[m]}
:ArgumentTypes:	{Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mA00c
:Pattern:	A00[m_?c]
:Arguments:	{N[Re[m]], N[Im[m]]}
:ArgumentTypes:	{Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mB0i
:Pattern:	B0i[id_, p_?r, m1_?r, m2_?r]
:Arguments:	{Bid[id], N[p], N[m1], N[m2]}
:ArgumentTypes:	{Integer, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mB0ic
:Pattern:	B0i[id_, p_?c, m1_?c, m2_?c]
:Arguments:	{Bid[id], N[Re[p]], N[Im[p]],
		  N[Re[m1]], N[Im[m1]], N[Re[m2]], N[Im[m2]]}
:ArgumentTypes:	{Integer, Real, Real, Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mBget
:Pattern:	Bget[p_?r, m1_?r, m2_?r]
:Arguments:	{N[p], N[m1], N[m2]}
:ArgumentTypes:	{Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mBgetc
:Pattern:	Bget[p_?c, m1_?c, m2_?c]
:Arguments:	{N[Re[p]], N[Im[p]],
		  N[Re[m1]], N[Im[m1]], N[Re[m2]], N[Im[m2]]}
:ArgumentTypes:	{Real, Real, Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mC0
:Pattern:	C0[p1_?r, p2_?r, p1p2_?r, m1_?r, m2_?r, m3_?r]
:Arguments:	{N[p1], N[p2], N[p1p2], N[m1], N[m2], N[m3]}
:ArgumentTypes:	{Real, Real, Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mC0c
:Pattern:	C0[p1_?c, p2_?c, p1p2_?c, m1_?c, m2_?c, m3_?c]
:Arguments:	{N[Re[p1]], N[Im[p1]], N[Re[p2]], N[Im[p2]],
		  N[Re[p1p2]], N[Im[p1p2]], N[Re[m1]], N[Im[m1]],
		  N[Re[m2]], N[Im[m2]], N[Re[m3]], N[Im[m3]]}
:ArgumentTypes:	{Real, Real, Real, Real, Real, Real,
		  Real, Real, Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mC0i
:Pattern:	C0i[id_, p1_?r, p2_?r, p1p2_?r, m1_?r, m2_?r, m3_?r]
:Arguments:	{Cid[id],
		  N[p1], N[p2], N[p1p2], N[m1], N[m2], N[m3]}
:ArgumentTypes:	{Integer, Real, Real, Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mC0ic
:Pattern:	C0i[id_, p1_?c, p2_?c, p1p2_?c, m1_?c, m2_?c, m3_?c]
:Arguments:	{Cid[id],
		  N[Re[p1]], N[Im[p1]], N[Re[p2]], N[Im[p2]],
		  N[Re[p1p2]], N[Im[p1p2]], N[Re[m1]], N[Im[m1]],
		  N[Re[m2]], N[Im[m2]], N[Re[m3]], N[Im[m3]]}
:ArgumentTypes:	{Integer,
		  Real, Real, Real, Real, Real, Real,
		  Real, Real, Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mCget
:Pattern:	Cget[p1_?r, p2_?r, p1p2_?r, m1_?r, m2_?r, m3_?r]
:Arguments:	{N[p1], N[p2], N[p1p2], N[m1], N[m2], N[m3]}
:ArgumentTypes:	{Real, Real, Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mCgetc
:Pattern:	Cget[p1_?c, p2_?c, p1p2_?c, m1_?c, m2_?c, m3_?c]
:Arguments:	{N[Re[p1]], N[Im[p1]], N[Re[p2]], N[Im[p2]],
		  N[Re[p1p2]], N[Im[p1p2]], N[Re[m1]], N[Im[m1]],
		  N[Re[m2]], N[Im[m2]], N[Re[m3]], N[Im[m3]]}
:ArgumentTypes:	{Real, Real, Real, Real, Real, Real,
		  Real, Real, Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mD0
:Pattern:	D0[p1_?r, p2_?r, p3_?r, p4_?r, p1p2_?r, p2p3_?r,
		  m1_?r, m2_?r, m3_?r, m4_?r]
:Arguments:	{N[p1], N[p2], N[p3], N[p4], N[p1p2], N[p2p3],
		  N[m1], N[m2], N[m3], N[m4]}
:ArgumentTypes:	{Real, Real, Real, Real, Real, Real,
		  Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mD0c
:Pattern:	D0[p1_?c, p2_?c, p3_?c, p4_?c, p1p2_?c, p2p3_?c,
		  m1_?c, m2_?c, m3_?c, m4_?c]
:Arguments:	{N[Re[p1]], N[Im[p1]], N[Re[p2]], N[Im[p2]],
		  N[Re[p3]], N[Im[p3]], N[Re[p4]], N[Im[p4]],
		  N[Re[p1p2]], N[Im[p1p2]], N[Re[p2p3]], N[Im[p2p3]],
		  N[Re[m1]], N[Im[m1]], N[Re[m2]], N[Im[m2]],
		  N[Re[m3]], N[Im[m3]], N[Re[m4]], N[Im[m4]]}
:ArgumentTypes:	{Real, Real, Real, Real, Real, Real, Real, Real,
		  Real, Real, Real, Real, Real, Real, Real, Real,
		  Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mD0i
:Pattern:	D0i[id_, p1_?r, p2_?r, p3_?r, p4_?r, p1p2_?r, p2p3_?r,
		  m1_?r, m2_?r, m3_?r, m4_?r]
:Arguments:	{Did[id],
		  N[p1], N[p2], N[p3], N[p4], N[p1p2], N[p2p3],
		  N[m1], N[m2], N[m3], N[m4]}
:ArgumentTypes:	{Integer,
		  Real, Real, Real, Real, Real, Real,
		  Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mD0ic
:Pattern:	D0i[id_, p1_?c, p2_?c, p3_?c, p4_?c, p1p2_?c, p2p3_?c,
		  m1_?c, m2_?c, m3_?c, m4_?c]
:Arguments:	{Did[id],
		  N[Re[p1]], N[Im[p1]], N[Re[p2]], N[Im[p2]],
		  N[Re[p3]], N[Im[p3]], N[Re[p4]], N[Im[p4]],
		  N[Re[p1p2]], N[Im[p1p2]], N[Re[p2p3]], N[Im[p2p3]],
		  N[Re[m1]], N[Im[m1]], N[Re[m2]], N[Im[m2]],
		  N[Re[m3]], N[Im[m3]], N[Re[m4]], N[Im[m4]]}
:ArgumentTypes:	{Integer,
		  Real, Real, Real, Real, Real, Real, Real, Real,
		  Real, Real, Real, Real, Real, Real, Real, Real,
		  Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mDget
:Pattern:	Dget[p1_?r, p2_?r, p3_?r, p4_?r, p1p2_?r, p2p3_?r,
		  m1_?r, m2_?r, m3_?r, m4_?r]
:Arguments:	{N[p1], N[p2], N[p3], N[p4], N[p1p2], N[p2p3],
		  N[m1], N[m2], N[m3], N[m4]}
:ArgumentTypes:	{Real, Real, Real, Real, Real, Real, Real,
		  Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mDgetc
:Pattern:	Dget[p1_?c, p2_?c, p3_?c, p4_?c, p1p2_?c, p2p3_?c,
		  m1_?c, m2_?c, m3_?c, m4_?c]
:Arguments:	{N[Re[p1]], N[Im[p1]], N[Re[p2]], N[Im[p2]],
		  N[Re[p3]], N[Im[p3]], N[Re[p4]], N[Im[p4]],
		  N[Re[p1p2]], N[Im[p1p2]], N[Re[p2p3]], N[Im[p2p3]],
		  N[Re[m1]], N[Im[m1]], N[Re[m2]], N[Im[m2]],
		  N[Re[m3]], N[Im[m3]], N[Re[m4]], N[Im[m4]]}
:ArgumentTypes:	{Real, Real, Real, Real, Real, Real, Real, Real,
		  Real, Real, Real, Real, Real, Real, Real, Real,
		  Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mE0
:Pattern:	E0[p1_?r, p2_?r, p3_?r, p4_?r, p5_?r,
		  p1p2_?r, p2p3_?r, p3p4_?r, p4p5_?r, p5p1_?r,
		  m1_?r, m2_?r, m3_?r, m4_?r, m5_?r]
:Arguments:	{N[p1], N[p2], N[p3], N[p4], N[p5], 
		  N[p1p2], N[p2p3], N[p3p4], N[p4p5], N[p5p1],
		  N[m1], N[m2], N[m3], N[m4], N[m5]}
:ArgumentTypes:	{Real, Real, Real, Real, Real, 
		  Real, Real, Real, Real, Real,
		  Real, Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mE0c
:Pattern:	E0[p1_?c, p2_?c, p3_?c, p4_?c, p5_?c,
		  p1p2_?c, p2p3_?c, p3p4_?c, p4p5_?c, p5p1_?c,
		  m1_?c, m2_?c, m3_?c, m4_?c, m5_?c]
:Arguments:	{N[Re[p1]], N[Im[p1]], N[Re[p2]], N[Im[p2]],
		  N[Re[p3]], N[Im[p3]], N[Re[p4]], N[Im[p4]],
		  N[Re[p5]], N[Im[p5]], N[Re[p1p2]], N[Im[p1p2]], 
		  N[Re[p2p3]], N[Im[p2p3]], N[Re[p3p4]], N[Im[p3p4]],
		  N[Re[p4p5]], N[Im[p4p5]], N[Re[p5p1]], N[Im[p5p1]],
		  N[Re[m1]], N[Im[m1]], N[Re[m2]], N[Im[m2]],
		  N[Re[m3]], N[Im[m3]], N[Re[m4]], N[Im[m4]],
		  N[Re[m5]], N[Im[m5]]}
:ArgumentTypes:	{Real, Real, Real, Real, Real, 
		 Real, Real, Real, Real, Real, 
		 Real, Real, Real, Real, Real, 
		 Real, Real, Real, Real, Real, 
		 Real, Real, Real, Real, Real, 
		 Real, Real, Real, Real, Real} 
:ReturnType:	Manual
:End:

:Begin:
:Function:	mE0i
:Pattern:	E0i[id_, p1_?r, p2_?r, p3_?r, p4_?r, p5_?r,
		  p1p2_?r, p2p3_?r, p3p4_?r, p4p5_?r, p5p1_?r,
		  m1_?r, m2_?r, m3_?r, m4_?r, m5_?r]
:Arguments:	{Eid[id],
		  N[p1], N[p2], N[p3], N[p4], N[p5], 
		  N[p1p2], N[p2p3], N[p3p4], N[p4p5], N[p5p1],
		  N[m1], N[m2], N[m3], N[m4], N[m5]}
:ArgumentTypes:	{Integer,
		  Real, Real, Real, Real, Real, 
		  Real, Real, Real, Real, Real,
		  Real, Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mE0ic
:Pattern:	E0i[id_, p1_?c, p2_?c, p3_?c, p4_?c, p5_?c,
		  p1p2_?c, p2p3_?c, p3p4_?c, p4p5_?c, p5p1_?c,
		  m1_?c, m2_?c, m3_?c, m4_?c, m5_?c]
:Arguments:	{Eid[id],
		  N[Re[p1]], N[Im[p1]], N[Re[p2]], N[Im[p2]],
		  N[Re[p3]], N[Im[p3]], N[Re[p4]], N[Im[p4]],
		  N[Re[p5]], N[Im[p5]], N[Re[p1p2]], N[Im[p1p2]], 
		  N[Re[p2p3]], N[Im[p2p3]], N[Re[p3p4]], N[Im[p3p4]],
		  N[Re[p4p5]], N[Im[p4p5]], N[Re[p5p1]], N[Im[p5p1]],
		  N[Re[m1]], N[Im[m1]], N[Re[m2]], N[Im[m2]],
		  N[Re[m3]], N[Im[m3]], N[Re[m4]], N[Im[m4]],
		  N[Re[m5]], N[Im[m5]]}
:ArgumentTypes:	{Integer, Real, Real, Real, Real, Real, 
		  Real, Real, Real, Real, Real, 
		  Real, Real, Real, Real, Real, 
		  Real, Real, Real, Real, Real, 
		  Real, Real, Real, Real, Real, 
		  Real, Real, Real, Real, Real} 
:ReturnType:	Manual
:End:

:Begin:
:Function:	mEget
:Pattern:	Eget[p1_?r, p2_?r, p3_?r, p4_?r, p5_?r,
		  p1p2_?r, p2p3_?r, p3p4_?r, p4p5_?r, p5p1_?r,
		  m1_?r, m2_?r, m3_?r, m4_?r, m5_?r]
:Arguments:	{N[p1], N[p2], N[p3], N[p4], N[p5],
		  N[p1p2], N[p2p3], N[p3p4], N[p4p5], N[p5p1],
		  N[m1], N[m2], N[m3], N[m4], N[m5]}
:ArgumentTypes:	{Real, Real, Real, Real, Real, 
		  Real, Real, Real, Real, Real,
		  Real, Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mEgetc
:Pattern:	Eget[p1_?c, p2_?c, p3_?c, p4_?c, p5_?c,
		  p1p2_?c, p2p3_?c, p3p4_?c,  p4p5_?c, p5p1_?c,
		  m1_?c, m2_?c, m3_?c, m4_?c, m5_?c]
:Arguments:	{N[Re[p1]], N[Im[p1]], N[Re[p2]], N[Im[p2]],
		  N[Re[p3]], N[Im[p3]], N[Re[p4]], N[Im[p4]],
		  N[Re[p5]], N[Im[p5]], N[Re[p1p2]], N[Im[p1p2]], 
		  N[Re[p2p3]], N[Im[p2p3]], N[Re[p3p4]], N[Im[p3p4]],
		  N[Re[p4p5]], N[Im[p4p5]], N[Re[p5p1]], N[Im[p5p1]],
		  N[Re[m1]], N[Im[m1]], N[Re[m2]], N[Im[m2]],
		  N[Re[m3]], N[Im[m3]], N[Re[m4]], N[Im[m4]],
		  N[Re[m5]], N[Im[m5]]}
:ArgumentTypes:	{Real, Real, Real, Real, Real,
		  Real, Real, Real, Real, Real,
		  Real, Real, Real, Real, Real,
		  Real, Real, Real, Real, Real,
		  Real, Real, Real, Real, Real,
		  Real, Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mLi2
:Pattern:	Li2[x_?r]
:Arguments:	{N[x]}
:ArgumentTypes:	{Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mLi2c
:Pattern:	Li2[x_?c]
:Arguments:	{N[Re[x]], N[Im[x]]}
:ArgumentTypes:	{Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mLi2omx
:Pattern:	Li2omx[x_?r]
:Arguments:	{N[x]}
:ArgumentTypes:	{Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mLi2omxc
:Pattern:	Li2omx[x_?c]
:Arguments:	{N[Re[x]], N[Im[x]]}
:ArgumentTypes:	{Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	msetmudim
:Pattern:	SetMudim[mudim_?r]
:Arguments:	{N[mudim]}
:ArgumentTypes:	{Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mgetmudim
:Pattern:	GetMudim[]
:Arguments:	{}
:ArgumentTypes:	{}
:ReturnType:	Real
:End:

:Begin:
:Function:	msetdelta
:Pattern:	SetDelta[delta_?r]
:Arguments:	{N[delta]}
:ArgumentTypes:	{Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mgetdelta
:Pattern:	GetDelta[]
:Arguments:	{}
:ArgumentTypes:	{}
:ReturnType:	Real
:End:

:Begin:
:Function:	msetlambda
:Pattern:	SetLambda[lambda_?r]
:Arguments:	{N[lambda]}
:ArgumentTypes:	{Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mgetlambda
:Pattern:	GetLambda[]
:Arguments:	{}
:ArgumentTypes:	{}
:ReturnType:	Real
:End:

:Begin:
:Function:	msetminmass
:Pattern:	SetMinMass[minmass_?r]
:Arguments:	{N[minmass]}
:ArgumentTypes:	{Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mgetminmass
:Pattern:	GetMinMass[]
:Arguments:	{}
:ArgumentTypes:	{}
:ReturnType:	Real
:End:

:Begin:
:Function:	mclearcache
:Pattern:	ClearCache[]
:Arguments:	{}
:ArgumentTypes:	{}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mmarkcache
:Pattern:	MarkCache[]
:Arguments:	{}
:ArgumentTypes:	{}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mrestorecache
:Pattern:	RestoreCache[]
:Arguments:	{}
:ArgumentTypes:	{}
:ReturnType:	Manual
:End:

:Begin:
:Function:	msetmaxdev
:Pattern:	SetMaxDev[maxdev_?r]
:Arguments:	{N[maxdev]}
:ArgumentTypes:	{Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mgetmaxdev
:Pattern:	GetMaxDev[]
:Arguments:	{}
:ArgumentTypes:	{}
:ReturnType:	Real
:End:

:Begin:
:Function:	msetwarndigits
:Pattern:	SetWarnDigits[warndigits_Integer]
:Arguments:	{warndigits}
:ArgumentTypes:	{Integer}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mgetwarndigits
:Pattern:	GetWarnDigits[]
:Arguments:	{}
:ArgumentTypes:	{}
:ReturnType:	Integer
:End:

:Begin:
:Function:	mseterrdigits
:Pattern:	SetErrDigits[errdigits_Integer]
:Arguments:	{errdigits}
:ArgumentTypes:	{Integer}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mgeterrdigits
:Pattern:	GetErrDigits[]
:Arguments:	{}
:ArgumentTypes:	{}
:ReturnType:	Integer
:End:

:Begin:
:Function:	msetversionkey
:Pattern:	SetVersionKey[versionkey_Integer]
:Arguments:	{versionkey}
:ArgumentTypes:	{Integer}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mgetversionkey
:Pattern:	GetVersionKey[]
:Arguments:	{}
:ArgumentTypes:	{}
:ReturnType:	Integer
:End:

:Begin:
:Function:	msetdebugkey
:Pattern:	SetDebugKey[debugkey_Integer]
:Arguments:	{debugkey}
:ArgumentTypes:	{Integer}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mgetdebugkey
:Pattern:	GetDebugKey[]
:Arguments:	{}
:ArgumentTypes:	{}
:ReturnType:	Integer
:End:

:Begin:
:Function:	msetdebugrange
:Pattern:	SetDebugRange[debugfrom_Integer, debugto_Integer]
:Arguments:	{debugfrom, debugto}
:ArgumentTypes:	{Integer, Integer}
:ReturnType:	Manual
:End:

:Begin:
:Function:	msetcmpbits
:Pattern:	SetCmpBits[cmpbits_Integer]
:Arguments:	{cmpbits}
:ArgumentTypes:	{Integer}
:ReturnType:	Manual
:End:

:Begin:
:Function:	mgetcmpbits
:Pattern:	GetCmpBits[]
:Arguments:	{}
:ArgumentTypes:	{}
:ReturnType:	Integer
:End:

:Evaluate: r = Head[# + 1.] === Real &

:Evaluate: c = Head[# + 1. I] === Complex &

:Evaluate: A0[0] = 0

:Evaluate: MapThread[
	(Derivative[0,1,0,0][B0i][#1, args__] := B0i[#2, args])&,
	{{bb0, bb1, bb00, bb11}, {dbb0, dbb1, dbb00, dbb11}} ]

:Evaluate: PaVe[i__Integer, {p__}, {m__}] :=
	ToExpression[#1 <> "0i"][
	  ToExpression[#2 <> #2 <> ToString/@ Sort[{i}]], p, m ]&[
	  FromCharacterCode[Length[m] + 64],
	  FromCharacterCode[Length[m] + 96] ]

:Evaluate: MapIndexed[(Bid[#1] = #2[[1]] - 1)&, Bids];
	MapIndexed[(Cid[#1] = #2[[1]] - 1)&, Cids];
	MapIndexed[(Did[#1] = #2[[1]] - 1)&, Dids];
	MapIndexed[(Eid[#1] = #2[[1]] - 1)&, Eids];
	Bid[x_] = Cid[x_] = Did[x_] = Eid[x_] = x

:Evaluate: idlist[ids_, x_] := MapThread[Rule,
	  {ids, Apply[Complex, Partition[Chop[x, 10^-14], 2], 1]}];

:Evaluate: End[]

:Evaluate: EndPackage[]


/*
	LoopTools.tm
		provides the LoopTools functions in Mathematica
		this file is part of LoopTools
		last modified 20 Jul 11 th
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <pthread.h>
#include <sched.h>

#include "mathlink.h"
#ifndef MLCONST
#define MLCONST
#endif

#include "clooptools.h"

#ifdef __MINGW32__
#include <io.h>
#include <fcntl.h>
#define pipe(fds) _pipe(fds, 4096, _O_BINARY)
#endif

typedef unsigned char byte;
typedef MLCONST char cchar;
typedef const int cint;
typedef const long clong;

#if QUAD
#define MLPutREAL MLPutReal128
static inline void MLPutREALList(MLINK mlp, CREAL *s, long n)
{
  Real d[n];
  int i;
  for( i = 0; i < n; ++i ) d[i] = ToReal(s[i]);
  MLPutReal128List(mlp, d, n);
}
#else
#define MLPutREAL MLPutReal64
#define MLPutREALList MLPutReal64List
#endif

extern void FORTRAN(fortranflush)();

#define Flush() \
  FORTRAN(fortranflush)(); \
  fflush(stdout)

/******************************************************************/

static int forcestderr = 0;
static int stdoutorig;
static int stdoutpipe[2];
static pthread_t stdouttid;
static int stdoutthr;

static void *MLstdout(void *fd)
{
  static byte *buf = NULL;
  static long size = 0;
  enum { unit = 10240 };
  long len = 0, n = 0;

  do {
    len += n;
    if( size - len < 128 ) buf = realloc(buf, size += unit);
    n = read(*(int *)fd, buf + len, size - len);
  } while( n > 0 );

  if( len ) {
    MLPutFunction(stdlink, "EvaluatePacket", 1);
    MLPutFunction(stdlink, "WriteString", 2);
    MLPutString(stdlink, "stdout");
    MLPutByteString(stdlink, buf, len);
    MLEndPacket(stdlink);

    MLNextPacket(stdlink);
    MLNewPacket(stdlink);
  }

  return NULL;
}

/******************************************************************/

static inline void BeginRedirect()
{
  stdoutthr = forcestderr == 0 &&
    pipe(stdoutpipe) != -1 &&
    pthread_create(&stdouttid, NULL, MLstdout, stdoutpipe) == 0;
  if( !stdoutthr ) stdoutpipe[1] = 2;

  dup2(stdoutpipe[1], 1);
  close(stdoutpipe[1]);
}

/******************************************************************/

static void EndRedirect()
{
  void *ret;

  Flush();
  dup2(stdoutorig, 1);
  if( stdoutthr ) pthread_join(stdouttid, &ret);
}

/******************************************************************/

#define ReturnComplex(expr) \
  Complex result; \
  BeginRedirect(); \
  result = expr; \
  EndRedirect(); \
  MLPutComplex(stdlink, result); \
  MLEndPacket(stdlink)

#define ReturnList(ids, expr, n) \
  COMPLEX *list; \
  BeginRedirect(); \
  list = expr; \
  EndRedirect(); \
  MLPutList(stdlink, ids, list, n); \
  MLEndPacket(stdlink)

#define ReturnVoid() \
  MLPutSymbol(stdlink, "Null"); \
  MLEndPacket(stdlink)

#define _Mr_(v) cReal v
#define _Mri_(v) cReal re_##v, cReal im_##v
#define _Mc_(v) ToComplex2(re_##v, im_##v)

/******************************************************************/

static inline void MLPutComplex(MLINK mlp, cComplex c)
{
  if( Im(c) == 0 ) MLPutREAL(mlp, Re(c));
  else {
    MLPutFunction(mlp, "Complex", 2);
    MLPutREAL(mlp, Re(c));
    MLPutREAL(mlp, Im(c));
  }
}

/******************************************************************/

static inline void MLPutList(MLINK mlp, cchar *ids, COMPLEX *list, cint n)
{
  MLPutFunction(mlp, "LoopTools`Private`idlist", 2);
  MLPutSymbol(mlp, ids);
  MLPutREALList(mlp, (REAL *)list, 2*n);
}

/******************************************************************/

static void mA0(AARGS(_Mr_))
{
  ReturnComplex(A0(AARGS(_Id_)));
}

static void mA0c(AARGS(_Mri_))
{
  ReturnComplex(A0C(AARGS(_Mc_)));
}

/******************************************************************/

static void mA00(AARGS(_Mr_))
{
  ReturnComplex(A00(AARGS(_Id_)));
}

static void mA00c(AARGS(_Mri_))
{
  ReturnComplex(A00C(AARGS(_Mc_)));
}

/******************************************************************/

static void mB0i(cint i, BARGS(_Mr_))
{
  ReturnComplex(B0i(i, BARGS(_Id_)));
}

static void mB0ic(cint i, BARGS(_Mri_))
{
  ReturnComplex(B0iC(i, BARGS(_Mc_)));
}

/******************************************************************/

static void mBget(BARGS(_Mr_))
{
  ReturnList("Bids", Bcache(Bget(BARGS(_Id_))), Nbb);
}

static void mBgetc(BARGS(_Mri_))
{
  ReturnList("Bids", BcacheC(BgetC(BARGS(_Mc_))), Nbb);
}

/******************************************************************/

static void mC0(CARGS(_Mr_))
{
  ReturnComplex(C0(CARGS(_Id_)));
}

static void mC0c(CARGS(_Mri_))
{
  ReturnComplex(C0C(CARGS(_Mc_)));
}

/******************************************************************/

static void mC0i(cint i, CARGS(_Mr_))
{
  ReturnComplex(C0i(i, CARGS(_Id_)));
}

static void mC0ic(cint i, CARGS(_Mri_))
{
  ReturnComplex(C0iC(i, CARGS(_Mc_)));
}

/******************************************************************/

static void mCget(CARGS(_Mr_))
{
  ReturnList("Cids", Ccache(Cget(CARGS(_Id_))), Ncc);
}

static void mCgetc(CARGS(_Mri_))
{
  ReturnList("Cids", CcacheC(CgetC(CARGS(_Mc_))), Ncc);
}

/******************************************************************/

static void mD0(DARGS(_Mr_))
{
  ReturnComplex(D0(DARGS()));
}

static void mD0c(DARGS(_Mri_))
{
  ReturnComplex(D0C(DARGS(_Mc_)));
}

/******************************************************************/

static void mD0i(cint i, DARGS(_Mr_))
{
  ReturnComplex(D0i(i, DARGS(_Id_)));
}

static void mD0ic(cint i, DARGS(_Mri_))
{
  ReturnComplex(D0iC(i, DARGS(_Mc_)));
}

/******************************************************************/

static void mDget(DARGS(_Mr_))
{
  ReturnList("Dids", Dcache(Dget(DARGS(_Id_))), Ndd);
}

static void mDgetc(DARGS(_Mri_))
{
  ReturnList("Dids", DcacheC(DgetC(DARGS(_Mc_))), Ndd);
}

/******************************************************************/

static void mE0(EARGS(_Mr_))
{
  ReturnComplex(E0(EARGS(_Id_)));
}

static void mE0c(EARGS(_Mri_))
{
  ReturnComplex(E0C(EARGS(_Mc_)));
}

/******************************************************************/

static void mE0i(cint i, EARGS(_Mr_))
{
  ReturnComplex(E0i(i, EARGS(_Id_)));
}

static void mE0ic(cint i, EARGS(_Mri_))
{
  ReturnComplex(E0iC(i, EARGS(_Mc_)));
}

/******************************************************************/

static void mEget(EARGS(_Mr_))
{
  ReturnList("Eids", Ecache(Eget(EARGS(_Id_))), Nee);
}

static void mEgetc(EARGS(_Mri_))
{
  ReturnList("Eids", EcacheC(EgetC(EARGS(_Mc_))), Nee);
}

/******************************************************************/

static void mLi2(XARGS(_Mr_))
{
  ReturnComplex(Li2(XARGS(_Id_)));
}

static void mLi2c(XARGS(_Mri_))
{
  ReturnComplex(Li2C(XARGS(_Mc_)));
}

static void mLi2omx(XARGS(_Mr_))
{
  ReturnComplex(Li2omx(XARGS(_Id_)));
}

static void mLi2omxc(XARGS(_Mri_))
{
  ReturnComplex(Li2omxC(XARGS(_Mc_)));
}

/******************************************************************/

static void mclearcache(void)
{
  clearcache();
  ReturnVoid();
}

static void mmarkcache(void)
{
  markcache();
  ReturnVoid();
}

static void mrestorecache(void)
{
  restorecache();
  ReturnVoid();
}

/******************************************************************/

static void msetmudim(cReal mudim)
{
  setmudim(mudim);
  ReturnVoid();
}

static Real mgetmudim(void)
{
  return getmudim();
}

/******************************************************************/

static void msetdelta(cReal delta)
{
  setdelta(delta);
  ReturnVoid();
}

static Real mgetdelta(void)
{
  return getdelta();
}

/******************************************************************/

static void msetlambda(cReal lambda)
{
  setlambda(lambda);
  ReturnVoid();
}

static Real mgetlambda(void)
{
  return getlambda();
}

/******************************************************************/

static void msetminmass(cReal minmass)
{
  setminmass(minmass);
  ReturnVoid();
}

static Real mgetminmass(void)
{
  return getminmass();
}

/******************************************************************/

static void msetmaxdev(cReal maxdev)
{
  setmaxdev(maxdev);
  ReturnVoid();
}

static Real mgetmaxdev(void)
{
  return getmaxdev();
}

/******************************************************************/

static void msetwarndigits(cint warndigits)
{
  setwarndigits(warndigits);
  ReturnVoid();
}

static int mgetwarndigits(void)
{
  return getwarndigits();
}

/******************************************************************/

static void mseterrdigits(cint errdigits)
{
  seterrdigits(errdigits);
  ReturnVoid();
}

static int mgeterrdigits(void)
{
  return geterrdigits();
}

/******************************************************************/

static void msetversionkey(cint versionkey)
{
  setversionkey(versionkey);
  ReturnVoid();
}

static int mgetversionkey(void)
{
  return getversionkey();
}

/******************************************************************/

static void msetdebugkey(cint debugkey)
{
  setdebugkey(debugkey);
  ReturnVoid();
}

static int mgetdebugkey(void)
{
  return getdebugkey();
}

/******************************************************************/

static void msetdebugrange(cint debugfrom, cint debugto)
{
  setdebugrange(debugfrom, debugto);
  ReturnVoid();
}

/******************************************************************/

static void msetcmpbits(cint cmpbits)
{
  setcmpbits(cmpbits);
  ReturnVoid();
}

static int mgetcmpbits(void)
{
  return getcmpbits();
}

/******************************************************************/

int main(int argc, char **argv)
{
  int fd, ret;

	/* make sure a pipe will not overlap with 0, 1, 2 */
  do fd = open("/dev/null", O_WRONLY); while( fd <= 2 );
  close(fd);

  if( getenv("LTFORCESTDERR") ) forcestderr = 1;
  stdoutorig = dup(1);

  dup2(2, 1);
  ltini();
  Flush();

  dup2(stdoutorig, 1);

  ret = MLMain(argc, argv);

  dup2(2, 1);
  ltexi();
  Flush();

  return ret;
}

