(*
	MSSM.mod
		Classes model file for the MSSM
		by Thomas Hahn
		based on the Feynman rules of the MSSM by Arnd Kraft
		last modified 13 Aug 09 by Thomas Hahn

This file contains the definition of the minimal supersymmetric standard
model for FeynArts.  It needs the Generic model file Lorentz.gen.

When you change things, remember:

-- All particles are arranged in classes.  For single particle
   model definitions each particle lives in its own class.

-- For each class the common SelfConjugate behaviour and the
   IndexRange MUST be present in the definitions.

-- IMPORTANT: The coupling matrices MUST be declared in the
   SAME order as the Generic coupling.

This file introduces the following symbols:

	coupling constants and masses:
	------------------------------
	EL:		electron charge (Thomson limit)
	CW, SW:		cosine and sine of weak mixing angle

	MW, MZ:		W, and Z masses
	Mh0, MHH, MA0, MHp: the Higgs masses

	MLE:		lepton class mass
	ME, MM, ML:	lepton masses (e, mu, tau)

	MQU:		u-type quark class mass
	MU, MC, MT:	u-type quark masses (up, charm, top)

	MQD:		d-type quark class mass
	MD, MS, MB:	d-type quark masses (down, strange, bottom)

	MSf:		sfermion mass

	CKM:		quark mixing matrix
			(set CKM = IndexDelta for no quark-mixing)

	CA, SA:		{Cos, Sin}[alpha]
	CB, SB, TB:	{Cos, Sin, Tan}[beta]
	C2A, S2A:	{Cos, Sin}[2 alpha]
	CAB, SAB:	{Cos, Sin}[alpha + beta]
	CBA, SBA:	{Cos, Sin}[beta - alpha]
			where alpha is the (h0, H0) mixing angle
			and tan[beta] is the ratio of the VEVs of
			the two Higgs doublets

	ZNeu:		neutralino mixing matrix (4x4)
	UCha, VCha:	chargino mixing matrices (2x2)
	USf[t]:		t-type sfermion 1-2 mixing matrices (2x2)

	Af[t, i]:	soft breaking parameters
	MUE:		the H1-H2 mixing parameter
*)


(* $HKSign is the sign in the SU(2) covariant derivative,
   i.e. D_\mu = \partial_\mu + $HKSign I g A^a_\mu \tau^a,
   so 1 = Haber-Kane, -1 = Denner conventions *)

If[ !ValueQ[$HKSign], $HKSign = 1 ]

IndexRange[ Index[Generation] ] = Range[3];
IndexRange[ Index[Colour] ] = NoUnfold[Range[3]];
IndexRange[ Index[Sfermion] ] = Range[2];
IndexRange[ Index[Chargino] ] = Range[2];
IndexRange[ Index[Neutralino] ] = Range[4]

IndexStyle[ Index[Generation | Chargino | Neutralino, i_Integer] ] :=
  Alph[i + 8] 

IndexStyle[ Index[Sfermion, i_Integer] ] := Alph[i + 18]

M$ClassesDescription = {

	(* Neutrinos: I_3 = +1/2, Q = 0 *)
  F[1] == {
	SelfConjugate -> False,
	Indices -> {Index[Generation]},
	Mass -> 0,
	PropagatorLabel -> ComposedChar["\\nu", Index[Generation]],
	PropagatorType -> Straight,
	PropagatorArrow -> Forward },

	(* massive Leptons: I_3 = -1/2, Q = -1 *)
  F[2] == {
	SelfConjugate -> False,
	Indices -> {Index[Generation]},
	Mass -> MLE,
	PropagatorLabel -> ComposedChar["e", Index[Generation]],
	PropagatorType -> Straight,
	PropagatorArrow -> Forward },

	(* Quarks (u): I_3 = +1/2, Q = +2/3 *)
  F[3] == {
	SelfConjugate -> False,
	Indices -> {Index[Generation], Index[Colour]},
	Mass -> MQU,
	PropagatorLabel -> ComposedChar["u", Index[Generation]],
	PropagatorType -> Straight,
	PropagatorArrow -> Forward },

	(* Quarks (d): I_3 = -1/2, Q = -1/3 *)
  F[4] == {
	SelfConjugate -> False,
	Indices -> {Index[Generation], Index[Colour]},
	Mass -> MQD,
	PropagatorLabel -> ComposedChar["d", Index[Generation]],
	PropagatorType -> Straight, 
	PropagatorArrow -> Forward },

	(* Neutralinos *)
  F[11] == {
	SelfConjugate -> True,
	Indices -> {Index[Neutralino]},
	Mass -> MNeu,
	PropagatorLabel ->
	  ComposedChar["\\chi", Index[Neutralino], "0", "\\tilde"],
	PropagatorType -> Straight,
	PropagatorArrow -> None },

	(* Charginos *)
  F[12] == {
	SelfConjugate -> False,
	Indices -> {Index[Chargino]},
	Mass -> MCha,
	PropagatorLabel ->
	  ComposedChar["\\chi", Index[Chargino], Null, "\\tilde"],
	PropagatorType -> Straight,
	PropagatorArrow -> Forward },

	(* Gauge bosons: Q = 0 *)
  V[1] == {
	SelfConjugate -> True,
	Indices -> {},
	Mass -> 0,
	PropagatorLabel -> "\\gamma",
	PropagatorType -> Sine,
	PropagatorArrow -> None },

  V[2] == {
	SelfConjugate -> True, 
	Indices -> {},
	Mass -> MZ,
	PropagatorLabel -> "Z",
	PropagatorType -> Sine,
	PropagatorArrow -> None },

	(* Gauge bosons: Q = -1 *)
  V[3] == {
	SelfConjugate -> False,
	Indices -> {},
	Mass -> MW,
	PropagatorLabel -> "W",
	PropagatorType -> Sine,
	PropagatorArrow -> Forward },

	(* CP-even Higgs doublet: Q = 0 *)
  S[1] == {
	SelfConjugate -> True,
	Indices -> {},
	Mass -> Mh0,
	PropagatorLabel -> ComposedChar["h", Null, "0"],
	PropagatorType -> ScalarDash,
	PropagatorArrow -> None },

  S[2] == {
	SelfConjugate -> True,
	Indices -> {},
	Mass -> MHH,
	PropagatorLabel -> ComposedChar["H", Null, "0"],
	PropagatorType -> ScalarDash,
	PropagatorArrow -> None },

	(* CP-odd Higgs doublet: Q = 0 *)
  S[3] == {
	SelfConjugate -> True,
	Indices -> {},
	Mass -> MA0,
	PropagatorLabel -> ComposedChar["A", Null, "0"],
	PropagatorType -> ScalarDash,
	PropagatorArrow -> None },

  S[4] == {
	SelfConjugate -> True,
	Indices -> {},
	Mass -> MZ,
	PropagatorLabel -> ComposedChar["G", Null, "0"],
	PropagatorType -> ScalarDash,
	PropagatorArrow -> None },

	(* charged Higgs doublet: Q = -1 *)
  S[5] == {
	SelfConjugate -> False,
	Indices -> {},
	Mass -> MHp,
	PropagatorLabel -> "H",
	PropagatorType -> ScalarDash,
	PropagatorArrow -> Forward },

  S[6] == {
	SelfConjugate -> False,
	Indices -> {},
	Mass -> MW,
	PropagatorLabel -> "G",
	PropagatorType -> ScalarDash,
	PropagatorArrow -> Forward },

	(* Sneutrinos: Q = 0 *)
  S[11] == {
	SelfConjugate -> False,
	Indices -> {Index[Generation]},
	PropagatorLabel ->
	  ComposedChar["\\nu", Index[Generation], Null, "\\tilde"],
	PropagatorType -> ScalarDash,
	PropagatorArrow -> Forward },

	(* Sleptons: Q = -1 *)
  S[12] == {
	SelfConjugate -> False,
	Indices -> {Index[Sfermion], Index[Generation]},
	PropagatorLabel ->
	  ComposedChar["e", Index[Generation], Index[Sfermion], "\\tilde"],
	PropagatorType -> ScalarDash,
	PropagatorArrow -> Forward },

	(* Squarks (u): Q = +2/3 *)
  S[13] == {
	SelfConjugate -> False,
	Indices -> {Index[Sfermion], Index[Generation], Index[Colour]},
	PropagatorLabel ->
	  ComposedChar["u", Index[Generation], Index[Sfermion], "\\tilde"],
	PropagatorType -> ScalarDash,
	PropagatorArrow -> Forward },

	(* Squarks (d): Q = -1/3 *)
  S[14] == {
	SelfConjugate -> False,
	Indices -> {Index[Sfermion], Index[Generation], Index[Colour]},
	PropagatorLabel ->
	  ComposedChar["d", Index[Generation], Index[Sfermion], "\\tilde"],
	PropagatorType -> ScalarDash,
	PropagatorArrow -> Forward },

	(* Ghosts: Q = 0 *)
  U[1] == {
	SelfConjugate -> False,
	Indices -> {},
	Mass -> 0,
	PropagatorLabel -> ComposedChar["u", "\\gamma"],
	PropagatorType -> GhostDash,
	PropagatorArrow -> Forward },

  U[2] == {
	SelfConjugate -> False,
	Indices -> {},
	Mass -> MZ,
	PropagatorLabel -> ComposedChar["u", "Z"],
	PropagatorType -> GhostDash,
	PropagatorArrow -> Forward },

	(* Ghosts: Q = -1 *)
  U[3] == {
	SelfConjugate -> False,
	Indices -> {},
	Mass -> MW,
	PropagatorLabel -> ComposedChar["u", "-"],
	PropagatorType -> GhostDash,
	PropagatorArrow -> Forward },

  U[4] == {
	SelfConjugate -> False,
	Indices -> {},
	Mass -> MW,
	PropagatorLabel -> ComposedChar["u", "+"],
	PropagatorType -> GhostDash,
	PropagatorArrow -> Forward }
}

MLE[1] = ME;
MLE[2] = MM;
MLE[3] = ML;
MQU[1] = MU;
MQU[2] = MC;
MQU[3] = MT;
MQD[1] = MD;
MQD[2] = MS;
MQD[3] = MB;
MQU[gen_, _] = MQU[gen];
MQD[gen_, _] = MQD[gen]

TheLabel[ F[1, {1}] ] = ComposedChar["\\nu", "e"]; 
TheLabel[ F[1, {2}] ] = ComposedChar["\\nu", "\\mu"]; 
TheLabel[ F[1, {3}] ] = ComposedChar["\\nu", "\\tau"]; 
TheLabel[ F[2, {1}] ] = "e"; 
TheLabel[ F[2, {2}] ] = "\\mu"; 
TheLabel[ F[2, {3}] ] = "\\tau";
TheLabel[ F[3, {1, ___}] ] = "u"; 
TheLabel[ F[3, {2, ___}] ] = "c";
TheLabel[ F[3, {3, ___}] ] = "t";
TheLabel[ F[4, {1, ___}] ] = "d"; 
TheLabel[ F[4, {2, ___}] ] = "s";
TheLabel[ F[4, {3, ___}] ] = "b"

TheMass[ S[11, {gen_, ___}] ] = MSf[1, 1, gen];
TheMass[ S[typ:12 | 13 | 14, {sf_, gen_, ___}] ] := MSf[sf, typ - 10, gen]

TheLabel[ S[11, {1}] ] = ComposedChar["\\nu", "e", Null, "\\tilde"];
TheLabel[ S[11, {2}] ] = ComposedChar["\\nu", "\\mu", Null, "\\tilde"];
TheLabel[ S[11, {3}] ] = ComposedChar["\\nu", "\\tau", Null, "\\tilde"];
TheLabel[ S[12, {sf_, 1}] ] :=
  ComposedChar["e", Null, IndexStyle[sf], "\\tilde"];
TheLabel[ S[12, {sf_, 2}] ] :=
  ComposedChar["\\mu", Null, IndexStyle[sf], "\\tilde"];
TheLabel[ S[12, {sf_, 3}] ] :=
  ComposedChar["\\tau", Null, IndexStyle[sf], "\\tilde"];
TheLabel[ S[13, {sf_, 1, ___}] ] :=
  ComposedChar["u", Null, IndexStyle[sf], "\\tilde"];
TheLabel[ S[13, {sf_, 2, ___}] ] :=
  ComposedChar["c", Null, IndexStyle[sf], "\\tilde"];
TheLabel[ S[13, {sf_, 3, ___}] ] :=
  ComposedChar["t", Null, IndexStyle[sf], "\\tilde"];
TheLabel[ S[14, {sf_, 1, ___}] ] :=
  ComposedChar["d", Null, IndexStyle[sf], "\\tilde"];
TheLabel[ S[14, {sf_, 2, ___}] ] :=
  ComposedChar["s", Null, IndexStyle[sf], "\\tilde"];
TheLabel[ S[14, {sf_, 3, ___}] ] :=
  ComposedChar["b", Null, IndexStyle[sf], "\\tilde"]

GaugeXi[ V[1] ] = GaugeXi[A];
GaugeXi[ V[2] ] = GaugeXi[Z];
GaugeXi[ V[3] ] = GaugeXi[W];
GaugeXi[ U[1] ] = GaugeXi[A];
GaugeXi[ U[2] ] = GaugeXi[Z];
GaugeXi[ U[3] ] = GaugeXi[W];
GaugeXi[ U[4] ] = GaugeXi[W];
GaugeXi[ S[4] ] = GaugeXi[Z];
GaugeXi[ S[6] ] = GaugeXi[W];
GaugeXi[ S[_Integer, ___] ] = 1

M$LastModelRules = {}


(* some short-hands for excluding classes of particles *)

NoGeneration1 = ExcludeParticles ->
  {F[1|2|3|4, {1, ___}], S[11, {1, ___}], S[12|13|14, {_, 1, ___}]}

NoGeneration2 = ExcludeParticles ->
  {F[1|2|3|4, {2, ___}], S[11, {2, ___}], S[12|13|14, {_, 2, ___}]}

NoGeneration3 = ExcludeParticles ->
  {F[1|2|3|4, {3, ___}], S[11, {3, ___}], S[12|13|14, {_, 3, ___}]}

NoSUSYParticles = ExcludeParticles ->
  {S[11], S[12], S[13], S[14], S[2], S[3], S[5], F[11], F[12]}

THDMParticles = ExcludeParticles ->
  {S[11], S[12], S[13], S[14], F[11], F[12]}

NoElectronHCoupling =
  ExcludeFieldPoints -> {
    FieldPoint[0][-F[2, {1}], F[2, {1}], S],
    FieldPoint[0][-F[2, {1}], F[1, {1}], S] }

NoLightFHCoupling =
  ExcludeFieldPoints -> {
    FieldPoint[_][-F[2], F[2], S],
    FieldPoint[_][-F[2], F[1], S],
    FieldPoint[_][-F[3, {1, ___}], F[3, {1, ___}], S],
    FieldPoint[_][-F[3, {2, ___}], F[3, {2, ___}], S],
    FieldPoint[_][-F[4], F[4], S],
    FieldPoint[_][-F[4], F[3, {1, ___}], S],
    FieldPoint[_][-F[4], F[3, {2, ___}], S] }

M$CouplingMatrices = {C[S[6], -S[6], V[1]] == {{I*EL}}, C[S[6], -S[6], V[2]] == 
  {{((I/2)*EL*(CW^2 - SW^2)*$HKSign)/(CW*SW)}}, 
 C[S[4], S[6], -V[3]] == {{(EL*$HKSign)/(2*SW)}}, 
 C[S[4], -S[6], V[3]] == {{(EL*$HKSign)/(2*SW)}}, 
 C[S[6], V[1], -V[3]] == {{I*EL*MW*$HKSign}}, 
 C[-S[6], V[1], V[3]] == {{I*EL*MW*$HKSign}}, 
 C[S[6], V[2], -V[3]] == {{((-I)*EL*MW*SW)/CW}}, 
 C[-S[6], V[2], V[3]] == {{((-I)*EL*MW*SW)/CW}}, 
 C[V[1], -V[3], V[3]] == {{(-I)*EL}}, C[V[2], -V[3], V[3]] == 
  {{((-I)*CW*EL*$HKSign)/SW}}, C[S[4], U[3], -U[3]] == 
  {{-(EL*MW*GaugeXi[W])/(2*SW)}}, C[S[4], U[4], -U[4]] == 
  {{(EL*MW*GaugeXi[W])/(2*SW)}}, C[S[6], U[1], -U[3]] == 
  {{(-I)*EL*MW*$HKSign*GaugeXi[W]}}, C[-S[6], U[1], -U[4]] == 
  {{(-I)*EL*MW*$HKSign*GaugeXi[W]}}, C[S[6], U[2], -U[3]] == 
  {{((I/2)*EL*MW*(-CW^2 + SW^2)*GaugeXi[W])/(CW*SW)}}, 
 C[-S[6], U[2], -U[4]] == {{((I/2)*EL*MW*(-CW^2 + SW^2)*GaugeXi[W])/
     (CW*SW)}}, C[S[6], U[4], -U[2]] == {{((I/2)*EL*MW*GaugeXi[Z])/(CW*SW)}}, 
 C[-S[6], U[3], -U[2]] == {{((I/2)*EL*MW*GaugeXi[Z])/(CW*SW)}}, 
 C[-U[3], U[3], V[1]] == {{(-I)*EL}, {0}}, C[-U[4], U[4], V[1]] == 
  {{I*EL}, {0}}, C[-U[3], U[3], V[2]] == {{((-I)*CW*EL*$HKSign)/SW}, {0}}, 
 C[-U[4], U[4], V[2]] == {{(I*CW*EL*$HKSign)/SW}, {0}}, 
 C[-U[3], U[1], V[3]] == {{I*EL}, {0}}, C[-U[4], U[1], -V[3]] == 
  {{(-I)*EL}, {0}}, C[-U[1], U[4], V[3]] == {{(-I)*EL}, {0}}, 
 C[-U[1], U[3], -V[3]] == {{I*EL}, {0}}, C[-U[3], U[2], V[3]] == 
  {{(I*CW*EL*$HKSign)/SW}, {0}}, C[-U[4], U[2], -V[3]] == 
  {{((-I)*CW*EL*$HKSign)/SW}, {0}}, C[-U[2], U[4], V[3]] == 
  {{((-I)*CW*EL*$HKSign)/SW}, {0}}, C[-U[2], U[3], -V[3]] == 
  {{(I*CW*EL*$HKSign)/SW}, {0}}, C[S[1], S[1], V[2], V[2]] == 
  {{((I/2)*EL^2)/(CW^2*SW^2)}}, C[S[1], S[1], V[3], -V[3]] == 
  {{((I/2)*EL^2)/SW^2}}, C[S[4], S[4], V[2], V[2]] == 
  {{((I/2)*EL^2)/(CW^2*SW^2)}}, C[S[4], S[4], V[3], -V[3]] == 
  {{((I/2)*EL^2)/SW^2}}, C[S[6], -S[6], V[1], V[1]] == {{(2*I)*EL^2}}, 
 C[S[6], -S[6], V[1], V[2]] == {{(I*EL^2*(CW^2 - SW^2)*$HKSign)/(CW*SW)}}, 
 C[S[6], -S[6], V[2], V[2]] == {{((I/2)*EL^2*(CW^2 - SW^2)^2)/(CW^2*SW^2)}}, 
 C[S[6], -S[6], V[3], -V[3]] == {{((I/2)*EL^2)/SW^2}}, 
 C[V[1], V[1], V[3], -V[3]] == {{(-2*I)*EL^2}, {I*EL^2}, {I*EL^2}}, 
 C[V[1], V[2], V[3], -V[3]] == {{((-2*I)*CW*EL^2*$HKSign)/SW}, 
   {(I*CW*EL^2*$HKSign)/SW}, {(I*CW*EL^2*$HKSign)/SW}}, 
 C[V[2], V[2], V[3], -V[3]] == {{((-2*I)*CW^2*EL^2)/SW^2}, 
   {(I*CW^2*EL^2)/SW^2}, {(I*CW^2*EL^2)/SW^2}}, 
 C[V[3], V[3], -V[3], -V[3]] == {{((2*I)*EL^2)/SW^2}, {((-I)*EL^2)/SW^2}, 
   {((-I)*EL^2)/SW^2}}, C[S[1], S[1], S[1]] == 
  {{(((-3*I)/2)*C2A*EL*MW*SAB)/(CW^2*SW)}}, C[S[1], S[1], S[2]] == 
  {{((I/2)*EL*MW*(C2A*CAB - 2*S2A*SAB))/(CW^2*SW)}}, 
 C[S[1], S[2], S[2]] == {{((I/2)*EL*MW*(2*CAB*S2A + C2A*SAB))/(CW^2*SW)}}, 
 C[S[2], S[2], S[2]] == {{(((-3*I)/2)*C2A*CAB*EL*MW)/(CW^2*SW)}}, 
 C[S[1], S[3], S[3]] == {{((-I/2)*C2B*EL*MW*SAB)/(CW^2*SW)}}, 
 C[S[1], S[3], S[4]] == {{((-I/2)*EL*MW*S2B*SAB)/(CW^2*SW)}}, 
 C[S[1], S[4], S[4]] == {{((I/2)*C2B*EL*MW*SAB)/(CW^2*SW)}}, 
 C[S[2], S[3], S[3]] == {{((I/2)*C2B*CAB*EL*MW)/(CW^2*SW)}}, 
 C[S[2], S[3], S[4]] == {{((I/2)*CAB*EL*MW*S2B)/(CW^2*SW)}}, 
 C[S[2], S[4], S[4]] == {{((-I/2)*C2B*CAB*EL*MW)/(CW^2*SW)}}, 
 C[S[1], S[5], -S[5]] == {{((-I)*EL*MW*((C2B*SAB)/(2*CW^2) + SBA))/SW}}, 
 C[S[1], S[5], -S[6]] == {{((I/2)*EL*MW*(CBA - (S2B*SAB)/CW^2))/SW}}, 
 C[S[1], S[6], -S[5]] == {{((I/2)*EL*MW*(CBA - (S2B*SAB)/CW^2))/SW}}, 
 C[S[1], S[6], -S[6]] == {{((I/2)*C2B*EL*MW*SAB)/(CW^2*SW)}}, 
 C[S[2], S[5], -S[5]] == {{((-I)*(CBA - (C2B*CAB)/(2*CW^2))*EL*MW)/SW}}, 
 C[S[2], S[5], -S[6]] == {{((-I/2)*EL*MW*(-((CAB*S2B)/CW^2) + SBA))/SW}}, 
 C[S[2], S[6], -S[5]] == {{((-I/2)*EL*MW*(-((CAB*S2B)/CW^2) + SBA))/SW}}, 
 C[S[2], S[6], -S[6]] == {{((-I/2)*C2B*CAB*EL*MW)/(CW^2*SW)}}, 
 C[S[3], S[5], -S[6]] == {{-(EL*MW)/(2*SW)}}, 
 C[S[3], S[6], -S[5]] == {{(EL*MW)/(2*SW)}}, 
 C[S[1], S[3], V[2]] == {{(CBA*EL*$HKSign)/(2*CW*SW)}}, 
 C[S[1], S[4], V[2]] == {{(EL*SBA*$HKSign)/(2*CW*SW)}}, 
 C[S[2], S[3], V[2]] == {{-(EL*SBA*$HKSign)/(2*CW*SW)}}, 
 C[S[2], S[4], V[2]] == {{(CBA*EL*$HKSign)/(2*CW*SW)}}, 
 C[S[5], -S[5], V[1]] == {{I*EL}}, C[S[5], -S[5], V[2]] == 
  {{((I/2)*EL*(CW^2 - SW^2)*$HKSign)/(CW*SW)}}, 
 C[S[1], S[5], -V[3]] == {{((-I/2)*CBA*EL*$HKSign)/SW}}, 
 C[S[1], S[6], -V[3]] == {{((-I/2)*EL*SBA*$HKSign)/SW}}, 
 C[S[2], S[5], -V[3]] == {{((I/2)*EL*SBA*$HKSign)/SW}}, 
 C[S[2], S[6], -V[3]] == {{((-I/2)*CBA*EL*$HKSign)/SW}}, 
 C[S[1], -S[5], V[3]] == {{((I/2)*CBA*EL*$HKSign)/SW}}, 
 C[S[1], -S[6], V[3]] == {{((I/2)*EL*SBA*$HKSign)/SW}}, 
 C[S[2], -S[5], V[3]] == {{((-I/2)*EL*SBA*$HKSign)/SW}}, 
 C[S[2], -S[6], V[3]] == {{((I/2)*CBA*EL*$HKSign)/SW}}, 
 C[S[3], S[5], -V[3]] == {{(EL*$HKSign)/(2*SW)}}, 
 C[S[3], -S[5], V[3]] == {{(EL*$HKSign)/(2*SW)}}, 
 C[S[1], V[2], V[2]] == {{(I*EL*MW*SBA)/(CW^2*SW)}}, 
 C[S[2], V[2], V[2]] == {{(I*CBA*EL*MW)/(CW^2*SW)}}, 
 C[S[1], V[3], -V[3]] == {{(I*EL*MW*SBA)/SW}}, 
 C[S[2], V[3], -V[3]] == {{(I*CBA*EL*MW)/SW}}, 
 C[S[1], U[2], -U[2]] == {{((-I/2)*EL*MW*SBA*GaugeXi[Z])/(CW^2*SW)}}, 
 C[S[2], U[2], -U[2]] == {{((-I/2)*CBA*EL*MW*GaugeXi[Z])/(CW^2*SW)}}, 
 C[S[1], U[3], -U[3]] == {{((-I/2)*EL*MW*SBA*GaugeXi[W])/SW}}, 
 C[S[2], U[3], -U[3]] == {{((-I/2)*CBA*EL*MW*GaugeXi[W])/SW}}, 
 C[S[1], U[4], -U[4]] == {{((-I/2)*EL*MW*SBA*GaugeXi[W])/SW}}, 
 C[S[2], U[4], -U[4]] == {{((-I/2)*CBA*EL*MW*GaugeXi[W])/SW}}, 
 C[S[1], S[1], S[1], S[1]] == {{(((-3*I)/4)*C2A^2*EL^2)/(CW^2*SW^2)}}, 
 C[S[1], S[1], S[1], S[2]] == {{(((-3*I)/4)*C2A*EL^2*S2A)/(CW^2*SW^2)}}, 
 C[S[1], S[1], S[2], S[2]] == {{((-I/4)*EL^2*(-1 + 3*S2A^2))/(CW^2*SW^2)}}, 
 C[S[1], S[2], S[2], S[2]] == {{(((3*I)/4)*C2A*EL^2*S2A)/(CW^2*SW^2)}}, 
 C[S[2], S[2], S[2], S[2]] == {{(((-3*I)/4)*C2A^2*EL^2)/(CW^2*SW^2)}}, 
 C[S[1], S[1], S[3], S[3]] == {{((-I/4)*C2A*C2B*EL^2)/(CW^2*SW^2)}}, 
 C[S[1], S[1], S[3], S[4]] == {{((-I/4)*C2A*EL^2*S2B)/(CW^2*SW^2)}}, 
 C[S[1], S[1], S[4], S[4]] == {{((I/4)*C2A*C2B*EL^2)/(CW^2*SW^2)}}, 
 C[S[1], S[2], S[3], S[3]] == {{((-I/4)*C2B*EL^2*S2A)/(CW^2*SW^2)}}, 
 C[S[1], S[2], S[3], S[4]] == {{((-I/4)*EL^2*S2A*S2B)/(CW^2*SW^2)}}, 
 C[S[1], S[2], S[4], S[4]] == {{((I/4)*C2B*EL^2*S2A)/(CW^2*SW^2)}}, 
 C[S[2], S[2], S[3], S[3]] == {{((I/4)*C2A*C2B*EL^2)/(CW^2*SW^2)}}, 
 C[S[2], S[2], S[3], S[4]] == {{((I/4)*C2A*EL^2*S2B)/(CW^2*SW^2)}}, 
 C[S[2], S[2], S[4], S[4]] == {{((-I/4)*C2A*C2B*EL^2)/(CW^2*SW^2)}}, 
 C[S[1], S[1], S[5], -S[5]] == 
  {{((-I/4)*EL^2*(1 - S2A*S2B + (C2A*C2B*SW^2)/CW^2))/SW^2}}, 
 C[S[1], S[1], S[5], -S[6]] == 
  {{((-I/4)*EL^2*(C2B*S2A + (C2A*S2B*SW^2)/CW^2))/SW^2}}, 
 C[S[1], S[1], S[6], -S[5]] == 
  {{((-I/4)*EL^2*(C2B*S2A + (C2A*S2B*SW^2)/CW^2))/SW^2}}, 
 C[S[1], S[1], S[6], -S[6]] == 
  {{((-I/4)*EL^2*(1 + S2A*S2B - (C2A*C2B*SW^2)/CW^2))/SW^2}}, 
 C[S[1], S[2], S[5], -S[5]] == 
  {{((-I/4)*EL^2*(C2A*S2B + (C2B*S2A*SW^2)/CW^2))/SW^2}}, 
 C[S[1], S[2], S[5], -S[6]] == 
  {{((I/4)*EL^2*(C2A*C2B - (S2A*S2B*SW^2)/CW^2))/SW^2}}, 
 C[S[1], S[2], S[6], -S[5]] == 
  {{((I/4)*EL^2*(C2A*C2B - (S2A*S2B*SW^2)/CW^2))/SW^2}}, 
 C[S[1], S[2], S[6], -S[6]] == 
  {{((I/4)*EL^2*(C2A*S2B + (C2B*S2A*SW^2)/CW^2))/SW^2}}, 
 C[S[2], S[2], S[5], -S[5]] == 
  {{((-I/4)*EL^2*(1 + S2A*S2B - (C2A*C2B*SW^2)/CW^2))/SW^2}}, 
 C[S[2], S[2], S[5], -S[6]] == 
  {{((I/4)*EL^2*(C2B*S2A + (C2A*S2B*SW^2)/CW^2))/SW^2}}, 
 C[S[2], S[2], S[6], -S[5]] == 
  {{((I/4)*EL^2*(C2B*S2A + (C2A*S2B*SW^2)/CW^2))/SW^2}}, 
 C[S[2], S[2], S[6], -S[6]] == 
  {{((-I/4)*EL^2*(1 - S2A*S2B + (C2A*C2B*SW^2)/CW^2))/SW^2}}, 
 C[S[1], S[3], S[5], -S[6]] == {{-(EL^2*SBA)/(4*SW^2)}}, 
 C[S[1], S[3], S[6], -S[5]] == {{(EL^2*SBA)/(4*SW^2)}}, 
 C[S[1], S[4], S[5], -S[6]] == {{(CBA*EL^2)/(4*SW^2)}}, 
 C[S[1], S[4], S[6], -S[5]] == {{-(CBA*EL^2)/(4*SW^2)}}, 
 C[S[2], S[3], S[5], -S[6]] == {{-(CBA*EL^2)/(4*SW^2)}}, 
 C[S[2], S[3], S[6], -S[5]] == {{(CBA*EL^2)/(4*SW^2)}}, 
 C[S[2], S[4], S[5], -S[6]] == {{-(EL^2*SBA)/(4*SW^2)}}, 
 C[S[2], S[4], S[6], -S[5]] == {{(EL^2*SBA)/(4*SW^2)}}, 
 C[S[3], S[3], S[3], S[3]] == {{(((-3*I)/4)*C2B^2*EL^2)/(CW^2*SW^2)}}, 
 C[S[3], S[3], S[3], S[4]] == {{(((-3*I)/4)*C2B*EL^2*S2B)/(CW^2*SW^2)}}, 
 C[S[3], S[3], S[4], S[4]] == {{((-I/4)*EL^2*(-1 + 3*S2B^2))/(CW^2*SW^2)}}, 
 C[S[3], S[4], S[4], S[4]] == {{(((3*I)/4)*C2B*EL^2*S2B)/(CW^2*SW^2)}}, 
 C[S[4], S[4], S[4], S[4]] == {{(((-3*I)/4)*C2B^2*EL^2)/(CW^2*SW^2)}}, 
 C[S[3], S[3], S[5], -S[5]] == {{((-I/4)*C2B^2*EL^2)/(CW^2*SW^2)}}, 
 C[S[3], S[3], S[5], -S[6]] == {{((-I/4)*C2B*EL^2*S2B)/(CW^2*SW^2)}}, 
 C[S[3], S[3], S[6], -S[5]] == {{((-I/4)*C2B*EL^2*S2B)/(CW^2*SW^2)}}, 
 C[S[3], S[3], S[6], -S[6]] == 
  {{((-I/4)*EL^2*(1 + S2B^2 - (C2B^2*SW^2)/CW^2))/SW^2}}, 
 C[S[3], S[4], S[5], -S[5]] == {{((-I/4)*C2B*EL^2*S2B)/(CW^2*SW^2)}}, 
 C[S[3], S[4], S[5], -S[6]] == 
  {{((I/4)*EL^2*(C2B^2 - (S2B^2*SW^2)/CW^2))/SW^2}}, 
 C[S[3], S[4], S[6], -S[5]] == 
  {{((I/4)*EL^2*(C2B^2 - (S2B^2*SW^2)/CW^2))/SW^2}}, 
 C[S[3], S[4], S[6], -S[6]] == {{((I/4)*C2B*EL^2*S2B)/(CW^2*SW^2)}}, 
 C[S[4], S[4], S[5], -S[5]] == 
  {{((-I/4)*EL^2*(1 + S2B^2 - (C2B^2*SW^2)/CW^2))/SW^2}}, 
 C[S[4], S[4], S[5], -S[6]] == {{((I/4)*C2B*EL^2*S2B)/(CW^2*SW^2)}}, 
 C[S[4], S[4], S[6], -S[5]] == {{((I/4)*C2B*EL^2*S2B)/(CW^2*SW^2)}}, 
 C[S[4], S[4], S[6], -S[6]] == {{((-I/4)*C2B^2*EL^2)/(CW^2*SW^2)}}, 
 C[S[5], S[5], -S[5], -S[5]] == {{((-I/2)*C2B^2*EL^2)/(CW^2*SW^2)}}, 
 C[S[5], S[5], -S[5], -S[6]] == {{((-I/2)*C2B*EL^2*S2B)/(CW^2*SW^2)}}, 
 C[S[5], S[5], -S[6], -S[6]] == {{((-I/2)*EL^2*S2B^2)/(CW^2*SW^2)}}, 
 C[S[5], S[6], -S[5], -S[5]] == {{((-I/2)*C2B*EL^2*S2B)/(CW^2*SW^2)}}, 
 C[S[5], S[6], -S[5], -S[6]] == {{((I/4)*EL^2*(C2B^2 - S2B^2))/(CW^2*SW^2)}}, 
 C[S[5], S[6], -S[6], -S[6]] == {{((I/2)*C2B*EL^2*S2B)/(CW^2*SW^2)}}, 
 C[S[6], S[6], -S[5], -S[5]] == {{((-I/2)*EL^2*S2B^2)/(CW^2*SW^2)}}, 
 C[S[6], S[6], -S[5], -S[6]] == {{((I/2)*C2B*EL^2*S2B)/(CW^2*SW^2)}}, 
 C[S[6], S[6], -S[6], -S[6]] == {{((-I/2)*C2B^2*EL^2)/(CW^2*SW^2)}}, 
 C[S[1], S[5], V[1], -V[3]] == {{((I/2)*CBA*EL^2*$HKSign)/SW}}, 
 C[S[1], S[5], V[2], -V[3]] == {{((-I/2)*CBA*EL^2)/CW}}, 
 C[S[1], S[6], V[1], -V[3]] == {{((I/2)*EL^2*SBA*$HKSign)/SW}}, 
 C[S[1], S[6], V[2], -V[3]] == {{((-I/2)*EL^2*SBA)/CW}}, 
 C[S[1], -S[5], V[1], V[3]] == {{((I/2)*CBA*EL^2*$HKSign)/SW}}, 
 C[S[1], -S[5], V[2], V[3]] == {{((-I/2)*CBA*EL^2)/CW}}, 
 C[S[1], -S[6], V[1], V[3]] == {{((I/2)*EL^2*SBA*$HKSign)/SW}}, 
 C[S[1], -S[6], V[2], V[3]] == {{((-I/2)*EL^2*SBA)/CW}}, 
 C[S[2], S[2], V[2], V[2]] == {{((I/2)*EL^2)/(CW^2*SW^2)}}, 
 C[S[2], S[2], V[3], -V[3]] == {{((I/2)*EL^2)/SW^2}}, 
 C[S[2], S[5], V[1], -V[3]] == {{((-I/2)*EL^2*SBA*$HKSign)/SW}}, 
 C[S[2], S[5], V[2], -V[3]] == {{((I/2)*EL^2*SBA)/CW}}, 
 C[S[2], S[6], V[1], -V[3]] == {{((I/2)*CBA*EL^2*$HKSign)/SW}}, 
 C[S[2], S[6], V[2], -V[3]] == {{((-I/2)*CBA*EL^2)/CW}}, 
 C[S[2], -S[5], V[1], V[3]] == {{((-I/2)*EL^2*SBA*$HKSign)/SW}}, 
 C[S[2], -S[5], V[2], V[3]] == {{((I/2)*EL^2*SBA)/CW}}, 
 C[S[2], -S[6], V[1], V[3]] == {{((I/2)*CBA*EL^2*$HKSign)/SW}}, 
 C[S[2], -S[6], V[2], V[3]] == {{((-I/2)*CBA*EL^2)/CW}}, 
 C[S[3], S[3], V[2], V[2]] == {{((I/2)*EL^2)/(CW^2*SW^2)}}, 
 C[S[3], S[3], V[3], -V[3]] == {{((I/2)*EL^2)/SW^2}}, 
 C[S[3], S[5], V[1], -V[3]] == {{-(EL^2*$HKSign)/(2*SW)}}, 
 C[S[3], S[5], V[2], -V[3]] == {{EL^2/(2*CW)}}, 
 C[S[3], -S[5], V[1], V[3]] == {{(EL^2*$HKSign)/(2*SW)}}, 
 C[S[3], -S[5], V[2], V[3]] == {{-EL^2/(2*CW)}}, 
 C[S[4], S[6], V[1], -V[3]] == {{-(EL^2*$HKSign)/(2*SW)}}, 
 C[S[4], S[6], V[2], -V[3]] == {{EL^2/(2*CW)}}, 
 C[S[4], -S[6], V[1], V[3]] == {{(EL^2*$HKSign)/(2*SW)}}, 
 C[S[4], -S[6], V[2], V[3]] == {{-EL^2/(2*CW)}}, 
 C[S[5], -S[5], V[1], V[1]] == {{(2*I)*EL^2}}, 
 C[S[5], -S[5], V[1], V[2]] == {{(I*EL^2*(CW^2 - SW^2)*$HKSign)/(CW*SW)}}, 
 C[S[5], -S[5], V[2], V[2]] == {{((I/2)*EL^2*(CW^2 - SW^2)^2)/(CW^2*SW^2)}}, 
 C[S[5], -S[5], V[3], -V[3]] == {{((I/2)*EL^2)/SW^2}}, 
 C[F[2, {j1}], -F[2, {j2}], S[1]] == 
  {{((I/2)*EL*SA*IndexDelta[j1, j2]*Mass[F[2, {j1}]])/(CB*MW*SW)}, 
   {((I/2)*EL*SA*IndexDelta[j1, j2]*Mass[F[2, {j1}]])/(CB*MW*SW)}}, 
 C[F[3, {j1, o1}], -F[3, {j2, o2}], S[1]] == 
  {{((-I/2)*CA*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[3, {j1, o1}]])/
     (MW*SB*SW)}, {((-I/2)*CA*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      Mass[F[3, {j1, o1}]])/(MW*SB*SW)}}, 
 C[F[4, {j1, o1}], -F[4, {j2, o2}], S[1]] == 
  {{((I/2)*EL*SA*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[4, {j1, o1}]])/
     (CB*MW*SW)}, {((I/2)*EL*SA*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      Mass[F[4, {j1, o1}]])/(CB*MW*SW)}}, C[F[2, {j1}], -F[2, {j2}], S[4]] == 
  {{-(EL*IndexDelta[j1, j2]*Mass[F[2, {j1}]])/(2*MW*SW)}, 
   {(EL*IndexDelta[j1, j2]*Mass[F[2, {j1}]])/(2*MW*SW)}}, 
 C[F[3, {j1, o1}], -F[3, {j2, o2}], S[4]] == 
  {{(EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[3, {j1, o1}]])/
     (2*MW*SW)}, {-(EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
       Mass[F[3, {j1, o1}]])/(2*MW*SW)}}, 
 C[F[4, {j1, o1}], -F[4, {j2, o2}], S[4]] == 
  {{-(EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[4, {j1, o1}]])/
     (2*MW*SW)}, {(EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      Mass[F[4, {j1, o1}]])/(2*MW*SW)}}, C[-F[2, {j2}], F[2, {j1}], V[1]] == 
  {{I*EL*IndexDelta[j1, j2]}, {I*EL*IndexDelta[j1, j2]}}, 
 C[-F[3, {j2, o1}], F[3, {j1, o2}], V[1]] == 
  {{((-2*I)/3)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]}, 
   {((-2*I)/3)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]}}, 
 C[-F[4, {j2, o1}], F[4, {j1, o2}], V[1]] == 
  {{(I/3)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]}, 
   {(I/3)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]}}, 
 C[-F[1, {j2}], F[1, {j1}], V[2]] == 
  {{((-I/2)*EL*$HKSign*IndexDelta[j1, j2])/(CW*SW)}, {0}}, 
 C[-F[2, {j2}], F[2, {j1}], V[2]] == 
  {{((-I)*EL*(-1/2 + SW^2)*$HKSign*IndexDelta[j1, j2])/(CW*SW)}, 
   {((-I)*EL*SW*$HKSign*IndexDelta[j1, j2])/CW}}, 
 C[-F[3, {j2, o1}], F[3, {j1, o2}], V[2]] == 
  {{((I/6)*EL*(-3 + 4*SW^2)*$HKSign*IndexDelta[j1, j2]*IndexDelta[o1, o2])/
     (CW*SW)}, {(((2*I)/3)*EL*SW*$HKSign*IndexDelta[j1, j2]*
      IndexDelta[o1, o2])/CW}}, C[-F[4, {j2, o1}], F[4, {j1, o2}], V[2]] == 
  {{((-I/6)*EL*(-3 + 2*SW^2)*$HKSign*IndexDelta[j1, j2]*IndexDelta[o1, o2])/
     (CW*SW)}, {((-I/3)*EL*SW*$HKSign*IndexDelta[j1, j2]*IndexDelta[o1, o2])/
     CW}}, C[F[2, {j1}], -F[2, {j2}], S[2]] == 
  {{((-I/2)*CA*EL*IndexDelta[j1, j2]*Mass[F[2, {j1}]])/(CB*MW*SW)}, 
   {((-I/2)*CA*EL*IndexDelta[j1, j2]*Mass[F[2, {j1}]])/(CB*MW*SW)}}, 
 C[F[3, {j1, o1}], -F[3, {j2, o2}], S[2]] == 
  {{((-I/2)*EL*SA*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[3, {j1, o1}]])/
     (MW*SB*SW)}, {((-I/2)*EL*SA*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      Mass[F[3, {j1, o1}]])/(MW*SB*SW)}}, 
 C[F[4, {j1, o1}], -F[4, {j2, o2}], S[2]] == 
  {{((-I/2)*CA*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[4, {j1, o1}]])/
     (CB*MW*SW)}, {((-I/2)*CA*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      Mass[F[4, {j1, o1}]])/(CB*MW*SW)}}, C[F[2, {j1}], -F[2, {j2}], S[3]] == 
  {{(EL*TB*IndexDelta[j1, j2]*Mass[F[2, {j1}]])/(2*MW*SW)}, 
   {-(EL*TB*IndexDelta[j1, j2]*Mass[F[2, {j1}]])/(2*MW*SW)}}, 
 C[F[3, {j1, o1}], -F[3, {j2, o2}], S[3]] == 
  {{(EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[3, {j1, o1}]])/
     (2*MW*SW*TB)}, {-(EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
       Mass[F[3, {j1, o1}]])/(2*MW*SW*TB)}}, 
 C[F[4, {j1, o1}], -F[4, {j2, o2}], S[3]] == 
  {{(EL*TB*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[4, {j1, o1}]])/
     (2*MW*SW)}, {-(EL*TB*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
       Mass[F[4, {j1, o1}]])/(2*MW*SW)}}, C[F[1, {j1}], -F[2, {j2}], S[6]] == 
  {{((-I)*EL*IndexDelta[j1, j2]*Mass[F[2, {j2}]])/(Sqrt[2]*MW*SW)}, {0}}, 
 C[F[2, {j1}], -F[1, {j2}], -S[6]] == 
  {{0}, {((-I)*EL*IndexDelta[j1, j2]*Mass[F[2, {j1}]])/(Sqrt[2]*MW*SW)}}, 
 C[-F[2, {j2}], F[1, {j1}], V[3]] == 
  {{((-I)*EL*$HKSign*IndexDelta[j1, j2])/(Sqrt[2]*SW)}, {0}}, 
 C[-F[1, {j2}], F[2, {j1}], -V[3]] == 
  {{((-I)*EL*$HKSign*IndexDelta[j1, j2])/(Sqrt[2]*SW)}, {0}}, 
 C[F[1, {j1}], -F[2, {j2}], S[5]] == 
  {{(I*EL*TB*IndexDelta[j1, j2]*Mass[F[2, {j2}]])/(Sqrt[2]*MW*SW)}, {0}}, 
 C[F[2, {j1}], -F[1, {j2}], -S[5]] == 
  {{0}, {(I*EL*TB*IndexDelta[j1, j2]*Mass[F[2, {j1}]])/(Sqrt[2]*MW*SW)}}, 
 C[F[3, {j1, o1}], -F[4, {j2, o2}], S[6]] == 
  {{((-I)*EL*Conjugate[CKM[j1, j2]]*IndexDelta[o1, o2]*Mass[F[4, {j2, o1}]])/
     (Sqrt[2]*MW*SW)}, {(I*EL*Conjugate[CKM[j1, j2]]*IndexDelta[o1, o2]*
      Mass[F[3, {j1, o1}]])/(Sqrt[2]*MW*SW)}}, 
 C[F[4, {j1, o1}], -F[3, {j2, o2}], -S[6]] == 
  {{(I*EL*CKM[j2, j1]*IndexDelta[o1, o2]*Mass[F[3, {j2, o1}]])/
     (Sqrt[2]*MW*SW)}, {((-I)*EL*CKM[j2, j1]*IndexDelta[o1, o2]*
      Mass[F[4, {j1, o1}]])/(Sqrt[2]*MW*SW)}}, 
 C[-F[4, {j2, o1}], F[3, {j1, o2}], V[3]] == 
  {{((-I)*EL*$HKSign*Conjugate[CKM[j1, j2]]*IndexDelta[o1, o2])/
     (Sqrt[2]*SW)}, {0}}, C[-F[3, {j2, o1}], F[4, {j1, o2}], -V[3]] == 
  {{((-I)*EL*$HKSign*CKM[j2, j1]*IndexDelta[o1, o2])/(Sqrt[2]*SW)}, {0}}, 
 C[F[3, {j1, o1}], -F[4, {j2, o2}], S[5]] == 
  {{(I*EL*TB*Conjugate[CKM[j1, j2]]*IndexDelta[o1, o2]*Mass[F[4, {j2, o1}]])/
     (Sqrt[2]*MW*SW)}, {(I*EL*Conjugate[CKM[j1, j2]]*IndexDelta[o1, o2]*
      Mass[F[3, {j1, o1}]])/(Sqrt[2]*MW*SW*TB)}}, 
 C[F[4, {j1, o1}], -F[3, {j2, o2}], -S[5]] == 
  {{(I*EL*CKM[j2, j1]*IndexDelta[o1, o2]*Mass[F[3, {j2, o1}]])/
     (Sqrt[2]*MW*SW*TB)}, {(I*EL*TB*CKM[j2, j1]*IndexDelta[o1, o2]*
      Mass[F[4, {j1, o1}]])/(Sqrt[2]*MW*SW)}}, 
 C[S[3], S[12, {s1, j1}], -S[12, {s2, j2}]] == 
  {{-(EL*IndexDelta[j1, j2]*Mass[F[2, {j1}]]*
       ((MUE + TB*Conjugate[Af[2, j1, j1]])*Conjugate[USf[2, j1][s1, 2]]*
         USf[2, j1][s2, 1] - (TB*Af[2, j1, j1] + Conjugate[MUE])*
         Conjugate[USf[2, j1][s1, 1]]*USf[2, j1][s2, 2]))/(2*MW*SW)}}, 
 C[S[4], S[12, {s1, j1}], -S[12, {s2, j2}]] == 
  {{-(EL*IndexDelta[j1, j2]*Mass[F[2, {j1}]]*
       ((MUE*TB - Conjugate[Af[2, j1, j1]])*Conjugate[USf[2, j1][s1, 2]]*
         USf[2, j1][s2, 1] + (Af[2, j1, j1] - TB*Conjugate[MUE])*
         Conjugate[USf[2, j1][s1, 1]]*USf[2, j1][s2, 2]))/(2*MW*SW)}}, 
 C[S[3], S[13, {s1, j1, o1}], -S[13, {s2, j2, o2}]] == 
  {{-(EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[3, {j1}]]*
       ((MUE*TB + Conjugate[Af[3, j1, j1]])*Conjugate[USf[3, j1][s1, 2]]*
         USf[3, j1][s2, 1] - (Af[3, j1, j1] + TB*Conjugate[MUE])*
         Conjugate[USf[3, j1][s1, 1]]*USf[3, j1][s2, 2]))/(2*MW*SW*TB)}}, 
 C[S[4], S[13, {s1, j1, o1}], -S[13, {s2, j2, o2}]] == 
  {{(EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[3, {j1}]]*
      ((MUE - TB*Conjugate[Af[3, j1, j1]])*Conjugate[USf[3, j1][s1, 2]]*
        USf[3, j1][s2, 1] + (TB*Af[3, j1, j1] - Conjugate[MUE])*
        Conjugate[USf[3, j1][s1, 1]]*USf[3, j1][s2, 2]))/(2*MW*SW*TB)}}, 
 C[S[3], S[14, {s1, j1, o1}], -S[14, {s2, j2, o2}]] == 
  {{-(EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[4, {j1}]]*
       ((MUE + TB*Conjugate[Af[4, j1, j1]])*Conjugate[USf[4, j1][s1, 2]]*
         USf[4, j1][s2, 1] - (TB*Af[4, j1, j1] + Conjugate[MUE])*
         Conjugate[USf[4, j1][s1, 1]]*USf[4, j1][s2, 2]))/(2*MW*SW)}}, 
 C[S[4], S[14, {s1, j1, o1}], -S[14, {s2, j2, o2}]] == 
  {{-(EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[4, {j1}]]*
       ((MUE*TB - Conjugate[Af[4, j1, j1]])*Conjugate[USf[4, j1][s1, 2]]*
         USf[4, j1][s2, 1] + (Af[4, j1, j1] - TB*Conjugate[MUE])*
         Conjugate[USf[4, j1][s1, 1]]*USf[4, j1][s2, 2]))/(2*MW*SW)}}, 
 C[S[1], S[11, {j1}], -S[11, {j2}]] == 
  {{((I/2)*EL*MZ*SAB*IndexDelta[j1, j2])/(CW*SW)}}, 
 C[S[2], S[11, {j1}], -S[11, {j2}]] == 
  {{((-I/2)*CAB*EL*MZ*IndexDelta[j1, j2])/(CW*SW)}}, 
 C[S[1], S[12, {s1, j1}], -S[12, {s2, j2}]] == 
  {{((I/2)*EL*IndexDelta[j1, j2]*(Conjugate[USf[2, j1][s1, 1]]*
        ((CB*MW*MZ*SAB*(-1 + 2*SW^2) + 2*CW*SA*Mass[F[2, {j1}]]^2)*
          USf[2, j1][s2, 1] + CW*(SA*Af[2, j1, j1] + CA*Conjugate[MUE])*
          Mass[F[2, {j1}]]*USf[2, j1][s2, 2]) + Conjugate[USf[2, j1][s1, 2]]*
        (CW*(CA*MUE + SA*Conjugate[Af[2, j1, j1]])*Mass[F[2, {j1}]]*
          USf[2, j1][s2, 1] - 2*CB*MW*MZ*SAB*SW^2*USf[2, j1][s2, 2] + 
         2*CW*SA*Mass[F[2, {j1}]]^2*USf[2, j1][s2, 2])))/(CB*CW*MW*SW)}}, 
 C[S[2], S[12, {s1, j1}], -S[12, {s2, j2}]] == 
  {{((I/2)*EL*IndexDelta[j1, j2]*(Conjugate[USf[2, j1][s1, 1]]*
        ((CAB*CB*MW*MZ*(1 - 2*SW^2) - 2*CA*CW*Mass[F[2, {j1}]]^2)*
          USf[2, j1][s2, 1] + CW*(-(CA*Af[2, j1, j1]) + SA*Conjugate[MUE])*
          Mass[F[2, {j1}]]*USf[2, j1][s2, 2]) + Conjugate[USf[2, j1][s1, 2]]*
        (CW*(MUE*SA - CA*Conjugate[Af[2, j1, j1]])*Mass[F[2, {j1}]]*
          USf[2, j1][s2, 1] + 2*CAB*CB*MW*MZ*SW^2*USf[2, j1][s2, 2] - 
         2*CA*CW*Mass[F[2, {j1}]]^2*USf[2, j1][s2, 2])))/(CB*CW*MW*SW)}}, 
 C[S[1], S[13, {s1, j1, o1}], -S[13, {s2, j2, o2}]] == 
  {{((-I/6)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      (Conjugate[USf[3, j1][s1, 1]]*((MW*MZ*SAB*SB*(-3 + 4*SW^2) + 
           6*CA*CW*Mass[F[3, {j1}]]^2)*USf[3, j1][s2, 1] + 
         3*CW*(CA*Af[3, j1, j1] + SA*Conjugate[MUE])*Mass[F[3, {j1}]]*
          USf[3, j1][s2, 2]) + Conjugate[USf[3, j1][s1, 2]]*
        (3*CW*(MUE*SA + CA*Conjugate[Af[3, j1, j1]])*Mass[F[3, {j1}]]*
          USf[3, j1][s2, 1] - 4*MW*MZ*SAB*SB*SW^2*USf[3, j1][s2, 2] + 
         6*CA*CW*Mass[F[3, {j1}]]^2*USf[3, j1][s2, 2])))/(CW*MW*SB*SW)}}, 
 C[S[2], S[13, {s1, j1, o1}], -S[13, {s2, j2, o2}]] == 
  {{((-I/6)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      (Conjugate[USf[3, j1][s1, 1]]*((CAB*MW*MZ*SB*(3 - 4*SW^2) + 
           6*CW*SA*Mass[F[3, {j1}]]^2)*USf[3, j1][s2, 1] + 
         3*CW*(SA*Af[3, j1, j1] - CA*Conjugate[MUE])*Mass[F[3, {j1}]]*
          USf[3, j1][s2, 2]) + Conjugate[USf[3, j1][s1, 2]]*
        (3*CW*(-(CA*MUE) + SA*Conjugate[Af[3, j1, j1]])*Mass[F[3, {j1}]]*
          USf[3, j1][s2, 1] + 4*CAB*MW*MZ*SB*SW^2*USf[3, j1][s2, 2] + 
         6*CW*SA*Mass[F[3, {j1}]]^2*USf[3, j1][s2, 2])))/(CW*MW*SB*SW)}}, 
 C[S[1], S[14, {s1, j1, o1}], -S[14, {s2, j2, o2}]] == 
  {{((I/6)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      (Conjugate[USf[4, j1][s1, 1]]*((CB*MW*MZ*SAB*(-3 + 2*SW^2) + 
           6*CW*SA*Mass[F[4, {j1}]]^2)*USf[4, j1][s2, 1] + 
         3*CW*(SA*Af[4, j1, j1] + CA*Conjugate[MUE])*Mass[F[4, {j1}]]*
          USf[4, j1][s2, 2]) + Conjugate[USf[4, j1][s1, 2]]*
        (3*CW*(CA*MUE + SA*Conjugate[Af[4, j1, j1]])*Mass[F[4, {j1}]]*
          USf[4, j1][s2, 1] - 2*CB*MW*MZ*SAB*SW^2*USf[4, j1][s2, 2] + 
         6*CW*SA*Mass[F[4, {j1}]]^2*USf[4, j1][s2, 2])))/(CB*CW*MW*SW)}}, 
 C[S[2], S[14, {s1, j1, o1}], -S[14, {s2, j2, o2}]] == 
  {{((-I/6)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      (Conjugate[USf[4, j1][s1, 1]]*((CAB*CB*MW*MZ*(-3 + 2*SW^2) + 
           6*CA*CW*Mass[F[4, {j1}]]^2)*USf[4, j1][s2, 1] + 
         3*CW*(CA*Af[4, j1, j1] - SA*Conjugate[MUE])*Mass[F[4, {j1}]]*
          USf[4, j1][s2, 2]) + Conjugate[USf[4, j1][s1, 2]]*
        (3*CW*(-(MUE*SA) + CA*Conjugate[Af[4, j1, j1]])*Mass[F[4, {j1}]]*
          USf[4, j1][s2, 1] - 2*CAB*CB*MW*MZ*SW^2*USf[4, j1][s2, 2] + 
         6*CA*CW*Mass[F[4, {j1}]]^2*USf[4, j1][s2, 2])))/(CB*CW*MW*SW)}}, 
 C[-S[5], S[14, {s2, j2, o1}], -S[13, {s1, j1, o2}]] == 
  {{((-I)*EL*CKM[j1, j2]*IndexDelta[o1, o2]*
      (-(Conjugate[USf[4, j2][s2, 2]]*Mass[F[4, {j2}]]*
         (TB*(MUE + TB*Conjugate[Af[4, j2, j2]])*USf[3, j1][s1, 1] + 
          (1 + TB^2)*Mass[F[3, {j1}]]*USf[3, j1][s1, 2])) + 
       Conjugate[USf[4, j2][s2, 1]]*
        (-((Mass[F[3, {j1}]]^2 + TB*(-(MW^2*S2B) + TB*Mass[F[4, {j2}]]^2))*
           USf[3, j1][s1, 1]) - (Af[3, j1, j1] + TB*Conjugate[MUE])*
          Mass[F[3, {j1}]]*USf[3, j1][s1, 2])))/(Sqrt[2]*MW*SW*TB)}}, 
 C[S[5], S[13, {s1, j1, o1}], -S[14, {s2, j2, o2}]] == 
  {{((-I)*EL*Conjugate[CKM[j1, j2]]*IndexDelta[o1, o2]*
      (-(Conjugate[USf[3, j1][s1, 2]]*Mass[F[3, {j1}]]*
         ((MUE*TB + Conjugate[Af[3, j1, j1]])*USf[4, j2][s2, 1] + 
          (1 + TB^2)*Mass[F[4, {j2}]]*USf[4, j2][s2, 2])) + 
       Conjugate[USf[3, j1][s1, 1]]*
        (-((Mass[F[3, {j1}]]^2 + TB*(-(MW^2*S2B) + TB*Mass[F[4, {j2}]]^2))*
           USf[4, j2][s2, 1]) - TB*(TB*Af[4, j2, j2] + Conjugate[MUE])*
          Mass[F[4, {j2}]]*USf[4, j2][s2, 2])))/(Sqrt[2]*MW*SW*TB)}}, 
 C[-S[5], S[12, {s2, j2}], -S[11, {j1}]] == 
  {{((-I)*EL*IndexDelta[j1, j2]*(-((MUE + TB*Conjugate[Af[2, j1, j1]])*
         Conjugate[USf[2, j1][s2, 2]]*Mass[F[2, {j1}]]) + 
       Conjugate[USf[2, j1][s2, 1]]*(MW^2*S2B - TB*Mass[F[2, {j1}]]^2)))/
     (Sqrt[2]*MW*SW)}}, C[S[5], S[11, {j1}], -S[12, {s2, j2}]] == 
  {{((-I)*EL*IndexDelta[j1, j2]*((MW^2*S2B - TB*Mass[F[2, {j1}]]^2)*
        USf[2, j1][s2, 1] - (TB*Af[2, j1, j1] + Conjugate[MUE])*
        Mass[F[2, {j1}]]*USf[2, j1][s2, 2]))/(Sqrt[2]*MW*SW)}}, 
 C[-S[6], S[14, {s2, j2, o1}], -S[13, {s1, j1, o2}]] == 
  {{(I*EL*CKM[j1, j2]*IndexDelta[o1, o2]*
      (TB*(MUE*TB - Conjugate[Af[4, j2, j2]])*Conjugate[USf[4, j2][s2, 2]]*
        Mass[F[4, {j2}]]*USf[3, j1][s1, 1] + Conjugate[USf[4, j2][s2, 1]]*
        (TB*(C2B*MW^2 + Mass[F[3, {j1}]]^2 - Mass[F[4, {j2}]]^2)*
          USf[3, j1][s1, 1] + (TB*Af[3, j1, j1] - Conjugate[MUE])*
          Mass[F[3, {j1}]]*USf[3, j1][s1, 2])))/(Sqrt[2]*MW*SW*TB)}}, 
 C[S[6], S[13, {s1, j1, o1}], -S[14, {s2, j2, o2}]] == 
  {{(I*EL*Conjugate[CKM[j1, j2]]*IndexDelta[o1, o2]*
      ((-MUE + TB*Conjugate[Af[3, j1, j1]])*Conjugate[USf[3, j1][s1, 2]]*
        Mass[F[3, {j1}]]*USf[4, j2][s2, 1] + TB*Conjugate[USf[3, j1][s1, 1]]*
        ((C2B*MW^2 + Mass[F[3, {j1}]]^2 - Mass[F[4, {j2}]]^2)*
          USf[4, j2][s2, 1] + (-Af[4, j2, j2] + TB*Conjugate[MUE])*
          Mass[F[4, {j2}]]*USf[4, j2][s2, 2])))/(Sqrt[2]*MW*SW*TB)}}, 
 C[-S[6], S[12, {s2, j2}], -S[11, {j1}]] == 
  {{(I*EL*IndexDelta[j1, j2]*((MUE*TB - Conjugate[Af[2, j1, j1]])*
        Conjugate[USf[2, j1][s2, 2]]*Mass[F[2, {j1}]] + 
       Conjugate[USf[2, j1][s2, 1]]*(C2B*MW^2 - Mass[F[2, {j1}]]^2)))/
     (Sqrt[2]*MW*SW)}}, C[S[6], S[11, {j1}], -S[12, {s2, j2}]] == 
  {{(I*EL*IndexDelta[j1, j2]*((C2B*MW^2 - Mass[F[2, {j1}]]^2)*
        USf[2, j1][s2, 1] + (-Af[2, j1, j1] + TB*Conjugate[MUE])*
        Mass[F[2, {j1}]]*USf[2, j1][s2, 2]))/(Sqrt[2]*MW*SW)}}, 
 C[S[11, {j1}], -S[11, {j2}], V[2]] == 
  {{((-I/2)*EL*$HKSign*IndexDelta[j1, j2])/(CW*SW)}}, 
 C[S[12, {s1, j1}], -S[12, {s2, j2}], V[1]] == 
  {{I*EL*IndexDelta[j1, j2]*IndexDelta[s1, s2]}}, 
 C[S[12, {s1, j1}], -S[12, {s2, j2}], V[2]] == 
  {{((-I/2)*EL*$HKSign*IndexDelta[j1, j2]*
      ((-1 + 2*SW^2)*Conjugate[USf[2, j1][s1, 1]]*USf[2, j1][s2, 1] + 
       2*SW^2*Conjugate[USf[2, j1][s1, 2]]*USf[2, j1][s2, 2]))/(CW*SW)}}, 
 C[S[13, {s1, j1, o1}], -S[13, {s2, j2, o2}], V[1]] == 
  {{((-2*I)/3)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*IndexDelta[s1, s2]}}, 
 C[S[13, {s1, j1, o1}], -S[13, {s2, j2, o2}], V[2]] == 
  {{((I/6)*EL*$HKSign*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      ((-3 + 4*SW^2)*Conjugate[USf[3, j1][s1, 1]]*USf[3, j1][s2, 1] + 
       4*SW^2*Conjugate[USf[3, j1][s1, 2]]*USf[3, j1][s2, 2]))/(CW*SW)}}, 
 C[S[14, {s1, j1, o1}], -S[14, {s2, j2, o2}], V[1]] == 
  {{(I/3)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*IndexDelta[s1, s2]}}, 
 C[S[14, {s1, j1, o1}], -S[14, {s2, j2, o2}], V[2]] == 
  {{((-I/6)*EL*$HKSign*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      ((-3 + 2*SW^2)*Conjugate[USf[4, j1][s1, 1]]*USf[4, j1][s2, 1] + 
       2*SW^2*Conjugate[USf[4, j1][s1, 2]]*USf[4, j1][s2, 2]))/(CW*SW)}}, 
 C[S[13, {s1, j1, o1}], -S[14, {s2, j2, o2}], V[3]] == 
  {{((-I)*EL*$HKSign*Conjugate[CKM[j1, j2]]*Conjugate[USf[3, j1][s1, 1]]*
      IndexDelta[o1, o2]*USf[4, j2][s2, 1])/(Sqrt[2]*SW)}}, 
 C[S[14, {s2, j2, o1}], -S[13, {s1, j1, o2}], -V[3]] == 
  {{((-I)*EL*$HKSign*CKM[j1, j2]*Conjugate[USf[4, j2][s2, 1]]*
      IndexDelta[o1, o2]*USf[3, j1][s1, 1])/(Sqrt[2]*SW)}}, 
 C[S[11, {j1}], -S[12, {s2, j2}], V[3]] == 
  {{((-I)*EL*$HKSign*IndexDelta[j1, j2]*USf[2, j1][s2, 1])/(Sqrt[2]*SW)}}, 
 C[S[12, {s2, j2}], -S[11, {j1}], -V[3]] == 
  {{((-I)*EL*$HKSign*Conjugate[USf[2, j1][s2, 1]]*IndexDelta[j1, j2])/
     (Sqrt[2]*SW)}}, C[F[11, {n2}], F[11, {n1}], S[1]] == 
  {{((-I/2)*EL*$HKSign*(SA*Conjugate[ZNeu[n1, 3]]*
        (SW*Conjugate[ZNeu[n2, 1]] - CW*Conjugate[ZNeu[n2, 2]]) + 
       CA*Conjugate[ZNeu[n1, 4]]*(SW*Conjugate[ZNeu[n2, 1]] - 
         CW*Conjugate[ZNeu[n2, 2]]) + (SW*Conjugate[ZNeu[n1, 1]] - 
         CW*Conjugate[ZNeu[n1, 2]])*(SA*Conjugate[ZNeu[n2, 3]] + 
         CA*Conjugate[ZNeu[n2, 4]])))/(CW*SW)}, 
   {((I/2)*EL*$HKSign*(ZNeu[n1, 4]*(-(CA*SW*ZNeu[n2, 1]) + 
         CA*CW*ZNeu[n2, 2]) + ZNeu[n1, 3]*(-(SA*SW*ZNeu[n2, 1]) + 
         CW*SA*ZNeu[n2, 2]) - (SW*ZNeu[n1, 1] - CW*ZNeu[n1, 2])*
        (SA*ZNeu[n2, 3] + CA*ZNeu[n2, 4])))/(CW*SW)}}, 
 C[F[11, {n2}], F[11, {n1}], S[2]] == 
  {{((I/2)*EL*$HKSign*(CA*Conjugate[ZNeu[n1, 3]]*(SW*Conjugate[ZNeu[n2, 1]] - 
         CW*Conjugate[ZNeu[n2, 2]]) + Conjugate[ZNeu[n1, 4]]*
        (-(SA*SW*Conjugate[ZNeu[n2, 1]]) + CW*SA*Conjugate[ZNeu[n2, 2]]) + 
       (SW*Conjugate[ZNeu[n1, 1]] - CW*Conjugate[ZNeu[n1, 2]])*
        (CA*Conjugate[ZNeu[n2, 3]] - SA*Conjugate[ZNeu[n2, 4]])))/(CW*SW)}, 
   {((I/2)*EL*$HKSign*(CA*ZNeu[n1, 3]*(SW*ZNeu[n2, 1] - CW*ZNeu[n2, 2]) + 
       ZNeu[n1, 4]*(-(SA*SW*ZNeu[n2, 1]) + CW*SA*ZNeu[n2, 2]) + 
       (SW*ZNeu[n1, 1] - CW*ZNeu[n1, 2])*(CA*ZNeu[n2, 3] - SA*ZNeu[n2, 4])))/
     (CW*SW)}}, C[F[11, {n2}], F[11, {n1}], S[3]] == 
  {{(EL*$HKSign*(SB*Conjugate[ZNeu[n1, 3]]*(SW*Conjugate[ZNeu[n2, 1]] - 
         CW*Conjugate[ZNeu[n2, 2]]) + Conjugate[ZNeu[n1, 4]]*
        (-(CB*SW*Conjugate[ZNeu[n2, 1]]) + CB*CW*Conjugate[ZNeu[n2, 2]]) + 
       (SW*Conjugate[ZNeu[n1, 1]] - CW*Conjugate[ZNeu[n1, 2]])*
        (SB*Conjugate[ZNeu[n2, 3]] - CB*Conjugate[ZNeu[n2, 4]])))/(2*CW*SW)}, 
   {-(EL*$HKSign*(SB*ZNeu[n1, 3]*(SW*ZNeu[n2, 1] - CW*ZNeu[n2, 2]) + 
        ZNeu[n1, 4]*(-(CB*SW*ZNeu[n2, 1]) + CB*CW*ZNeu[n2, 2]) + 
        (SW*ZNeu[n1, 1] - CW*ZNeu[n1, 2])*(SB*ZNeu[n2, 3] - CB*ZNeu[n2, 4])))/
     (2*CW*SW)}}, C[F[11, {n2}], F[11, {n1}], S[4]] == 
  {{-(EL*$HKSign*(CB*Conjugate[ZNeu[n1, 3]]*(SW*Conjugate[ZNeu[n2, 1]] - 
          CW*Conjugate[ZNeu[n2, 2]]) + SB*Conjugate[ZNeu[n1, 4]]*
         (SW*Conjugate[ZNeu[n2, 1]] - CW*Conjugate[ZNeu[n2, 2]]) + 
        (SW*Conjugate[ZNeu[n1, 1]] - CW*Conjugate[ZNeu[n1, 2]])*
         (CB*Conjugate[ZNeu[n2, 3]] + SB*Conjugate[ZNeu[n2, 4]])))/
     (2*CW*SW)}, 
   {(EL*$HKSign*(CB*ZNeu[n1, 3]*(SW*ZNeu[n2, 1] - CW*ZNeu[n2, 2]) + 
       SB*ZNeu[n1, 4]*(SW*ZNeu[n2, 1] - CW*ZNeu[n2, 2]) + 
       (SW*ZNeu[n1, 1] - CW*ZNeu[n1, 2])*(CB*ZNeu[n2, 3] + SB*ZNeu[n2, 4])))/
     (2*CW*SW)}}, C[F[12, {c1}], -F[12, {c2}], S[1]] == 
  {{(I*EL*(SA*Conjugate[UCha[c1, 2]]*Conjugate[VCha[c2, 1]] - 
       CA*Conjugate[UCha[c1, 1]]*Conjugate[VCha[c2, 2]]))/(Sqrt[2]*SW)}, 
   {(I*EL*(SA*UCha[c2, 2]*VCha[c1, 1] - CA*UCha[c2, 1]*VCha[c1, 2]))/
     (Sqrt[2]*SW)}}, C[F[12, {c1}], -F[12, {c2}], S[2]] == 
  {{((-I)*EL*(CA*Conjugate[UCha[c1, 2]]*Conjugate[VCha[c2, 1]] + 
       SA*Conjugate[UCha[c1, 1]]*Conjugate[VCha[c2, 2]]))/(Sqrt[2]*SW)}, 
   {((-I)*EL*(CA*UCha[c2, 2]*VCha[c1, 1] + SA*UCha[c2, 1]*VCha[c1, 2]))/
     (Sqrt[2]*SW)}}, C[F[12, {c1}], -F[12, {c2}], S[3]] == 
  {{-((EL*(SB*Conjugate[UCha[c1, 2]]*Conjugate[VCha[c2, 1]] + 
        CB*Conjugate[UCha[c1, 1]]*Conjugate[VCha[c2, 2]]))/(Sqrt[2]*SW))}, 
   {(EL*(SB*UCha[c2, 2]*VCha[c1, 1] + CB*UCha[c2, 1]*VCha[c1, 2]))/
     (Sqrt[2]*SW)}}, C[F[12, {c1}], -F[12, {c2}], S[4]] == 
  {{(EL*(CB*Conjugate[UCha[c1, 2]]*Conjugate[VCha[c2, 1]] - 
       SB*Conjugate[UCha[c1, 1]]*Conjugate[VCha[c2, 2]]))/(Sqrt[2]*SW)}, 
   {(EL*(-(CB*UCha[c2, 2]*VCha[c1, 1]) + SB*UCha[c2, 1]*VCha[c1, 2]))/
     (Sqrt[2]*SW)}}, C[F[11, {n1}], -F[12, {c2}], S[5]] == 
  {{((-I)*CB*EL*$HKSign*((Conjugate[VCha[c2, 2]]*
         ((SW*Conjugate[ZNeu[n1, 1]])/CW + Conjugate[ZNeu[n1, 2]]))/Sqrt[2] + 
       Conjugate[VCha[c2, 1]]*Conjugate[ZNeu[n1, 4]]))/SW}, 
   {((-I)*EL*SB*$HKSign*(-((UCha[c2, 2]*((SW*ZNeu[n1, 1])/CW + ZNeu[n1, 2]))/
         Sqrt[2]) + UCha[c2, 1]*ZNeu[n1, 3]))/SW}}, 
 C[F[11, {n1}], -F[12, {c2}], S[6]] == 
  {{((-I)*EL*SB*$HKSign*((Conjugate[VCha[c2, 2]]*
         ((SW*Conjugate[ZNeu[n1, 1]])/CW + Conjugate[ZNeu[n1, 2]]))/Sqrt[2] + 
       Conjugate[VCha[c2, 1]]*Conjugate[ZNeu[n1, 4]]))/SW}, 
   {(I*CB*EL*$HKSign*(-((UCha[c2, 2]*((SW*ZNeu[n1, 1])/CW + ZNeu[n1, 2]))/
         Sqrt[2]) + UCha[c2, 1]*ZNeu[n1, 3]))/SW}}, 
 C[F[12, {c2}], F[11, {n1}], -S[5]] == 
  {{((-I)*EL*SB*$HKSign*(-((Conjugate[UCha[c2, 2]]*
          ((SW*Conjugate[ZNeu[n1, 1]])/CW + Conjugate[ZNeu[n1, 2]]))/
         Sqrt[2]) + Conjugate[UCha[c2, 1]]*Conjugate[ZNeu[n1, 3]]))/SW}, 
   {((-I)*CB*EL*$HKSign*((VCha[c2, 2]*((SW*ZNeu[n1, 1])/CW + ZNeu[n1, 2]))/
        Sqrt[2] + VCha[c2, 1]*ZNeu[n1, 4]))/SW}}, 
 C[F[12, {c2}], F[11, {n1}], -S[6]] == 
  {{(I*CB*EL*$HKSign*(-((Conjugate[UCha[c2, 2]]*((SW*Conjugate[ZNeu[n1, 1]])/
            CW + Conjugate[ZNeu[n1, 2]]))/Sqrt[2]) + Conjugate[UCha[c2, 1]]*
        Conjugate[ZNeu[n1, 3]]))/SW}, 
   {((-I)*EL*SB*$HKSign*((VCha[c2, 2]*((SW*ZNeu[n1, 1])/CW + ZNeu[n1, 2]))/
        Sqrt[2] + VCha[c2, 1]*ZNeu[n1, 4]))/SW}}, 
 C[F[11, {n1}], -F[1, {j1}], S[11, {j2}]] == 
  {{0}, {(I*EL*$HKSign*IndexDelta[j1, j2]*(SW*ZNeu[n1, 1] - CW*ZNeu[n1, 2]))/
     (Sqrt[2]*CW*SW)}}, C[F[11, {n1}], -F[2, {j1}], S[12, {s2, j2}]] == 
  {{((-I)*EL*$HKSign*IndexDelta[j1, j2]*(2*CB*MW*SW*Conjugate[ZNeu[n1, 1]]*
        Conjugate[USf[2, j1][s2, 2]] + CW*Conjugate[ZNeu[n1, 3]]*
        Conjugate[USf[2, j1][s2, 1]]*Mass[F[2, {j1}]]))/
     (Sqrt[2]*CB*CW*MW*SW)}, 
   {(I*EL*$HKSign*IndexDelta[j1, j2]*(CB*MW*Conjugate[USf[2, j1][s2, 1]]*
        (SW*ZNeu[n1, 1] + CW*ZNeu[n1, 2]) - CW*Conjugate[USf[2, j1][s2, 2]]*
        Mass[F[2, {j1}]]*ZNeu[n1, 3]))/(Sqrt[2]*CB*CW*MW*SW)}}, 
 C[F[11, {n1}], -F[3, {j1, o1}], S[13, {s2, j2, o2}]] == 
  {{((I/3)*EL*$HKSign*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      (4*MW*SB*SW*Conjugate[ZNeu[n1, 1]]*Conjugate[USf[3, j1][s2, 2]] - 
       3*CW*Conjugate[ZNeu[n1, 4]]*Conjugate[USf[3, j1][s2, 1]]*
        Mass[F[3, {j1, o1}]]))/(Sqrt[2]*CW*MW*SB*SW)}, 
   {((-I/3)*EL*$HKSign*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      (MW*SB*Conjugate[USf[3, j1][s2, 1]]*(SW*ZNeu[n1, 1] + 
         3*CW*ZNeu[n1, 2]) + 3*CW*Conjugate[USf[3, j1][s2, 2]]*
        Mass[F[3, {j1, o1}]]*ZNeu[n1, 4]))/(Sqrt[2]*CW*MW*SB*SW)}}, 
 C[F[11, {n1}], -F[4, {j1, o1}], S[14, {s2, j2, o2}]] == 
  {{((-I/3)*EL*$HKSign*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      (2*CB*MW*SW*Conjugate[ZNeu[n1, 1]]*Conjugate[USf[4, j1][s2, 2]] + 
       3*CW*Conjugate[ZNeu[n1, 3]]*Conjugate[USf[4, j1][s2, 1]]*
        Mass[F[4, {j1, o1}]]))/(Sqrt[2]*CB*CW*MW*SW)}, 
   {((I/3)*EL*$HKSign*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      (CB*MW*Conjugate[USf[4, j1][s2, 1]]*(-(SW*ZNeu[n1, 1]) + 
         3*CW*ZNeu[n1, 2]) - 3*CW*Conjugate[USf[4, j1][s2, 2]]*
        Mass[F[4, {j1, o1}]]*ZNeu[n1, 3]))/(Sqrt[2]*CB*CW*MW*SW)}}, 
 C[F[1, {j1}], F[11, {n1}], -S[11, {j2}]] == 
  {{(I*EL*$HKSign*(SW*Conjugate[ZNeu[n1, 1]] - CW*Conjugate[ZNeu[n1, 2]])*
      IndexDelta[j1, j2])/(Sqrt[2]*CW*SW)}, {0}}, 
 C[F[2, {j1}], F[11, {n1}], -S[12, {s2, j2}]] == 
  {{(I*EL*$HKSign*IndexDelta[j1, j2]*(CB*MW*SW*Conjugate[ZNeu[n1, 1]]*
        USf[2, j1][s2, 1] + CW*(CB*MW*Conjugate[ZNeu[n1, 2]]*
          USf[2, j1][s2, 1] - Conjugate[ZNeu[n1, 3]]*Mass[F[2, {j1}]]*
          USf[2, j1][s2, 2])))/(Sqrt[2]*CB*CW*MW*SW)}, 
   {((-I)*EL*$HKSign*IndexDelta[j1, j2]*(CW*Mass[F[2, {j1}]]*ZNeu[n1, 3]*
        USf[2, j1][s2, 1] + 2*CB*MW*SW*ZNeu[n1, 1]*USf[2, j1][s2, 2]))/
     (Sqrt[2]*CB*CW*MW*SW)}}, C[F[3, {j1, o1}], F[11, {n1}], 
   -S[13, {s2, j2, o2}]] == 
  {{((-I/3)*EL*$HKSign*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      (MW*SB*SW*Conjugate[ZNeu[n1, 1]]*USf[3, j1][s2, 1] + 
       3*CW*(MW*SB*Conjugate[ZNeu[n1, 2]]*USf[3, j1][s2, 1] + 
         Conjugate[ZNeu[n1, 4]]*Mass[F[3, {j1, o1}]]*USf[3, j1][s2, 2])))/
     (Sqrt[2]*CW*MW*SB*SW)}, {((I/3)*EL*$HKSign*IndexDelta[j1, j2]*
      IndexDelta[o1, o2]*(-3*CW*Mass[F[3, {j1, o1}]]*ZNeu[n1, 4]*
        USf[3, j1][s2, 1] + 4*MW*SB*SW*ZNeu[n1, 1]*USf[3, j1][s2, 2]))/
     (Sqrt[2]*CW*MW*SB*SW)}}, C[F[4, {j1, o1}], F[11, {n1}], 
   -S[14, {s2, j2, o2}]] == 
  {{((I/3)*EL*$HKSign*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      (-(CB*MW*SW*Conjugate[ZNeu[n1, 1]]*USf[4, j1][s2, 1]) + 
       3*CW*(CB*MW*Conjugate[ZNeu[n1, 2]]*USf[4, j1][s2, 1] - 
         Conjugate[ZNeu[n1, 3]]*Mass[F[4, {j1, o1}]]*USf[4, j1][s2, 2])))/
     (Sqrt[2]*CB*CW*MW*SW)}, {((-I/3)*EL*$HKSign*IndexDelta[j1, j2]*
      IndexDelta[o1, o2]*(3*CW*Mass[F[4, {j1, o1}]]*ZNeu[n1, 3]*
        USf[4, j1][s2, 1] + 2*CB*MW*SW*ZNeu[n1, 1]*USf[4, j1][s2, 2]))/
     (Sqrt[2]*CB*CW*MW*SW)}}, C[F[12, {c1}], -F[4, {j2, o1}], 
   S[13, {s1, j1, o2}]] == 
  {{(I*EL*$HKSign*Conjugate[CKM[j1, j2]]*Conjugate[UCha[c1, 2]]*
      Conjugate[USf[3, j1][s1, 1]]*IndexDelta[o1, o2]*Mass[F[4, {j2, o1}]])/
     (Sqrt[2]*CB*MW*SW)}, {((-I/2)*EL*$HKSign*Conjugate[CKM[j1, j2]]*
      IndexDelta[o1, o2]*(2*MW*SB*Conjugate[USf[3, j1][s1, 1]]*VCha[c1, 1] - 
       Sqrt[2]*Conjugate[USf[3, j1][s1, 2]]*Mass[F[3, {j1}]]*VCha[c1, 2]))/
     (MW*SB*SW)}}, C[-F[12, {c1}], -F[3, {j1, o1}], S[14, {s2, j2, o2}]] == 
  {{(I*EL*$HKSign*CKM[j1, j2]*Conjugate[VCha[c1, 2]]*
      Conjugate[USf[4, j2][s2, 1]]*IndexDelta[o1, o2]*Mass[F[3, {j1, o1}]])/
     (Sqrt[2]*MW*SB*SW)}, {((-I/2)*EL*$HKSign*CKM[j1, j2]*IndexDelta[o1, o2]*
      (2*CB*MW*Conjugate[USf[4, j2][s2, 1]]*UCha[c1, 1] - 
       Sqrt[2]*Conjugate[USf[4, j2][s2, 2]]*Mass[F[4, {j2}]]*UCha[c1, 2]))/
     (CB*MW*SW)}}, C[F[12, {c1}], -F[2, {j2}], S[11, {j1}]] == 
  {{(I*EL*$HKSign*Conjugate[UCha[c1, 2]]*IndexDelta[j1, j2]*Mass[F[2, {j1}]])/
     (Sqrt[2]*CB*MW*SW)}, {((-I)*EL*$HKSign*IndexDelta[j1, j2]*VCha[c1, 1])/
     SW}}, C[-F[12, {c1}], -F[1, {j1}], S[12, {s2, j2}]] == 
  {{0}, {((I/2)*EL*$HKSign*IndexDelta[j1, j2]*
      (-2*Conjugate[USf[2, j1][s2, 1]]*UCha[c1, 1] + 
       (Sqrt[2]*Conjugate[USf[2, j1][s2, 2]]*Mass[F[2, {j1}]]*UCha[c1, 2])/
        (CB*MW)))/SW}}, C[F[4, {j2, o1}], -F[12, {c1}], 
   -S[13, {s1, j1, o2}]] == 
  {{((-I/2)*EL*$HKSign*CKM[j1, j2]*IndexDelta[o1, o2]*
      (2*MW*SB*Conjugate[VCha[c1, 1]]*USf[3, j1][s1, 1] - 
       Sqrt[2]*Conjugate[VCha[c1, 2]]*Mass[F[3, {j1}]]*USf[3, j1][s1, 2]))/
     (MW*SB*SW)}, {(I*EL*$HKSign*CKM[j1, j2]*IndexDelta[o1, o2]*
      Mass[F[4, {j2, o1}]]*UCha[c1, 2]*USf[3, j1][s1, 1])/
     (Sqrt[2]*CB*MW*SW)}}, C[F[3, {j1, o1}], F[12, {c1}], 
   -S[14, {s2, j2, o2}]] == 
  {{((-I/2)*EL*$HKSign*Conjugate[CKM[j1, j2]]*IndexDelta[o1, o2]*
      (2*CB*MW*Conjugate[UCha[c1, 1]]*USf[4, j2][s2, 1] - 
       Sqrt[2]*Conjugate[UCha[c1, 2]]*Mass[F[4, {j2}]]*USf[4, j2][s2, 2]))/
     (CB*MW*SW)}, {(I*EL*$HKSign*Conjugate[CKM[j1, j2]]*IndexDelta[o1, o2]*
      Mass[F[3, {j1, o1}]]*VCha[c1, 2]*USf[4, j2][s2, 1])/
     (Sqrt[2]*MW*SB*SW)}}, C[F[2, {j2}], -F[12, {c1}], -S[11, {j1}]] == 
  {{((-I)*EL*$HKSign*Conjugate[VCha[c1, 1]]*IndexDelta[j1, j2])/SW}, 
   {(I*EL*$HKSign*IndexDelta[j1, j2]*Mass[F[2, {j1}]]*UCha[c1, 2])/
     (Sqrt[2]*CB*MW*SW)}}, C[F[1, {j1}], F[12, {c1}], -S[12, {s2, j2}]] == 
  {{((I/2)*EL*$HKSign*IndexDelta[j1, j2]*(-2*Conjugate[UCha[c1, 1]]*
        USf[2, j1][s2, 1] + (Sqrt[2]*Conjugate[UCha[c1, 2]]*Mass[F[2, {j1}]]*
         USf[2, j1][s2, 2])/(CB*MW)))/SW}, {0}}, 
 C[F[11, {n1}], F[11, {n2}], V[2]] == 
  {{((I/2)*EL*$HKSign*(-(Conjugate[ZNeu[n2, 3]]*ZNeu[n1, 3]) + 
       Conjugate[ZNeu[n2, 4]]*ZNeu[n1, 4]))/(CW*SW)}, 
   {((I/2)*EL*$HKSign*(Conjugate[ZNeu[n1, 3]]*ZNeu[n2, 3] - 
       Conjugate[ZNeu[n1, 4]]*ZNeu[n2, 4]))/(CW*SW)}}, 
 C[F[11, {n2}], -F[12, {c1}], V[3]] == 
  {{(I*EL*$HKSign*(Conjugate[VCha[c1, 1]]*ZNeu[n2, 2] - 
       (Conjugate[VCha[c1, 2]]*ZNeu[n2, 4])/Sqrt[2]))/SW}, 
   {(I*EL*$HKSign*(Conjugate[ZNeu[n2, 2]]*UCha[c1, 1] + 
       (Conjugate[ZNeu[n2, 3]]*UCha[c1, 2])/Sqrt[2]))/SW}}, 
 C[F[12, {c1}], F[11, {n2}], -V[3]] == 
  {{(I*EL*$HKSign*(Conjugate[ZNeu[n2, 2]]*VCha[c1, 1] - 
       (Conjugate[ZNeu[n2, 4]]*VCha[c1, 2])/Sqrt[2]))/SW}, 
   {(I*EL*$HKSign*(Conjugate[UCha[c1, 1]]*ZNeu[n2, 2] + 
       (Conjugate[UCha[c1, 2]]*ZNeu[n2, 3])/Sqrt[2]))/SW}}, 
 C[-F[12, {c2}], F[12, {c1}], V[1]] == {{I*EL*IndexDelta[c1, c2]}, 
   {I*EL*IndexDelta[c1, c2]}}, C[-F[12, {c2}], F[12, {c1}], V[2]] == 
  {{((-I)*EL*$HKSign*(SW^2*IndexDelta[c1, c2] - Conjugate[UCha[c1, 1]]*
        UCha[c2, 1] - (Conjugate[UCha[c1, 2]]*UCha[c2, 2])/2))/(CW*SW)}, 
   {((-I)*EL*$HKSign*(SW^2*IndexDelta[c1, c2] - Conjugate[VCha[c2, 1]]*
        VCha[c1, 1] - (Conjugate[VCha[c2, 2]]*VCha[c1, 2])/2))/(CW*SW)}}, 
 C[S[1], S[1], S[11, {j2}], -S[11, {j1}]] == 
  {{((I/4)*C2A*EL^2*IndexDelta[j1, j2])/(CW^2*SW^2)}}, 
 C[S[1], S[1], S[12, {s2, j2}], -S[12, {s1, j1}]] == 
  {{((I/4)*EL^2*IndexDelta[j1, j2]*(Conjugate[USf[2, j1][s2, 1]]*
        (C2A*CB^2*MW^2*(-1 + 2*SW^2) - 2*CW^2*SA^2*Mass[F[2, {j1}]]^2)*
        USf[2, j1][s1, 1] - 2*Conjugate[USf[2, j1][s2, 2]]*
        (C2A*CB^2*MW^2*SW^2 + CW^2*SA^2*Mass[F[2, {j1}]]^2)*
        USf[2, j1][s1, 2]))/(CB^2*CW^2*MW^2*SW^2)}}, 
 C[S[1], S[1], S[13, {s2, j2, o1}], -S[13, {s1, j1, o2}]] == 
  {{((-I/12)*EL^2*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      (Conjugate[USf[3, j1][s2, 1]]*(C2A*MW^2*SB^2*(-3 + 4*SW^2) + 
         6*CA^2*CW^2*Mass[F[3, {j1}]]^2)*USf[3, j1][s1, 1] + 
       2*Conjugate[USf[3, j1][s2, 2]]*(-2*C2A*MW^2*SB^2*SW^2 + 
         3*CA^2*CW^2*Mass[F[3, {j1}]]^2)*USf[3, j1][s1, 2]))/
     (CW^2*MW^2*SB^2*SW^2)}}, C[S[1], S[1], S[14, {s2, j2, o1}], 
   -S[14, {s1, j1, o2}]] == 
  {{((I/12)*EL^2*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      (Conjugate[USf[4, j1][s2, 1]]*(C2A*CB^2*MW^2*(-3 + 2*SW^2) - 
         6*CW^2*SA^2*Mass[F[4, {j1}]]^2)*USf[4, j1][s1, 1] - 
       2*Conjugate[USf[4, j1][s2, 2]]*(C2A*CB^2*MW^2*SW^2 + 
         3*CW^2*SA^2*Mass[F[4, {j1}]]^2)*USf[4, j1][s1, 2]))/
     (CB^2*CW^2*MW^2*SW^2)}}, C[S[2], S[2], S[11, {j2}], -S[11, {j1}]] == 
  {{((-I/4)*C2A*EL^2*IndexDelta[j1, j2])/(CW^2*SW^2)}}, 
 C[S[2], S[2], S[12, {s2, j2}], -S[12, {s1, j1}]] == 
  {{((-I/4)*EL^2*IndexDelta[j1, j2]*(Conjugate[USf[2, j1][s2, 1]]*
        (C2A*CB^2*MW^2*(-1 + 2*SW^2) + 2*CA^2*CW^2*Mass[F[2, {j1}]]^2)*
        USf[2, j1][s1, 1] + 2*Conjugate[USf[2, j1][s2, 2]]*
        (-(C2A*CB^2*MW^2*SW^2) + CA^2*CW^2*Mass[F[2, {j1}]]^2)*
        USf[2, j1][s1, 2]))/(CB^2*CW^2*MW^2*SW^2)}}, 
 C[S[2], S[2], S[13, {s2, j2, o1}], -S[13, {s1, j1, o2}]] == 
  {{((-I/12)*EL^2*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      (Conjugate[USf[3, j1][s2, 1]]*(C2A*MW^2*SB^2*(3 - 4*SW^2) + 
         6*CW^2*SA^2*Mass[F[3, {j1}]]^2)*USf[3, j1][s1, 1] + 
       2*Conjugate[USf[3, j1][s2, 2]]*(2*C2A*MW^2*SB^2*SW^2 + 
         3*CW^2*SA^2*Mass[F[3, {j1}]]^2)*USf[3, j1][s1, 2]))/
     (CW^2*MW^2*SB^2*SW^2)}}, C[S[2], S[2], S[14, {s2, j2, o1}], 
   -S[14, {s1, j1, o2}]] == 
  {{((-I/12)*EL^2*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      (Conjugate[USf[4, j1][s2, 1]]*(C2A*CB^2*MW^2*(-3 + 2*SW^2) + 
         6*CA^2*CW^2*Mass[F[4, {j1}]]^2)*USf[4, j1][s1, 1] - 
       2*Conjugate[USf[4, j1][s2, 2]]*(C2A*CB^2*MW^2*SW^2 - 
         3*CA^2*CW^2*Mass[F[4, {j1}]]^2)*USf[4, j1][s1, 2]))/
     (CB^2*CW^2*MW^2*SW^2)}}, C[S[3], S[3], S[11, {j2}], -S[11, {j1}]] == 
  {{((I/4)*C2B*EL^2*IndexDelta[j1, j2])/(CW^2*SW^2)}}, 
 C[S[3], S[3], S[12, {s2, j2}], -S[12, {s1, j1}]] == 
  {{((-I/4)*EL^2*IndexDelta[j1, j2]*(Conjugate[USf[2, j1][s2, 1]]*
        (C2B*MW^2*(1 - 2*SW^2) + 2*CW^2*TB^2*Mass[F[2, {j1}]]^2)*
        USf[2, j1][s1, 1] + 2*Conjugate[USf[2, j1][s2, 2]]*
        (C2B*MW^2*SW^2 + CW^2*TB^2*Mass[F[2, {j1}]]^2)*USf[2, j1][s1, 2]))/
     (CW^2*MW^2*SW^2)}}, C[S[3], S[3], S[13, {s2, j2, o1}], 
   -S[13, {s1, j1, o2}]] == 
  {{((-I/12)*EL^2*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      (Conjugate[USf[3, j1][s2, 1]]*(C2B*MW^2*(-3 + 4*SW^2)*TB^2 + 
         6*CW^2*Mass[F[3, {j1}]]^2)*USf[3, j1][s1, 1] + 
       2*Conjugate[USf[3, j1][s2, 2]]*(-2*C2B*MW^2*SW^2*TB^2 + 
         3*CW^2*Mass[F[3, {j1}]]^2)*USf[3, j1][s1, 2]))/
     (CW^2*MW^2*SW^2*TB^2)}}, C[S[3], S[3], S[14, {s2, j2, o1}], 
   -S[14, {s1, j1, o2}]] == 
  {{((-I/12)*EL^2*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      (Conjugate[USf[4, j1][s2, 1]]*(C2B*MW^2*(3 - 2*SW^2) + 
         6*CW^2*TB^2*Mass[F[4, {j1}]]^2)*USf[4, j1][s1, 1] + 
       2*Conjugate[USf[4, j1][s2, 2]]*(C2B*MW^2*SW^2 + 
         3*CW^2*TB^2*Mass[F[4, {j1}]]^2)*USf[4, j1][s1, 2]))/
     (CW^2*MW^2*SW^2)}}, C[S[4], S[4], S[11, {j2}], -S[11, {j1}]] == 
  {{((-I/4)*C2B*EL^2*IndexDelta[j1, j2])/(CW^2*SW^2)}}, 
 C[S[4], S[4], S[12, {s2, j2}], -S[12, {s1, j1}]] == 
  {{((-I/4)*EL^2*IndexDelta[j1, j2]*(Conjugate[USf[2, j1][s2, 1]]*
        (C2B*MW^2*(-1 + 2*SW^2) + 2*CW^2*Mass[F[2, {j1}]]^2)*
        USf[2, j1][s1, 1] + 2*Conjugate[USf[2, j1][s2, 2]]*
        (-(C2B*MW^2*SW^2) + CW^2*Mass[F[2, {j1}]]^2)*USf[2, j1][s1, 2]))/
     (CW^2*MW^2*SW^2)}}, C[S[4], S[4], S[13, {s2, j2, o1}], 
   -S[13, {s1, j1, o2}]] == 
  {{((-I/12)*EL^2*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      (Conjugate[USf[3, j1][s2, 1]]*(C2B*MW^2*(3 - 4*SW^2) + 
         6*CW^2*Mass[F[3, {j1}]]^2)*USf[3, j1][s1, 1] + 
       2*Conjugate[USf[3, j1][s2, 2]]*(2*C2B*MW^2*SW^2 + 
         3*CW^2*Mass[F[3, {j1}]]^2)*USf[3, j1][s1, 2]))/(CW^2*MW^2*SW^2)}}, 
 C[S[4], S[4], S[14, {s2, j2, o1}], -S[14, {s1, j1, o2}]] == 
  {{((-I/12)*EL^2*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      (Conjugate[USf[4, j1][s2, 1]]*(C2B*MW^2*(-3 + 2*SW^2) + 
         6*CW^2*Mass[F[4, {j1}]]^2)*USf[4, j1][s1, 1] - 
       2*Conjugate[USf[4, j1][s2, 2]]*(C2B*MW^2*SW^2 - 
         3*CW^2*Mass[F[4, {j1}]]^2)*USf[4, j1][s1, 2]))/(CW^2*MW^2*SW^2)}}, 
 C[S[1], S[2], S[11, {j2}], -S[11, {j1}]] == 
  {{((I/4)*EL^2*S2A*IndexDelta[j1, j2])/(CW^2*SW^2)}}, 
 C[S[1], S[2], S[12, {s2, j2}], -S[12, {s1, j1}]] == 
  {{((I/4)*EL^2*S2A*IndexDelta[j1, j2]*(Conjugate[USf[2, j1][s2, 1]]*
        (CB^2*MW^2*(-1 + 2*SW^2) + CW^2*Mass[F[2, {j1}]]^2)*
        USf[2, j1][s1, 1] + Conjugate[USf[2, j1][s2, 2]]*
        (-2*CB^2*MW^2*SW^2 + CW^2*Mass[F[2, {j1}]]^2)*USf[2, j1][s1, 2]))/
     (CB^2*CW^2*MW^2*SW^2)}}, C[S[3], S[4], S[11, {j2}], -S[11, {j1}]] == 
  {{((I/4)*EL^2*S2B*IndexDelta[j1, j2])/(CW^2*SW^2)}}, 
 C[S[3], S[4], S[12, {s2, j2}], -S[12, {s1, j1}]] == 
  {{((I/4)*EL^2*S2B*IndexDelta[j1, j2]*(Conjugate[USf[2, j1][s2, 1]]*
        (CB^2*MW^2*(-1 + 2*SW^2) + CW^2*Mass[F[2, {j1}]]^2)*
        USf[2, j1][s1, 1] + Conjugate[USf[2, j1][s2, 2]]*
        (-2*CB^2*MW^2*SW^2 + CW^2*Mass[F[2, {j1}]]^2)*USf[2, j1][s1, 2]))/
     (CB^2*CW^2*MW^2*SW^2)}}, C[S[1], S[2], S[13, {s2, j2, o1}], 
   -S[13, {s1, j1, o2}]] == 
  {{((-I/12)*EL^2*S2A*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      (Conjugate[USf[3, j1][s2, 1]]*(MW^2*SB^2*(-3 + 4*SW^2) + 
         3*CW^2*Mass[F[3, {j1}]]^2)*USf[3, j1][s1, 1] + 
       Conjugate[USf[3, j1][s2, 2]]*(-4*MW^2*SB^2*SW^2 + 
         3*CW^2*Mass[F[3, {j1}]]^2)*USf[3, j1][s1, 2]))/
     (CW^2*MW^2*SB^2*SW^2)}}, C[S[1], S[2], S[14, {s2, j2, o1}], 
   -S[14, {s1, j1, o2}]] == 
  {{((I/12)*EL^2*S2A*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      (Conjugate[USf[4, j1][s2, 1]]*(CB^2*MW^2*(-3 + 2*SW^2) + 
         3*CW^2*Mass[F[4, {j1}]]^2)*USf[4, j1][s1, 1] + 
       Conjugate[USf[4, j1][s2, 2]]*(-2*CB^2*MW^2*SW^2 + 
         3*CW^2*Mass[F[4, {j1}]]^2)*USf[4, j1][s1, 2]))/
     (CB^2*CW^2*MW^2*SW^2)}}, C[S[3], S[4], S[13, {s2, j2, o1}], 
   -S[13, {s1, j1, o2}]] == 
  {{((-I/12)*EL^2*S2B*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      (Conjugate[USf[3, j1][s2, 1]]*(MW^2*SB^2*(-3 + 4*SW^2) + 
         3*CW^2*Mass[F[3, {j1}]]^2)*USf[3, j1][s1, 1] + 
       Conjugate[USf[3, j1][s2, 2]]*(-4*MW^2*SB^2*SW^2 + 
         3*CW^2*Mass[F[3, {j1}]]^2)*USf[3, j1][s1, 2]))/
     (CW^2*MW^2*SB^2*SW^2)}}, C[S[3], S[4], S[14, {s2, j2, o1}], 
   -S[14, {s1, j1, o2}]] == 
  {{((I/12)*EL^2*S2B*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      (Conjugate[USf[4, j1][s2, 1]]*(CB^2*MW^2*(-3 + 2*SW^2) + 
         3*CW^2*Mass[F[4, {j1}]]^2)*USf[4, j1][s1, 1] + 
       Conjugate[USf[4, j1][s2, 2]]*(-2*CB^2*MW^2*SW^2 + 
         3*CW^2*Mass[F[4, {j1}]]^2)*USf[4, j1][s1, 2]))/
     (CB^2*CW^2*MW^2*SW^2)}}, C[S[1], S[5], S[13, {s1, j1, o1}], 
   -S[14, {s2, j2, o2}]] == 
  {{((-I/2)*EL^2*Conjugate[CKM[j1, j2]]*IndexDelta[o1, o2]*
      (S2B*Conjugate[USf[3, j1][s1, 1]]*(-(CA*CB*Mass[F[3, {j1}]]^2) + 
         SB*(CAB*MW^2*SB + SA*TB^2*Mass[F[4, {j2}]]^2))*USf[4, j2][s2, 1] - 
       2*SB^2*SBA*Conjugate[USf[3, j1][s1, 2]]*Mass[F[3, {j1}]]*
        Mass[F[4, {j2}]]*USf[4, j2][s2, 2]))/(Sqrt[2]*MW^2*S2B*SB^2*SW^2)}}, 
 C[S[1], -S[5], S[14, {s2, j2, o1}], -S[13, {s1, j1, o2}]] == 
  {{((-I/2)*EL^2*CKM[j1, j2]*IndexDelta[o1, o2]*
      (S2B*Conjugate[USf[4, j2][s2, 1]]*(-(CA*CB*Mass[F[3, {j1}]]^2) + 
         SB*(CAB*MW^2*SB + SA*TB^2*Mass[F[4, {j2}]]^2))*USf[3, j1][s1, 1] - 
       2*SB^2*SBA*Conjugate[USf[4, j2][s2, 2]]*Mass[F[3, {j1}]]*
        Mass[F[4, {j2}]]*USf[3, j1][s1, 2]))/(Sqrt[2]*MW^2*S2B*SB^2*SW^2)}}, 
 C[S[1], S[6], S[13, {s1, j1, o1}], -S[14, {s2, j2, o2}]] == 
  {{((-I/2)*EL^2*Conjugate[CKM[j1, j2]]*IndexDelta[o1, o2]*
      (S2B*Conjugate[USf[3, j1][s1, 1]]*(CB*MW^2*SAB*SB - 
         CA*CB*Mass[F[3, {j1}]]^2 - SA*SB*Mass[F[4, {j2}]]^2)*
        USf[4, j2][s2, 1] + 2*CB*CBA*SB*Conjugate[USf[3, j1][s1, 2]]*
        Mass[F[3, {j1}]]*Mass[F[4, {j2}]]*USf[4, j2][s2, 2]))/
     (Sqrt[2]*CB*MW^2*S2B*SB*SW^2)}}, 
 C[S[1], -S[6], S[14, {s2, j2, o1}], -S[13, {s1, j1, o2}]] == 
  {{((-I/2)*EL^2*CKM[j1, j2]*IndexDelta[o1, o2]*
      (S2B*Conjugate[USf[4, j2][s2, 1]]*(CB*MW^2*SAB*SB - 
         CA*CB*Mass[F[3, {j1}]]^2 - SA*SB*Mass[F[4, {j2}]]^2)*
        USf[3, j1][s1, 1] + 2*CB*CBA*SB*Conjugate[USf[4, j2][s2, 2]]*
        Mass[F[3, {j1}]]*Mass[F[4, {j2}]]*USf[3, j1][s1, 2]))/
     (Sqrt[2]*CB*MW^2*S2B*SB*SW^2)}}, 
 C[S[3], S[5], S[13, {s1, j1, o1}], -S[14, {s2, j2, o2}]] == 
  {{(EL^2*Conjugate[CKM[j1, j2]]*Conjugate[USf[3, j1][s1, 1]]*
      IndexDelta[o1, o2]*(C2B - Mass[F[3, {j1}]]^2/(MW^2*TB^2) + 
       (TB^2*Mass[F[4, {j2}]]^2)/MW^2)*USf[4, j2][s2, 1])/(2*Sqrt[2]*SW^2)}}, 
 C[S[3], -S[5], S[14, {s2, j2, o1}], -S[13, {s1, j1, o2}]] == 
  {{-(EL^2*CKM[j1, j2]*Conjugate[USf[4, j2][s2, 1]]*IndexDelta[o1, o2]*
       (C2B - Mass[F[3, {j1}]]^2/(MW^2*TB^2) + (TB^2*Mass[F[4, {j2}]]^2)/
         MW^2)*USf[3, j1][s1, 1])/(2*Sqrt[2]*SW^2)}}, 
 C[S[3], S[6], S[13, {s1, j1, o1}], -S[14, {s2, j2, o2}]] == 
  {{(EL^2*Conjugate[CKM[j1, j2]]*IndexDelta[o1, o2]*
      (S2B*Conjugate[USf[3, j1][s1, 1]]*(-Mass[F[3, {j1}]]^2 + 
         TB*(MW^2*S2B - TB*Mass[F[4, {j2}]]^2))*USf[4, j2][s2, 1] - 
       2*TB*Conjugate[USf[3, j1][s1, 2]]*Mass[F[3, {j1}]]*Mass[F[4, {j2}]]*
        USf[4, j2][s2, 2]))/(2*Sqrt[2]*MW^2*S2B*SW^2*TB)}}, 
 C[S[3], -S[6], S[14, {s2, j2, o1}], -S[13, {s1, j1, o2}]] == 
  {{-(EL^2*CKM[j1, j2]*IndexDelta[o1, o2]*(S2B*Conjugate[USf[4, j2][s2, 1]]*
         (-Mass[F[3, {j1}]]^2 + TB*(MW^2*S2B - TB*Mass[F[4, {j2}]]^2))*
         USf[3, j1][s1, 1] - 2*TB*Conjugate[USf[4, j2][s2, 2]]*
         Mass[F[3, {j1}]]*Mass[F[4, {j2}]]*USf[3, j1][s1, 2]))/
     (2*Sqrt[2]*MW^2*S2B*SW^2*TB)}}, 
 C[S[1], S[5], S[11, {j1}], -S[12, {s2, j2}]] == 
  {{((-I/2)*EL^2*IndexDelta[j1, j2]*(CAB + (SA*TB*Mass[F[2, {j1}]]^2)/
        (CB*MW^2))*USf[2, j1][s2, 1])/(Sqrt[2]*SW^2)}}, 
 C[S[1], -S[5], S[12, {s2, j2}], -S[11, {j1}]] == 
  {{((-I/2)*EL^2*Conjugate[USf[2, j1][s2, 1]]*IndexDelta[j1, j2]*
      (CAB + (SA*TB*Mass[F[2, {j1}]]^2)/(CB*MW^2)))/(Sqrt[2]*SW^2)}}, 
 C[S[1], S[6], S[11, {j1}], -S[12, {s2, j2}]] == 
  {{((-I/2)*EL^2*IndexDelta[j1, j2]*(SAB - (SA*Mass[F[2, {j1}]]^2)/(CB*MW^2))*
      USf[2, j1][s2, 1])/(Sqrt[2]*SW^2)}}, 
 C[S[1], -S[6], S[12, {s2, j2}], -S[11, {j1}]] == 
  {{((-I/2)*EL^2*Conjugate[USf[2, j1][s2, 1]]*IndexDelta[j1, j2]*
      (SAB - (SA*Mass[F[2, {j1}]]^2)/(CB*MW^2)))/(Sqrt[2]*SW^2)}}, 
 C[S[3], S[5], S[11, {j1}], -S[12, {s2, j2}]] == 
  {{(EL^2*IndexDelta[j1, j2]*(C2B + (TB^2*Mass[F[2, {j1}]]^2)/MW^2)*
      USf[2, j1][s2, 1])/(2*Sqrt[2]*SW^2)}}, 
 C[S[3], -S[5], S[12, {s2, j2}], -S[11, {j1}]] == 
  {{-(EL^2*Conjugate[USf[2, j1][s2, 1]]*IndexDelta[j1, j2]*
       (C2B + (TB^2*Mass[F[2, {j1}]]^2)/MW^2))/(2*Sqrt[2]*SW^2)}}, 
 C[S[3], S[6], S[11, {j1}], -S[12, {s2, j2}]] == 
  {{(EL^2*IndexDelta[j1, j2]*(S2B - (TB*Mass[F[2, {j1}]]^2)/MW^2)*
      USf[2, j1][s2, 1])/(2*Sqrt[2]*SW^2)}}, 
 C[S[3], -S[6], S[12, {s2, j2}], -S[11, {j1}]] == 
  {{-(EL^2*Conjugate[USf[2, j1][s2, 1]]*IndexDelta[j1, j2]*
       (S2B - (TB*Mass[F[2, {j1}]]^2)/MW^2))/(2*Sqrt[2]*SW^2)}}, 
 C[S[2], S[5], S[13, {s1, j1, o1}], -S[14, {s2, j2, o2}]] == 
  {{((-I/2)*EL^2*Conjugate[CKM[j1, j2]]*IndexDelta[o1, o2]*
      (-(S2B*Conjugate[USf[3, j1][s1, 1]]*(CB*SA*Mass[F[3, {j1}]]^2 + 
          SB*(-(MW^2*SAB*SB) + CA*TB^2*Mass[F[4, {j2}]]^2))*
         USf[4, j2][s2, 1]) - 2*CBA*SB^2*Conjugate[USf[3, j1][s1, 2]]*
        Mass[F[3, {j1}]]*Mass[F[4, {j2}]]*USf[4, j2][s2, 2]))/
     (Sqrt[2]*MW^2*S2B*SB^2*SW^2)}}, 
 C[S[2], -S[5], S[14, {s2, j2, o1}], -S[13, {s1, j1, o2}]] == 
  {{((-I/2)*EL^2*CKM[j1, j2]*IndexDelta[o1, o2]*
      (-(S2B*Conjugate[USf[4, j2][s2, 1]]*(CB*SA*Mass[F[3, {j1}]]^2 + 
          SB*(-(MW^2*SAB*SB) + CA*TB^2*Mass[F[4, {j2}]]^2))*
         USf[3, j1][s1, 1]) - 2*CBA*SB^2*Conjugate[USf[4, j2][s2, 2]]*
        Mass[F[3, {j1}]]*Mass[F[4, {j2}]]*USf[3, j1][s1, 2]))/
     (Sqrt[2]*MW^2*S2B*SB^2*SW^2)}}, 
 C[S[2], S[6], S[13, {s1, j1, o1}], -S[14, {s2, j2, o2}]] == 
  {{((I/2)*EL^2*Conjugate[CKM[j1, j2]]*IndexDelta[o1, o2]*
      (S2B*Conjugate[USf[3, j1][s1, 1]]*(CAB*CB*MW^2*SB + 
         CB*SA*Mass[F[3, {j1}]]^2 - CA*SB*Mass[F[4, {j2}]]^2)*
        USf[4, j2][s2, 1] + 2*CB*SB*SBA*Conjugate[USf[3, j1][s1, 2]]*
        Mass[F[3, {j1}]]*Mass[F[4, {j2}]]*USf[4, j2][s2, 2]))/
     (Sqrt[2]*CB*MW^2*S2B*SB*SW^2)}}, 
 C[S[2], -S[6], S[14, {s2, j2, o1}], -S[13, {s1, j1, o2}]] == 
  {{((I/2)*EL^2*CKM[j1, j2]*IndexDelta[o1, o2]*
      (S2B*Conjugate[USf[4, j2][s2, 1]]*(CAB*CB*MW^2*SB + 
         CB*SA*Mass[F[3, {j1}]]^2 - CA*SB*Mass[F[4, {j2}]]^2)*
        USf[3, j1][s1, 1] + 2*CB*SB*SBA*Conjugate[USf[4, j2][s2, 2]]*
        Mass[F[3, {j1}]]*Mass[F[4, {j2}]]*USf[3, j1][s1, 2]))/
     (Sqrt[2]*CB*MW^2*S2B*SB*SW^2)}}, 
 C[S[4], S[5], S[13, {s1, j1, o1}], -S[14, {s2, j2, o2}]] == 
  {{(EL^2*Conjugate[CKM[j1, j2]]*IndexDelta[o1, o2]*
      (S2B*Conjugate[USf[3, j1][s1, 1]]*(-Mass[F[3, {j1}]]^2 + 
         TB*(MW^2*S2B - TB*Mass[F[4, {j2}]]^2))*USf[4, j2][s2, 1] + 
       2*TB*Conjugate[USf[3, j1][s1, 2]]*Mass[F[3, {j1}]]*Mass[F[4, {j2}]]*
        USf[4, j2][s2, 2]))/(2*Sqrt[2]*MW^2*S2B*SW^2*TB)}}, 
 C[S[4], -S[5], S[14, {s2, j2, o1}], -S[13, {s1, j1, o2}]] == 
  {{-(EL^2*CKM[j1, j2]*IndexDelta[o1, o2]*(S2B*Conjugate[USf[4, j2][s2, 1]]*
         (-Mass[F[3, {j1}]]^2 + TB*(MW^2*S2B - TB*Mass[F[4, {j2}]]^2))*
         USf[3, j1][s1, 1] + 2*TB*Conjugate[USf[4, j2][s2, 2]]*
         Mass[F[3, {j1}]]*Mass[F[4, {j2}]]*USf[3, j1][s1, 2]))/
     (2*Sqrt[2]*MW^2*S2B*SW^2*TB)}}, 
 C[S[4], S[6], S[13, {s1, j1, o1}], -S[14, {s2, j2, o2}]] == 
  {{-(EL^2*Conjugate[CKM[j1, j2]]*Conjugate[USf[3, j1][s1, 1]]*
       IndexDelta[o1, o2]*(C2B*MW^2 + Mass[F[3, {j1}]]^2 - 
        Mass[F[4, {j2}]]^2)*USf[4, j2][s2, 1])/(2*Sqrt[2]*MW^2*SW^2)}}, 
 C[S[4], -S[6], S[14, {s2, j2, o1}], -S[13, {s1, j1, o2}]] == 
  {{(EL^2*CKM[j1, j2]*Conjugate[USf[4, j2][s2, 1]]*IndexDelta[o1, o2]*
      (C2B*MW^2 + Mass[F[3, {j1}]]^2 - Mass[F[4, {j2}]]^2)*USf[3, j1][s1, 1])/
     (2*Sqrt[2]*MW^2*SW^2)}}, C[S[2], S[5], S[11, {j1}], -S[12, {s2, j2}]] == 
  {{((-I/2)*EL^2*IndexDelta[j1, j2]*(SAB - (CA*TB*Mass[F[2, {j1}]]^2)/
        (CB*MW^2))*USf[2, j1][s2, 1])/(Sqrt[2]*SW^2)}}, 
 C[S[2], -S[5], S[12, {s2, j2}], -S[11, {j1}]] == 
  {{((-I/2)*EL^2*Conjugate[USf[2, j1][s2, 1]]*IndexDelta[j1, j2]*
      (SAB - (CA*TB*Mass[F[2, {j1}]]^2)/(CB*MW^2)))/(Sqrt[2]*SW^2)}}, 
 C[S[2], S[6], S[11, {j1}], -S[12, {s2, j2}]] == 
  {{((-I/2)*EL^2*IndexDelta[j1, j2]*(-CAB + (CA*Mass[F[2, {j1}]]^2)/
        (CB*MW^2))*USf[2, j1][s2, 1])/(Sqrt[2]*SW^2)}}, 
 C[S[2], -S[6], S[12, {s2, j2}], -S[11, {j1}]] == 
  {{((-I/2)*EL^2*Conjugate[USf[2, j1][s2, 1]]*IndexDelta[j1, j2]*
      (-CAB + (CA*Mass[F[2, {j1}]]^2)/(CB*MW^2)))/(Sqrt[2]*SW^2)}}, 
 C[S[4], S[5], S[11, {j1}], -S[12, {s2, j2}]] == 
  {{(EL^2*IndexDelta[j1, j2]*(S2B - (TB*Mass[F[2, {j1}]]^2)/MW^2)*
      USf[2, j1][s2, 1])/(2*Sqrt[2]*SW^2)}}, 
 C[S[4], -S[5], S[12, {s2, j2}], -S[11, {j1}]] == 
  {{-(EL^2*Conjugate[USf[2, j1][s2, 1]]*IndexDelta[j1, j2]*
       (S2B - (TB*Mass[F[2, {j1}]]^2)/MW^2))/(2*Sqrt[2]*SW^2)}}, 
 C[S[4], S[6], S[11, {j1}], -S[12, {s2, j2}]] == 
  {{(EL^2*IndexDelta[j1, j2]*(-C2B + Mass[F[2, {j1}]]^2/MW^2)*
      USf[2, j1][s2, 1])/(2*Sqrt[2]*SW^2)}}, 
 C[S[4], -S[6], S[12, {s2, j2}], -S[11, {j1}]] == 
  {{-(EL^2*Conjugate[USf[2, j1][s2, 1]]*IndexDelta[j1, j2]*
       (-C2B + Mass[F[2, {j1}]]^2/MW^2))/(2*Sqrt[2]*SW^2)}}, 
 C[S[5], -S[5], S[11, {j1}], -S[11, {j2}]] == 
  {{((I/2)*EL^2*IndexDelta[j1, j2]*((C2B*(-2 + CW^(-2)))/2 - 
       (TB^2*Mass[F[2, {j1}]]^2)/MW^2))/SW^2}}, 
 C[S[5], -S[6], S[11, {j1}], -S[11, {j2}]] == 
  {{((I/2)*EL^2*IndexDelta[j1, j2]*(((-2 + CW^(-2))*S2B)/2 + 
       (TB*Mass[F[2, {j1}]]^2)/MW^2))/SW^2}}, 
 C[S[6], -S[5], S[11, {j1}], -S[11, {j2}]] == 
  {{((I/2)*EL^2*IndexDelta[j1, j2]*(((-2 + CW^(-2))*S2B)/2 + 
       (TB*Mass[F[2, {j1}]]^2)/MW^2))/SW^2}}, 
 C[S[5], -S[5], S[12, {s1, j1}], -S[12, {s2, j2}]] == 
  {{((I/4)*EL^2*IndexDelta[j1, j2]*(C2B*MW^2*Conjugate[USf[2, j1][s1, 1]]*
        USf[2, j1][s2, 1] - 2*Conjugate[USf[2, j1][s1, 2]]*
        (C2B*MW^2*SW^2 + CW^2*TB^2*Mass[F[2, {j1}]]^2)*USf[2, j1][s2, 2]))/
     (CW^2*MW^2*SW^2)}}, C[S[5], -S[6], S[12, {s1, j1}], -S[12, {s2, j2}]] == 
  {{((I/2)*EL^2*IndexDelta[j1, j2]*(S2B*(1 + (-1/2 + SW^2)/CW^2)*
        Conjugate[USf[2, j1][s1, 1]]*USf[2, j1][s2, 1] + 
       Conjugate[USf[2, j1][s1, 2]]*(-((S2B*SW^2)/CW^2) + 
         (TB*Mass[F[2, {j1}]]^2)/MW^2)*USf[2, j1][s2, 2]))/SW^2}}, 
 C[S[6], -S[5], S[12, {s1, j1}], -S[12, {s2, j2}]] == 
  {{((I/2)*EL^2*IndexDelta[j1, j2]*(S2B*(1 + (-1/2 + SW^2)/CW^2)*
        Conjugate[USf[2, j1][s1, 1]]*USf[2, j1][s2, 1] + 
       Conjugate[USf[2, j1][s1, 2]]*(-((S2B*SW^2)/CW^2) + 
         (TB*Mass[F[2, {j1}]]^2)/MW^2)*USf[2, j1][s2, 2]))/SW^2}}, 
 C[S[5], -S[5], S[13, {s1, j1, o1}], -S[13, {s2, j2, o2}]] == 
  {{((-I/12)*EL^2*IndexDelta[o1, o2]*(TB^2*Conjugate[USf[3, j1][s1, 1]]*
        (C2B*(1 + 2*CW^2)*MW^2*IndexDelta[j1, j2] + 6*CW^2*TB^2*
          IndexSum[CKM[j2, gn]*Conjugate[CKM[j1, gn]]*Mass[F[4, {gn}]]^2, 
           {gn, 3}])*USf[3, j2][s2, 1] + 2*Conjugate[USf[3, j1][s1, 2]]*
        IndexDelta[j1, j2]*(-2*C2B*MW^2*SW^2*TB^2 + 
         3*CW^2*Mass[F[3, {j1}]]^2)*USf[3, j2][s2, 2]))/
     (CW^2*MW^2*SW^2*TB^2)}}, C[S[5], -S[6], S[13, {s1, j1, o1}], 
   -S[13, {s2, j2, o2}]] == 
  {{((-I/12)*EL^2*IndexDelta[o1, o2]*(TB*Conjugate[USf[3, j1][s1, 1]]*
        ((1 + 2*CW^2)*MW^2*S2B*IndexDelta[j1, j2] - 
         6*CW^2*TB*IndexSum[CKM[j2, gn]*Conjugate[CKM[j1, gn]]*
            Mass[F[4, {gn}]]^2, {gn, 3}])*USf[3, j2][s2, 1] + 
       2*Conjugate[USf[3, j1][s1, 2]]*IndexDelta[j1, j2]*
        (-2*MW^2*S2B*SW^2*TB + 3*CW^2*Mass[F[3, {j1}]]^2)*USf[3, j2][s2, 2]))/
     (CW^2*MW^2*SW^2*TB)}}, C[S[6], -S[5], S[13, {s1, j1, o1}], 
   -S[13, {s2, j2, o2}]] == 
  {{((-I/12)*EL^2*IndexDelta[o1, o2]*(TB*Conjugate[USf[3, j1][s1, 1]]*
        ((1 + 2*CW^2)*MW^2*S2B*IndexDelta[j1, j2] - 
         6*CW^2*TB*IndexSum[CKM[j2, gn]*Conjugate[CKM[j1, gn]]*
            Mass[F[4, {gn}]]^2, {gn, 3}])*USf[3, j2][s2, 1] + 
       2*Conjugate[USf[3, j1][s1, 2]]*IndexDelta[j1, j2]*
        (-2*MW^2*S2B*SW^2*TB + 3*CW^2*Mass[F[3, {j1}]]^2)*USf[3, j2][s2, 2]))/
     (CW^2*MW^2*SW^2*TB)}}, C[S[5], -S[5], S[14, {s1, j1, o1}], 
   -S[14, {s2, j2, o2}]] == 
  {{((I/12)*EL^2*IndexDelta[o1, o2]*(Conjugate[USf[4, j1][s1, 1]]*
        (C2B*(-1 + 4*CW^2)*MW^2*TB^2*IndexDelta[j1, j2] - 
         6*CW^2*IndexSum[CKM[gn, j1]*Conjugate[CKM[gn, j2]]*
            Mass[F[3, {gn}]]^2, {gn, 3}])*USf[4, j2][s2, 1] - 
       2*TB^2*Conjugate[USf[4, j1][s1, 2]]*IndexDelta[j1, j2]*
        (C2B*MW^2*SW^2 + 3*CW^2*TB^2*Mass[F[4, {j1}]]^2)*USf[4, j2][s2, 2]))/
     (CW^2*MW^2*SW^2*TB^2)}}, C[S[5], -S[6], S[14, {s1, j1, o1}], 
   -S[14, {s2, j2, o2}]] == 
  {{((I/12)*EL^2*IndexDelta[o1, o2]*(Conjugate[USf[4, j1][s1, 1]]*
        ((-1 + 4*CW^2)*MW^2*S2B*TB*IndexDelta[j1, j2] - 
         6*CW^2*IndexSum[CKM[gn, j1]*Conjugate[CKM[gn, j2]]*
            Mass[F[3, {gn}]]^2, {gn, 3}])*USf[4, j2][s2, 1] + 
       2*TB*Conjugate[USf[4, j1][s1, 2]]*IndexDelta[j1, j2]*
        (-(MW^2*S2B*SW^2) + 3*CW^2*TB*Mass[F[4, {j1}]]^2)*USf[4, j2][s2, 2]))/
     (CW^2*MW^2*SW^2*TB)}}, C[S[6], -S[5], S[14, {s1, j1, o1}], 
   -S[14, {s2, j2, o2}]] == 
  {{((I/12)*EL^2*IndexDelta[o1, o2]*(Conjugate[USf[4, j1][s1, 1]]*
        ((-1 + 4*CW^2)*MW^2*S2B*TB*IndexDelta[j1, j2] - 
         6*CW^2*IndexSum[CKM[gn, j1]*Conjugate[CKM[gn, j2]]*
            Mass[F[3, {gn}]]^2, {gn, 3}])*USf[4, j2][s2, 1] + 
       2*TB*Conjugate[USf[4, j1][s1, 2]]*IndexDelta[j1, j2]*
        (-(MW^2*S2B*SW^2) + 3*CW^2*TB*Mass[F[4, {j1}]]^2)*USf[4, j2][s2, 2]))/
     (CW^2*MW^2*SW^2*TB)}}, C[S[6], -S[6], S[11, {j1}], -S[11, {j2}]] == 
  {{((I/4)*EL^2*IndexDelta[j1, j2]*(C2B*(-1 + 2*CW^2)*MW^2 - 
       2*CW^2*Mass[F[2, {j1}]]^2))/(CW^2*MW^2*SW^2)}}, 
 C[S[6], -S[6], S[12, {s1, j1}], -S[12, {s2, j2}]] == 
  {{((I/2)*EL^2*IndexDelta[j1, j2]*(-(C2B*(1 + (-1/2 + SW^2)/CW^2)*
         Conjugate[USf[2, j1][s1, 1]]*USf[2, j1][s2, 1]) + 
       Conjugate[USf[2, j1][s1, 2]]*((C2B*SW^2)/CW^2 - Mass[F[2, {j1}]]^2/
          MW^2)*USf[2, j1][s2, 2]))/SW^2}}, 
 C[S[6], -S[6], S[13, {s1, j1, o1}], -S[13, {s2, j2, o2}]] == 
  {{((I/12)*EL^2*IndexDelta[o1, o2]*(Conjugate[USf[3, j1][s1, 1]]*
        (C2B*(1 + 2*CW^2)*MW^2*IndexDelta[j1, j2] - 
         6*CW^2*IndexSum[CKM[j2, gn]*Conjugate[CKM[j1, gn]]*
            Mass[F[4, {gn}]]^2, {gn, 3}])*USf[3, j2][s2, 1] - 
       2*Conjugate[USf[3, j1][s1, 2]]*IndexDelta[j1, j2]*
        (2*C2B*MW^2*SW^2 + 3*CW^2*Mass[F[3, {j1}]]^2)*USf[3, j2][s2, 2]))/
     (CW^2*MW^2*SW^2)}}, C[S[6], -S[6], S[14, {s1, j1, o1}], 
   -S[14, {s2, j2, o2}]] == 
  {{((-I/12)*EL^2*IndexDelta[o1, o2]*(Conjugate[USf[4, j1][s1, 1]]*
        (C2B*(-1 + 4*CW^2)*MW^2*IndexDelta[j1, j2] + 
         6*CW^2*IndexSum[CKM[gn, j1]*Conjugate[CKM[gn, j2]]*
            Mass[F[3, {gn}]]^2, {gn, 3}])*USf[4, j2][s2, 1] - 
       2*Conjugate[USf[4, j1][s1, 2]]*IndexDelta[j1, j2]*
        (C2B*MW^2*SW^2 - 3*CW^2*Mass[F[4, {j1}]]^2)*USf[4, j2][s2, 2]))/
     (CW^2*MW^2*SW^2)}}, C[S[11, {j1}], -S[11, {j2}], V[2], V[2]] == 
  {{((I/2)*EL^2*IndexDelta[j1, j2])/(CW^2*SW^2)}}, 
 C[S[12, {s1, j1}], -S[12, {s2, j2}], V[1], V[1]] == 
  {{(2*I)*EL^2*IndexDelta[j1, j2]*IndexDelta[s1, s2]}}, 
 C[S[12, {s1, j1}], -S[12, {s2, j2}], V[1], V[2]] == 
  {{((-I)*EL^2*$HKSign*IndexDelta[j1, j2]*
      ((-1 + 2*SW^2)*Conjugate[USf[2, j1][s1, 1]]*USf[2, j1][s2, 1] + 
       2*SW^2*Conjugate[USf[2, j1][s1, 2]]*USf[2, j1][s2, 2]))/(CW*SW)}}, 
 C[S[12, {s1, j1}], -S[12, {s2, j2}], V[2], V[2]] == 
  {{((I/2)*EL^2*IndexDelta[j1, j2]*
      ((1 - 2*SW^2)^2*Conjugate[USf[2, j1][s1, 1]]*USf[2, j1][s2, 1] + 
       4*SW^4*Conjugate[USf[2, j1][s1, 2]]*USf[2, j1][s2, 2]))/(CW^2*SW^2)}}, 
 C[S[13, {s1, j1, o1}], -S[13, {s2, j2, o2}], V[1], V[1]] == 
  {{((8*I)/9)*EL^2*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
     IndexDelta[s1, s2]}}, C[S[13, {s1, j1, o1}], -S[13, {s2, j2, o2}], V[1], 
   V[2]] == {{(((-2*I)/9)*EL^2*$HKSign*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      ((-3 + 4*SW^2)*Conjugate[USf[3, j1][s1, 1]]*USf[3, j1][s2, 1] + 
       4*SW^2*Conjugate[USf[3, j1][s1, 2]]*USf[3, j1][s2, 2]))/(CW*SW)}}, 
 C[S[13, {s1, j1, o1}], -S[13, {s2, j2, o2}], V[2], V[2]] == 
  {{((I/18)*EL^2*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      ((3 - 4*SW^2)^2*Conjugate[USf[3, j1][s1, 1]]*USf[3, j1][s2, 1] + 
       16*SW^4*Conjugate[USf[3, j1][s1, 2]]*USf[3, j1][s2, 2]))/
     (CW^2*SW^2)}}, C[S[14, {s1, j1, o1}], -S[14, {s2, j2, o2}], V[1], 
   V[1]] == {{((2*I)/9)*EL^2*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
     IndexDelta[s1, s2]}}, C[S[14, {s1, j1, o1}], -S[14, {s2, j2, o2}], V[1], 
   V[2]] == {{((-I/9)*EL^2*$HKSign*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      ((-3 + 2*SW^2)*Conjugate[USf[4, j1][s1, 1]]*USf[4, j1][s2, 1] + 
       2*SW^2*Conjugate[USf[4, j1][s1, 2]]*USf[4, j1][s2, 2]))/(CW*SW)}}, 
 C[S[14, {s1, j1, o1}], -S[14, {s2, j2, o2}], V[2], V[2]] == 
  {{((I/18)*EL^2*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      ((3 - 2*SW^2)^2*Conjugate[USf[4, j1][s1, 1]]*USf[4, j1][s2, 1] + 
       4*SW^4*Conjugate[USf[4, j1][s1, 2]]*USf[4, j1][s2, 2]))/(CW^2*SW^2)}}, 
 C[S[13, {s1, j1, o1}], -S[14, {s2, j2, o2}], V[1], V[3]] == 
  {{((I/3)*EL^2*$HKSign*Conjugate[CKM[j1, j2]]*Conjugate[USf[3, j1][s1, 1]]*
      IndexDelta[o1, o2]*USf[4, j2][s2, 1])/(Sqrt[2]*SW)}}, 
 C[S[14, {s2, j2, o1}], -S[13, {s1, j1, o2}], V[1], -V[3]] == 
  {{((I/3)*EL^2*$HKSign*CKM[j1, j2]*Conjugate[USf[4, j2][s2, 1]]*
      IndexDelta[o1, o2]*USf[3, j1][s1, 1])/(Sqrt[2]*SW)}}, 
 C[S[11, {j1}], -S[12, {s2, j2}], V[1], V[3]] == 
  {{((-I)*EL^2*$HKSign*IndexDelta[j1, j2]*USf[2, j1][s2, 1])/(Sqrt[2]*SW)}}, 
 C[S[12, {s2, j2}], -S[11, {j1}], V[1], -V[3]] == 
  {{((-I)*EL^2*$HKSign*Conjugate[USf[2, j1][s2, 1]]*IndexDelta[j1, j2])/
     (Sqrt[2]*SW)}}, C[S[13, {s1, j1, o1}], -S[14, {s2, j2, o2}], V[2], 
   V[3]] == {{((-I/3)*EL^2*Conjugate[CKM[j1, j2]]*
      Conjugate[USf[3, j1][s1, 1]]*IndexDelta[o1, o2]*USf[4, j2][s2, 1])/
     (Sqrt[2]*CW)}}, C[S[14, {s2, j2, o1}], -S[13, {s1, j1, o2}], V[2], 
   -V[3]] == {{((-I/3)*EL^2*CKM[j1, j2]*Conjugate[USf[4, j2][s2, 1]]*
      IndexDelta[o1, o2]*USf[3, j1][s1, 1])/(Sqrt[2]*CW)}}, 
 C[S[11, {j1}], -S[12, {s2, j2}], V[2], V[3]] == 
  {{(I*EL^2*IndexDelta[j1, j2]*USf[2, j1][s2, 1])/(Sqrt[2]*CW)}}, 
 C[S[12, {s2, j2}], -S[11, {j1}], V[2], -V[3]] == 
  {{(I*EL^2*Conjugate[USf[2, j1][s2, 1]]*IndexDelta[j1, j2])/(Sqrt[2]*CW)}}, 
 C[S[11, {j1}], -S[11, {j2}], V[3], -V[3]] == 
  {{((I/2)*EL^2*IndexDelta[j1, j2])/SW^2}}, 
 C[S[12, {s1, j1}], -S[12, {s2, j2}], V[3], -V[3]] == 
  {{((I/2)*EL^2*Conjugate[USf[2, j1][s1, 1]]*IndexDelta[j1, j2]*
      USf[2, j1][s2, 1])/SW^2}}, C[S[13, {s1, j1, o1}], -S[13, {s2, j2, o2}], 
   V[3], -V[3]] == {{((I/2)*EL^2*Conjugate[USf[3, j1][s1, 1]]*
      IndexDelta[j1, j2]*IndexDelta[o1, o2]*USf[3, j1][s2, 1])/SW^2}}, 
 C[S[14, {s1, j1, o1}], -S[14, {s2, j2, o2}], V[3], -V[3]] == 
  {{((I/2)*EL^2*Conjugate[USf[4, j1][s1, 1]]*IndexDelta[j1, j2]*
      IndexDelta[o1, o2]*USf[4, j1][s2, 1])/SW^2}}, 
 C[S[14, {s1, j1, o1}], -S[14, {s2, j2, o2}], S[14, {s3, j3, o3}], 
   -S[14, {s4, j4, o4}]] == 
  {{IndexDelta[j1, j4]*IndexDelta[j2, j3]*((-I)*GS^2*SUNTSum[o2, o3, o4, o1]*
        (Conjugate[USf[4, j1][s1, 1]]*USf[4, j1][s4, 1] - 
         Conjugate[USf[4, j1][s1, 2]]*USf[4, j1][s4, 2])*
        (Conjugate[USf[4, j2][s3, 1]]*USf[4, j2][s2, 1] - 
         Conjugate[USf[4, j2][s3, 2]]*USf[4, j2][s2, 2]) - 
       ((I/36)*EL^2*IndexDelta[o1, o4]*IndexDelta[o2, o3]*
         (Conjugate[USf[4, j1][s1, 1]]*(CB^2*(1 + 8*CW^2)*MW^2*
             Conjugate[USf[4, j2][s3, 1]]*USf[4, j1][s4, 1]*
             USf[4, j2][s2, 1] + 2*Conjugate[USf[4, j2][s3, 2]]*
             (9*CW^2*MQD[j1]*MQD[j2]*USf[4, j1][s4, 2]*USf[4, j2][s2, 1] + 
              CB^2*MW^2*SW^2*USf[4, j1][s4, 1]*USf[4, j2][s2, 2])) + 
          2*Conjugate[USf[4, j1][s1, 2]]*(2*CB^2*MW^2*SW^2*
             Conjugate[USf[4, j2][s3, 2]]*USf[4, j1][s4, 2]*
             USf[4, j2][s2, 2] + Conjugate[USf[4, j2][s3, 1]]*
             (CB^2*MW^2*SW^2*USf[4, j1][s4, 2]*USf[4, j2][s2, 1] + 
              9*CW^2*MQD[j1]*MQD[j2]*USf[4, j1][s4, 1]*USf[4, j2][s2, 2]))))/
        (CB^2*CW^2*MW^2*SW^2)) + IndexDelta[j1, j2]*IndexDelta[j3, j4]*
      ((-I)*GS^2*SUNTSum[o2, o1, o4, o3]*(Conjugate[USf[4, j1][s1, 1]]*
          USf[4, j1][s2, 1] - Conjugate[USf[4, j1][s1, 2]]*USf[4, j1][s2, 2])*
        (Conjugate[USf[4, j3][s3, 1]]*USf[4, j3][s4, 1] - 
         Conjugate[USf[4, j3][s3, 2]]*USf[4, j3][s4, 2]) - 
       ((I/36)*EL^2*IndexDelta[o1, o2]*IndexDelta[o3, o4]*
         (Conjugate[USf[4, j1][s1, 1]]*(CB^2*(1 + 8*CW^2)*MW^2*
             Conjugate[USf[4, j3][s3, 1]]*USf[4, j1][s2, 1]*
             USf[4, j3][s4, 1] + 2*Conjugate[USf[4, j3][s3, 2]]*
             (9*CW^2*MQD[j1]*MQD[j3]*USf[4, j1][s2, 2]*USf[4, j3][s4, 1] + 
              CB^2*MW^2*SW^2*USf[4, j1][s2, 1]*USf[4, j3][s4, 2])) + 
          2*Conjugate[USf[4, j1][s1, 2]]*(2*CB^2*MW^2*SW^2*
             Conjugate[USf[4, j3][s3, 2]]*USf[4, j1][s2, 2]*
             USf[4, j3][s4, 2] + Conjugate[USf[4, j3][s3, 1]]*
             (CB^2*MW^2*SW^2*USf[4, j1][s2, 2]*USf[4, j3][s4, 1] + 
              9*CW^2*MQD[j1]*MQD[j3]*USf[4, j1][s2, 1]*USf[4, j3][s4, 2]))))/
        (CB^2*CW^2*MW^2*SW^2))}}, 
 C[S[14, {s1, j1, o1}], -S[14, {s2, j2, o2}], S[12, {s3, j3}], 
   -S[12, {s4, j4}]] == 
  {{((-I/12)*EL^2*IndexDelta[j1, j2]*IndexDelta[j3, j4]*IndexDelta[o1, o2]*
      (Conjugate[USf[2, j3][s3, 1]]*(CB^2*MW^2*(3*CW^2 - SW^2)*
          Conjugate[USf[4, j1][s1, 1]]*USf[2, j3][s4, 1]*USf[4, j1][s2, 1] + 
         2*Conjugate[USf[4, j1][s1, 2]]*(3*CW^2*MLE[j3]*MQD[j1]*
            USf[2, j3][s4, 2]*USf[4, j1][s2, 1] - CB^2*MW^2*SW^2*
            USf[2, j3][s4, 1]*USf[4, j1][s2, 2])) + 
       2*Conjugate[USf[2, j3][s3, 2]]*(2*CB^2*MW^2*SW^2*
          Conjugate[USf[4, j1][s1, 2]]*USf[2, j3][s4, 2]*USf[4, j1][s2, 2] + 
         Conjugate[USf[4, j1][s1, 1]]*(CB^2*MW^2*SW^2*USf[2, j3][s4, 2]*
            USf[4, j1][s2, 1] + 3*CW^2*MLE[j3]*MQD[j1]*USf[2, j3][s4, 1]*
            USf[4, j1][s2, 2]))))/(CB^2*CW^2*MW^2*SW^2)}}, 
 C[S[14, {s1, j1, o1}], -S[14, {s2, j2, o2}], S[11, {j3}], -S[11, {j4}]] == 
  {{((I/12)*EL^2*IndexDelta[j1, j2]*IndexDelta[j3, j4]*IndexDelta[o1, o2]*
      ((1 + 2*CW^2)*Conjugate[USf[4, j1][s1, 1]]*USf[4, j1][s2, 1] + 
       2*SW^2*Conjugate[USf[4, j1][s1, 2]]*USf[4, j1][s2, 2]))/(CW^2*SW^2)}}, 
 C[S[14, {s1, j1, o1}], -S[14, {s2, j2, o2}], S[13, {s3, j3, o3}], 
   -S[13, {s4, j4, o4}]] == 
  {{IndexDelta[j1, j2]*IndexDelta[j3, j4]*((-I)*GS^2*SUNTSum[o2, o1, o4, o3]*
        (Conjugate[USf[3, j3][s3, 1]]*USf[3, j3][s4, 1] - 
         Conjugate[USf[3, j3][s3, 2]]*USf[3, j3][s4, 2])*
        (Conjugate[USf[4, j1][s1, 1]]*USf[4, j1][s2, 1] - 
         Conjugate[USf[4, j1][s1, 2]]*USf[4, j1][s2, 2]) + 
       ((I/36)*EL^2*IndexDelta[o1, o2]*IndexDelta[o3, o4]*
         (4*SW^2*Conjugate[USf[3, j3][s3, 2]]*USf[3, j3][s4, 2]*
           (Conjugate[USf[4, j1][s1, 1]]*USf[4, j1][s2, 1] + 
            2*Conjugate[USf[4, j1][s1, 2]]*USf[4, j1][s2, 2]) + 
          Conjugate[USf[3, j3][s3, 1]]*USf[3, j3][s4, 1]*
           ((9*CW^2 - SW^2)*Conjugate[USf[4, j1][s1, 1]]*USf[4, j1][s2, 1] - 
            2*SW^2*Conjugate[USf[4, j1][s1, 2]]*USf[4, j1][s2, 2])))/
        (CW^2*SW^2)) - ((I/2)*EL^2*CKM[j4, j1]*Conjugate[CKM[j3, j2]]*
       IndexDelta[o1, o4]*IndexDelta[o2, o3]*
       (CB^2*Conjugate[USf[3, j3][s3, 2]]*Conjugate[USf[4, j1][s1, 1]]*
         MQU[j3]*MQU[j4]*USf[3, j4][s4, 2]*USf[4, j2][s2, 1] + 
        SB^2*Conjugate[USf[3, j3][s3, 1]]*USf[3, j4][s4, 1]*
         (CB^2*MW^2*Conjugate[USf[4, j1][s1, 1]]*USf[4, j2][s2, 1] + 
          Conjugate[USf[4, j1][s1, 2]]*MQD[j1]*MQD[j2]*USf[4, j2][s2, 2])))/
      (CB^2*MW^2*SB^2*SW^2)}}, C[S[14, {s1, j1, o1}], -S[12, {s2, j2}], 
   S[11, {j3}], -S[13, {s4, j4, o4}]] == 
  {{((-I/2)*EL^2*CKM[j4, j1]*IndexDelta[j2, j3]*IndexDelta[o1, o4]*
      (CB^2*MW^2*Conjugate[USf[4, j1][s1, 1]]*USf[2, j2][s2, 1] + 
       Conjugate[USf[4, j1][s1, 2]]*MLE[j2]*MQD[j1]*USf[2, j2][s2, 2])*
      USf[3, j4][s4, 1])/(CB^2*MW^2*SW^2)}}, 
 C[S[12, {s1, j1}], -S[14, {s2, j2, o2}], S[13, {s3, j3, o3}], 
   -S[11, {j4}]] == 
  {{((-I/2)*EL^2*Conjugate[CKM[j3, j2]]*Conjugate[USf[3, j3][s3, 1]]*
      IndexDelta[j1, j4]*IndexDelta[o2, o3]*
      (CB^2*MW^2*Conjugate[USf[2, j1][s1, 1]]*USf[4, j2][s2, 1] + 
       Conjugate[USf[2, j1][s1, 2]]*MLE[j1]*MQD[j2]*USf[4, j2][s2, 2]))/
     (CB^2*MW^2*SW^2)}}, C[S[12, {s1, j1}], -S[12, {s2, j2}], 
   S[12, {s3, j3}], -S[12, {s4, j4}]] == 
  {{((-I/4)*EL^2*(Conjugate[USf[2, j1][s1, 1]]*
        (CB^2*MW^2*Conjugate[USf[2, j2][s3, 1]]*IndexDelta[j1, j4]*
          IndexDelta[j2, j3]*USf[2, j1][s4, 1]*USf[2, j2][s2, 1] + 
         2*Conjugate[USf[2, j2][s3, 2]]*IndexDelta[j1, j4]*IndexDelta[j2, j3]*
          (CW^2*MLE[j1]*MLE[j2]*USf[2, j1][s4, 2]*USf[2, j2][s2, 1] - 
           CB^2*MW^2*SW^2*USf[2, j1][s4, 1]*USf[2, j2][s2, 2]) + 
         IndexDelta[j1, j2]*IndexDelta[j3, j4]*
          (CB^2*MW^2*Conjugate[USf[2, j3][s3, 1]]*USf[2, j1][s2, 1]*
            USf[2, j3][s4, 1] + 2*Conjugate[USf[2, j3][s3, 2]]*
            (CW^2*MLE[j1]*MLE[j3]*USf[2, j1][s2, 2]*USf[2, j3][s4, 1] - 
             CB^2*MW^2*SW^2*USf[2, j1][s2, 1]*USf[2, j3][s4, 2]))) + 
       2*Conjugate[USf[2, j1][s1, 2]]*(2*CB^2*MW^2*SW^2*
          Conjugate[USf[2, j2][s3, 2]]*IndexDelta[j1, j4]*IndexDelta[j2, j3]*
          USf[2, j1][s4, 2]*USf[2, j2][s2, 2] + Conjugate[USf[2, j2][s3, 1]]*
          IndexDelta[j1, j4]*IndexDelta[j2, j3]*
          (-(CB^2*MW^2*SW^2*USf[2, j1][s4, 2]*USf[2, j2][s2, 1]) + 
           CW^2*MLE[j1]*MLE[j2]*USf[2, j1][s4, 1]*USf[2, j2][s2, 2]) + 
         IndexDelta[j1, j2]*IndexDelta[j3, j4]*
          (2*CB^2*MW^2*SW^2*Conjugate[USf[2, j3][s3, 2]]*USf[2, j1][s2, 2]*
            USf[2, j3][s4, 2] + Conjugate[USf[2, j3][s3, 1]]*
            (-(CB^2*MW^2*SW^2*USf[2, j1][s2, 2]*USf[2, j3][s4, 1]) + 
             CW^2*MLE[j1]*MLE[j3]*USf[2, j1][s2, 1]*USf[2, j3][s4, 2])))))/
     (CB^2*CW^2*MW^2*SW^2)}}, C[S[12, {s1, j1}], -S[12, {s2, j2}], 
   S[11, {j3}], -S[11, {j4}]] == 
  {{((I/4)*EL^2*((IndexDelta[j1, j2]*IndexDelta[j3, j4]*
         ((CW^2 - SW^2)*Conjugate[USf[2, j1][s1, 1]]*USf[2, j1][s2, 1] + 
          2*SW^2*Conjugate[USf[2, j1][s1, 2]]*USf[2, j1][s2, 2]))/CW^2 - 
       (2*IndexDelta[j1, j4]*IndexDelta[j2, j3]*
         (CB^2*MW^2*Conjugate[USf[2, j1][s1, 1]]*USf[2, j2][s2, 1] + 
          Conjugate[USf[2, j1][s1, 2]]*MLE[j1]*MLE[j2]*USf[2, j2][s2, 2]))/
        (CB^2*MW^2)))/SW^2}}, C[S[12, {s1, j1}], -S[12, {s2, j2}], 
   S[13, {s3, j3, o3}], -S[13, {s4, j4, o4}]] == 
  {{((I/12)*EL^2*IndexDelta[j1, j2]*IndexDelta[j3, j4]*IndexDelta[o3, o4]*
      (-2*SW^2*Conjugate[USf[2, j1][s1, 2]]*USf[2, j1][s2, 2]*
        (Conjugate[USf[3, j3][s3, 1]]*USf[3, j3][s4, 1] - 
         4*Conjugate[USf[3, j3][s3, 2]]*USf[3, j3][s4, 2]) + 
       Conjugate[USf[2, j1][s1, 1]]*USf[2, j1][s2, 1]*
        ((1 + 2*CW^2)*Conjugate[USf[3, j3][s3, 1]]*USf[3, j3][s4, 1] - 
         4*SW^2*Conjugate[USf[3, j3][s3, 2]]*USf[3, j3][s4, 2])))/
     (CW^2*SW^2)}}, C[S[11, {j1}], -S[11, {j2}], S[11, {j3}], 
   -S[11, {j4}]] == 
  {{((-I/4)*EL^2*(IndexDelta[j1, j4]*IndexDelta[j2, j3] + 
       IndexDelta[j1, j2]*IndexDelta[j3, j4]))/(CW^2*SW^2)}}, 
 C[S[11, {j1}], -S[11, {j2}], S[13, {s3, j3, o3}], -S[13, {s4, j4, o4}]] == 
  {{((-I/12)*EL^2*IndexDelta[j1, j2]*IndexDelta[j3, j4]*IndexDelta[o3, o4]*
      ((3*CW^2 - SW^2)*Conjugate[USf[3, j3][s3, 1]]*USf[3, j3][s4, 1] + 
       4*SW^2*Conjugate[USf[3, j3][s3, 2]]*USf[3, j3][s4, 2]))/(CW^2*SW^2)}}, 
 C[S[13, {s1, j1, o1}], -S[13, {s2, j2, o2}], S[13, {s3, j3, o3}], 
   -S[13, {s4, j4, o4}]] == 
  {{IndexDelta[j1, j4]*IndexDelta[j2, j3]*((-I)*GS^2*SUNTSum[o2, o3, o4, o1]*
        (Conjugate[USf[3, j1][s1, 1]]*USf[3, j1][s4, 1] - 
         Conjugate[USf[3, j1][s1, 2]]*USf[3, j1][s4, 2])*
        (Conjugate[USf[3, j2][s3, 1]]*USf[3, j2][s2, 1] - 
         Conjugate[USf[3, j2][s3, 2]]*USf[3, j2][s2, 2]) - 
       ((I/36)*EL^2*IndexDelta[o1, o4]*IndexDelta[o2, o3]*
         (Conjugate[USf[3, j1][s1, 1]]*((1 + 8*CW^2)*MW^2*SB^2*
             Conjugate[USf[3, j2][s3, 1]]*USf[3, j1][s4, 1]*
             USf[3, j2][s2, 1] + 2*Conjugate[USf[3, j2][s3, 2]]*
             (9*CW^2*MQU[j1]*MQU[j2]*USf[3, j1][s4, 2]*USf[3, j2][s2, 1] - 
              2*MW^2*SB^2*SW^2*USf[3, j1][s4, 1]*USf[3, j2][s2, 2])) + 
          2*Conjugate[USf[3, j1][s1, 2]]*(8*MW^2*SB^2*SW^2*
             Conjugate[USf[3, j2][s3, 2]]*USf[3, j1][s4, 2]*
             USf[3, j2][s2, 2] + Conjugate[USf[3, j2][s3, 1]]*
             (-2*MW^2*SB^2*SW^2*USf[3, j1][s4, 2]*USf[3, j2][s2, 1] + 
              9*CW^2*MQU[j1]*MQU[j2]*USf[3, j1][s4, 1]*USf[3, j2][s2, 2]))))/
        (CW^2*MW^2*SB^2*SW^2)) + IndexDelta[j1, j2]*IndexDelta[j3, j4]*
      ((-I)*GS^2*SUNTSum[o2, o1, o4, o3]*(Conjugate[USf[3, j1][s1, 1]]*
          USf[3, j1][s2, 1] - Conjugate[USf[3, j1][s1, 2]]*USf[3, j1][s2, 2])*
        (Conjugate[USf[3, j3][s3, 1]]*USf[3, j3][s4, 1] - 
         Conjugate[USf[3, j3][s3, 2]]*USf[3, j3][s4, 2]) - 
       ((I/36)*EL^2*IndexDelta[o1, o2]*IndexDelta[o3, o4]*
         (Conjugate[USf[3, j1][s1, 1]]*((1 + 8*CW^2)*MW^2*SB^2*
             Conjugate[USf[3, j3][s3, 1]]*USf[3, j1][s2, 1]*
             USf[3, j3][s4, 1] + 2*Conjugate[USf[3, j3][s3, 2]]*
             (9*CW^2*MQU[j1]*MQU[j3]*USf[3, j1][s2, 2]*USf[3, j3][s4, 1] - 
              2*MW^2*SB^2*SW^2*USf[3, j1][s2, 1]*USf[3, j3][s4, 2])) + 
          2*Conjugate[USf[3, j1][s1, 2]]*(8*MW^2*SB^2*SW^2*
             Conjugate[USf[3, j3][s3, 2]]*USf[3, j1][s2, 2]*
             USf[3, j3][s4, 2] + Conjugate[USf[3, j3][s3, 1]]*
             (-2*MW^2*SB^2*SW^2*USf[3, j1][s2, 2]*USf[3, j3][s4, 1] + 
              9*CW^2*MQU[j1]*MQU[j3]*USf[3, j1][s2, 1]*USf[3, j3][s4, 2]))))/
        (CW^2*MW^2*SB^2*SW^2))}}}


(* The following definitions of renormalization constants are
   for the on-shell renormalization of the MSSM in a scheme
   similar to A. Denner, Fortschr. d. Physik, 41 (1993) 4.

   The renormalization constants are not directly used by
   FeynArts, and hence do not restrict the generation of diagrams
   and amplitudes in any way. *)

Clear[RenConst]

RenConst[ dMf1[type_, j1_] ] := MassRC[F[type, {j1}]]

RenConst[ dZfL1[type_, j1_, j2_] ] :=
  FieldRC[F[type, {j1}], F[type, {j2}]][[1]]

RenConst[ dZfR1[type_, j1_, j2_] ] :=
  FieldRC[F[type, {j1}], F[type, {j2}]][[2]]

RenConst[ dMZsq1 ] := MassRC[V[2]]

RenConst[ dMWsq1 ] := MassRC[V[3]]

RenConst[ dMHsq1 ] := MassRC[S[1]]

RenConst[ dZAA1 ] := FieldRC[V[1]]

RenConst[ dZAZ1 ] := FieldRC[V[1], V[2]]

RenConst[ dZZA1 ] := FieldRC[V[2], V[1]]

RenConst[ dZZZ1 ] := FieldRC[V[2]]

RenConst[ dZG01 ] := FieldRC[S[2]]

RenConst[ dZW1 ] := FieldRC[V[3]]

RenConst[ dZGp1 ] := FieldRC[S[3]]

RenConst[ dZH1 ] := FieldRC[S[1]]

RenConst[ dTH1 ] := TadpoleRC[S[1]]

RenConst[ dSW1 ] := -$HKSign CW^2/SW/2 (dMZsq1/MZ^2 - dMWsq1/MW^2)

RenConst[ dZe1 ] := -1/2 (dZAA1 - $HKSign SW/CW dZZA1)

