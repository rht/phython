(*
	SM.mod
		Classes model file for the Standard Model
		by Hagen Eck and Sepp Kueblbeck 1995
		last modified 4 May 10 by Thomas Hahn

This file contains the definition of a Classes model for FeynArts.
It needs the Generic model file Lorentz.gen.

When you change things, remember:

-- All particles are arranged in classes. For single particle
   model definitions each particle lives in its own class.

-- For each class the common SelfConjugate behaviour and the
   IndexRange MUST be present in the definitions.

-- IMPORTANT: The coupling matrices MUST be declared in the
   SAME order as the Generic coupling.

Reference:
	Ansgar Denner, "Techniques for the calculation of electroweak
	radiative corrections at one-loop level and results for
	W-physics at LEP200", Fortschr. d. Physik, 41 (1993) 4

Oct 95: one-loop counter terms added by Stefan Bauberger:
	Some corrections and addition of all one-loop counter terms
	according to A. Denner.  The gauge-fixing terms are assumed not
	to be renormalized.  The Denner conventions are extended to
	include field renormalization of the Goldstone bosons.
	The counter terms associated with quark mixing are not well
	tested yet.

Apr 99: Christian Schappacher added colour indices for the quarks

Apr 99:	Terms for ghost sector updated by Ayres Freitas.
	The gauge-fixing terms are still assumed not to be renormalized
	but the renormalized gauge parameters follow the R_xi-gauge.
	In addition, renormalization for the ghost fields is included.
	The 2-loop counter terms for vector-boson selfenergies and for
	the W-nu-l vertex have been added.
	Old versions of the changes of sbau are removed!

Apr 01:	Thomas Hahn added the definitions of the renormalization
	constants a la A. Denner.

This file introduces the following symbols:

	coupling constants and masses:
	------------------------------
	EL:		electron charge (Thomson limit)
	CW, SW:		cosine and sine of Weinberg angle

	MW, MZ, MH:	W, Z, Higgs masses

	MLE:		lepton class mass
	ME, MM, ML:	lepton masses (e, mu, tau)

	MQU:		u-type quark class mass
	MU, MC, MT:	u-type quark masses (up, charm, top)

	MQD:		d-type quark class mass
	MD, MS, MB:	d-type quark masses (down, strange, bottom)

	CKM:		quark mixing matrix
			(set CKM = IndexDelta for no quark-mixing)

	GaugeXi[A, W, Z]: photon, W, Z gauge parameters


	one-loop renormalization constants (RCs):
	-----------------------------------------
	dZe1:		electromagnetic charge RC
	dSW1, dCW1:	Weinberg angle sine/cosine RC

	dZH1, dMHsq1:	Higgs field and mass RC
	dZW1, dMWsq1:	W field and mass RC
	dMZsq1:		Z mass RC
	dZZZ1, dZZA1,
	dZAZ1, dZAA1:	Z and photon field RCs

	dMf1:		fermion mass RCs
	dZfL1, dZfR1:	fermion field RCs

	dCKM1:		quark mixing matrix RCs

	dZG01, dZGp1:	field RC for unphysical scalars
	dUZZ1, dUZA1,
	dUAZ1, dUAA1:	field RCs for photon and Z ghosts
	dUW1:		field RC for +/- ghosts

	two-loop renormalization constants:
	-----------------------------------
	dZe2:           electromagnetic charge RC
	dSW2:           weak mixing angle sine/cosine RC

	dZW2, dMWsq2:	W field and mass RC
	dMZsq2:		Z mass RC
	dZZZ2, dZZA2,
	dZAZ2, dZAA2:	Z and photon field RCs

	dZfL2:		fermion field RCs
*)


IndexRange[ Index[Generation] ] = Range[3]

IndexRange[ Index[Colour] ] = NoUnfold[Range[3]]

IndexStyle[ Index[Generation, i_Integer] ] := Alph[i + 8]

MaxGenerationIndex = 3


ViolatesQ[ q__ ] := Plus[q] =!= 0


mdZfLR1[ type_, j1_, j2_ ] :=
  Mass[F[type, {j1}]]/2 dZfL1[type, j1, j2] +
    Mass[F[type, {j2}]]/2 Conjugate[dZfR1[type, j2, j1]]

mdZfRL1[ type_, j1_, j2_ ] :=
  Mass[F[type, {j1}]]/2 dZfR1[type, j1, j2] +
    Mass[F[type, {j2}]]/2 Conjugate[dZfL1[type, j2, j1]]


(* the leptonic field RCs are diagonal: *)

dZfL1[ type:1 | 2, j1_, j2_ ] :=
  IndexDelta[j1, j2] dZfL1[type, j1, j1] /; j1 =!= j2

dZfR1[ type:1 | 2, j1_, j2_ ] :=
  IndexDelta[j1, j2] dZfR1[type, j1, j1] /; j1 =!= j2


(* some short-hands for fermionic couplings: *)

FermionCharge[1] = 0;
FermionCharge[2] = -1;
FermionCharge[3] = 2/3;
FermionCharge[4] = -1/3

gR[ type_ ] :=
  -SW/CW FermionCharge[type];
gL[ type_ ] :=
  (If[ OddQ[type], 1/2, -1/2 ] - SW^2 FermionCharge[type])/(SW CW);
dgR[ type_ ] :=
  gR[type] (dZe1 + 1/(CW^2 SW) dSW1);
dgL[ type_ ] :=
  If[ OddQ[type], 1/2, -1/2 ]/(SW CW) *
    (dZe1 + (SW^2 - CW^2)/(CW^2 SW) dSW1) + dgR[type]


M$ClassesDescription = {

	(* Leptons (neutrino): I_3 = +1/2, Q = 0 *)
  F[1] == {
	SelfConjugate -> False,
	Indices -> {Index[Generation]},
	Mass -> 0,
	QuantumNumbers -> LeptonNumber,
	PropagatorLabel -> ComposedChar["\\nu", Index[Generation]],
	PropagatorType -> Straight,
	PropagatorArrow -> Forward },

	(* Leptons (electron): I_3 = -1/2, Q = -1 *)
  F[2] == {
	SelfConjugate -> False,
	Indices -> {Index[Generation]},
	Mass -> MLE,
	QuantumNumbers -> {-Charge, LeptonNumber},
	PropagatorLabel -> ComposedChar["e", Index[Generation]],
	PropagatorType -> Straight,
	PropagatorArrow -> Forward },

	(* Quarks (u): I_3 = +1/2, Q = +2/3 *)
  F[3] == {
	SelfConjugate -> False,
	Indices -> {Index[Generation], Index[Colour]},
	Mass -> MQU,
	QuantumNumbers -> 2/3 Charge,
	PropagatorLabel -> ComposedChar["u", Index[Generation]],
	PropagatorType -> Straight,
	PropagatorArrow -> Forward },

	(* Quarks (d): I_3 = -1/2, Q = -1/3 *) 
  F[4] == {
	SelfConjugate -> False,
	Indices -> {Index[Generation], Index[Colour]},
	Mass -> MQD,
	QuantumNumbers -> -1/3 Charge,
	PropagatorLabel -> ComposedChar["d", Index[Generation]],
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
	QuantumNumbers -> -Charge,
	PropagatorLabel -> "W",
	PropagatorType -> Sine,
	PropagatorArrow -> Forward },

(*
  V[4] == {
	SelfConjugate -> True,
	Indices -> {},
	Mass -> MAZ,
	MixingPartners -> {V[1], V[2]},
	PropagatorLabel -> {"\\gamma", "Z"},
	PropagatorType -> Sine,
	PropagatorArrow -> None },
*)

	(* mixing Higgs gauge bosons: Q = 0 *) 
  SV[2] == {
	SelfConjugate -> True,
	Indices -> {},
	Mass -> MZ,
	MixingPartners -> {S[2], V[2]},
	PropagatorLabel -> {ComposedChar["G", Null, "0"], "Z"},
	PropagatorType -> {ScalarDash, Sine},
	PropagatorArrow -> None },

	(* mixing Higgs gauge bosons: charged *) 
  SV[3] == {
	SelfConjugate -> False,
	Indices -> {},
	Mass -> MW,
	QuantumNumbers -> -Charge,
	MixingPartners -> {S[3], V[3]},
	PropagatorLabel -> {"G", "W"},
	PropagatorType -> {ScalarDash, Sine},
	PropagatorArrow -> Forward },

	(* physical Higgs: Q = 0 *) 
  S[1] == {
	SelfConjugate -> True,
	Indices -> {},
	Mass -> MH,
	PropagatorLabel -> "H",
	PropagatorType -> ScalarDash,
	PropagatorArrow -> None },

	(* unphysical Higgs: neutral *) 
  S[2] == {
	SelfConjugate -> True,
	Indices -> {},
	Mass -> MZ,
	PropagatorLabel -> ComposedChar["G", Null, "0"],
	PropagatorType -> ScalarDash,
	PropagatorArrow -> None },

	(* unphysical Higgs: Q = -1 *)  
  S[3] == {
	SelfConjugate -> False,
	Indices -> {},
	Mass -> MW,
	QuantumNumbers -> -Charge,
	PropagatorLabel -> "G",
	PropagatorType -> ScalarDash,
	PropagatorArrow -> Forward },

	(* Ghosts: neutral *) 
  U[1] == {
	SelfConjugate -> False,
	Indices -> {},
	Mass -> 0,
	QuantumNumbers -> GhostNumber,
	PropagatorLabel -> ComposedChar["u", "\\gamma"],
	PropagatorType -> GhostDash,
	PropagatorArrow -> Forward },

  U[2] == {
	SelfConjugate -> False,
	Indices -> {},
	Mass -> MZ,
	QuantumNumbers -> GhostNumber,
	PropagatorLabel -> ComposedChar["u", "Z"],
	PropagatorType -> GhostDash,
	PropagatorArrow -> Forward },

	(* Ghosts: charged *) 
  U[3] == {
	SelfConjugate -> False,
	Indices -> {},
	Mass -> MW,
	QuantumNumbers -> {-Charge, GhostNumber},
	PropagatorLabel -> ComposedChar["u", "-"],
	PropagatorType -> GhostDash,
	PropagatorArrow -> Forward },

  U[4] == {
	SelfConjugate -> False,
	Indices -> {},
	Mass -> MW,
	QuantumNumbers -> {Charge, GhostNumber},
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
MQU[gen_, _] := MQU[gen];
MQD[gen_, _] := MQD[gen]

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

GaugeXi[ V[1] ] = GaugeXi[A];
GaugeXi[ V[2] ] = GaugeXi[Z];
GaugeXi[ V[3] ] = GaugeXi[W];
GaugeXi[ S[1] ] = 1;
GaugeXi[ S[2] ] = GaugeXi[Z];
GaugeXi[ S[3] ] = GaugeXi[W];
GaugeXi[ U[1] ] = GaugeXi[A];
GaugeXi[ U[2] ] = GaugeXi[Z];
GaugeXi[ U[3] ] = GaugeXi[W];
GaugeXi[ U[4] ] = GaugeXi[W]


M$CouplingMatrices = {

	(* V-V:  G(+) . { -g[mu, nu] mom^2, g[mu, nu], -mom[mu] mom[nu] } *)

  C[ -V[3], V[3] ] == I *
    { {0, dZW1, dZW2},
      {0, MW^2 dZW1 + dMWsq1, MW^2 dZW2 + dMWsq2 + dMWsq1 dZW1},
      {0, -dZW1, -dZW2} },

  C[ V[2], V[2] ] == I *
    { {0, dZZZ1, dZZZ2 + 1/4 dZAZ1^2},
      {0, MZ^2 dZZZ1 + dMZsq1, MZ^2 dZZZ2 + dMZsq2 + dMZsq1 dZZZ1},
      {0, -dZZZ1, -dZZZ2 - 1/4 dZAZ1^2} },

  C[ V[1], V[1] ] == I *
    { {0, dZAA1, dZAA2 + 1/4 dZZA1^2},
      {0, 0, 1/4 MZ^2 dZZA1^2},
      {0, -dZAA1, -dZAA2 - 1/4 dZZA1^2} },

  C[ V[1], V[2] ] == I *
    { {0, dZAZ1/2 + dZZA1/2,
    	  (dZAZ2 + dZZA2 + 1/2 dZZA1 dZZZ1 + 1/2 dZAZ1 dZAA1)/2},
      {0, MZ^2 dZZA1/2,
          (MZ^2 dZZA2 + 1/2 MZ^2 dZZZ1 dZZA1 + dMZsq1 dZZA1)/2},
      {0, -dZAZ1/2 - dZZA1/2,
          -(dZAZ2 + dZZA2 + 1/2 dZZA1 dZZZ1 + 1/2 dZAZ1 dZAA1)/2} },

	(* S-V:  G(+) . { mom1[mu], mom2[mu] } *)

  C[ S[3], -V[3] ] == I MW/4 *
    { {0, -dZW1 - dZGp1 - dMWsq1/MW^2},
      {0, dZW1 + dZGp1 + dMWsq1/MW^2} },

  C[ -S[3], V[3] ] == I MW/4 *
    { {0, dZW1 + dZGp1 + dMWsq1/MW^2},
      {0, -dZW1 - dZGp1 - dMWsq1/MW^2} },

  C[ S[2], V[2] ] == MZ/4 *
    { {0, dZZZ1 + dZG01 + dMZsq1/MZ^2},
      {0, -dZZZ1 - dZG01 - dMZsq1/MZ^2} },

  C[ S[2], V[1] ] == MZ/4 *
    { {0, dZZA1},
      {0, -dZZA1} },

	(* S-S:  G(+) . { -mom^2, 1 } *)

  C[ S[1], S[1] ] == -I *
    { {0, dZH1},
      {0, dMHsq1 + MH^2 dZH1} },

  C[ S[2], S[2] ] == -I *
    { {0, dZG01},
      {0, -EL/(2 MW SW) dTH1} },

  C[ S[3], -S[3] ] == -I *
    { {0, dZGp1},
      {0, -EL/(2 MW SW) dTH1} },

	(* U-U:  G(+) . { -mom^2, 1 } *)

  C[ U[1], -U[1] ] == -I *
    { {0, -dZAA1/2 + dUAA1},
      {0, 0} },
  
  C[ U[2], -U[2] ] == -I *
    { {0, -dZZZ1/2 + dUZZ1},
      {0, GaugeXi[Z] (MZ^2 (-dZG01/2 + dUZZ1) + dMZsq1/2) } },
  
  C[ U[2], -U[1] ] == -I *
    { {0, -dZAZ1/2 + dUAZ1},
      {0, 0} },
  
  C[ U[1], -U[2] ] == -I *
    { {0, -dZZA1/2 + dUZA1},
      {0, GaugeXi[Z] MZ^2 dUZA1} },
  
  C[ U[3], -U[3] ] == -I *
    { {0, -dZW1/2 + dUW1},
      {0, GaugeXi[W] (MW^2 (-dZGp1/2 + dUW1) + dMWsq1/2) } },

  C[ U[4], -U[4] ] == -I *
    { {0, -dZW1/2 + dUW1},
      {0, GaugeXi[W] (MW^2 (-dZGp1/2 + dUW1) + dMWsq1/2) } },

	(* F-F:  G(+) . { slash[mom1] omega[-], slash[mom2] omega[+],
	                  omega[-], omega[+] } *)

  C[ -F[1, {j1}], F[1, {j2}] ] == I *
    { {0, -AddHC[dZfL1[1, j1, j2]]},
      {0, AddHC[dZfR1[1, j1, j2]]},
      {0, 0},
      {0, 0} },

  C[ -F[2, {j1}], F[2, {j2}] ] == I *
    { {0, -AddHC[dZfL1[2, j1, j2]]},
      {0, AddHC[dZfR1[2, j1, j2]]},
      {0, -mdZfLR1[2, j1, j2] - IndexDelta[j1, j2] dMf1[2, j1]},
      {0, -mdZfRL1[2, j1, j2] - IndexDelta[j1, j2] dMf1[2, j1]} },

  C[ -F[3, {j1, o1}], F[3, {j2, o2}] ] == I IndexDelta[o1, o2] *
    { {0, -AddHC[dZfL1[3, j1, j2]]},
      {0, AddHC[dZfR1[3, j1, j2]]},
      {0, -mdZfLR1[3, j1, j2] - IndexDelta[j1, j2] dMf1[3, j1]},
      {0, -mdZfRL1[3, j1, j2] - IndexDelta[j1, j2] dMf1[3, j1]} },

  C[ -F[4, {j1, o1}], F[4, {j2, o2}] ] == I IndexDelta[o1, o2] *
    { {0, -AddHC[dZfL1[4, j1, j2]]},
      {0, AddHC[dZfR1[4, j1, j2]]},
      {0, -mdZfLR1[4, j1, j2] - IndexDelta[j1, j2] dMf1[4, j1]},
      {0, -mdZfRL1[4, j1, j2] - IndexDelta[j1, j2] dMf1[4, j1]} },

	(* V-V-V-V:  G(+) . { g[mu1, mu2] g[mu3, mu4],
	                      g[mu1, mu4] g[mu2, mu3],
	                      g[mu1, mu3] g[mu2, mu4] } *)

  C[ -V[3], -V[3], V[3], V[3] ] == I EL^2/SW^2 *
    { {2, 4 dZe1 - 4 dSW1/SW + 4 dZW1}, 
      {-1, -2 dZe1 + 2 dSW1/SW - 2 dZW1},
      {-1, -2 dZe1 + 2 dSW1/SW - 2*dZW1} },

  C[ -V[3], V[3], V[2], V[2] ] == -I EL^2 CW^2/SW^2 *
    { {2, 4 dZe1 - 4 dSW1/(SW CW^2) + 2 dZW1 + 2 dZZZ1 - 2 dZAZ1 SW/CW}, 
      {-1, -2 dZe1 + 2 dSW1/(SW CW^2) - dZW1 - dZZZ1 + dZAZ1 SW/CW},
      {-1, -2 dZe1 + 2 dSW1/(SW CW^2) - dZW1 - dZZZ1 + dZAZ1 SW/CW} },

  C[ -V[3], V[3], V[1], V[2] ] == I EL^2 CW/SW *
    { {2, 4 dZe1 - 2 dSW1/(SW CW^2) + 2 dZW1 +
            dZZZ1 + dZAA1 - SW/CW dZAZ1 - CW/SW dZZA1},
      {-1, -2 dZe1 + dSW1/(SW CW^2) - dZW1 -
            dZZZ1/2 - dZAA1/2 + SW/CW dZAZ1/2 + CW/SW dZZA1/2},
      {-1, -2 dZe1 + dSW1/(SW CW^2) - dZW1 -
            dZZZ1/2 - dZAA1/2 + SW/CW dZAZ1/2 + CW/SW dZZA1/2} },

  C[ -V[3], V[3], V[1], V[1] ] == -I EL^2 *
    { {2, 4 dZe1 + 2 dZW1 + 2 dZAA1 - 2 CW/SW dZZA1}, 
      {-1, -2 dZe1 - dZW1 - dZAA1 + CW/SW dZZA1},
      {-1, -2 dZe1 - dZW1 - dZAA1 + CW/SW dZZA1} },

	(* V-V-V:  G(-) . (g[mu1, mu2] (p2 - p1)_mu3 +
	                   g[mu2, mu3] (p3 - p2)_mu1 +
	                   g[mu3, mu1] (p1 - p3)_mu2) *)

  C[ V[1], -V[3], V[3] ] == -I EL *
    { {1, dZe1 + dZW1 + dZAA1/2 - CW/SW dZZA1/2} },

  C[ V[2], -V[3], V[3] ] == I EL CW/SW *
    { {1, dZe1 - dSW1/(SW CW^2) + dZW1 + dZZZ1/2 - SW/CW dZAZ1/2} },

	(* S-S-S-S:  G(+) . 1 *)

  C[ S[1], S[1], S[1], S[1] ] == -3 I EL^2 MH^2/(4 SW^2 MW^2) *
    { {1, 2 dZe1 - 2 dSW1/SW + dMHsq1/MH^2 + EL/(2 SW MW MH^2) dTH1 -
            dMWsq1/MW^2 + 2 dZH1} },

  C[ S[1], S[1], S[2], S[2] ] == -I EL^2 MH^2/(4 SW^2 MW^2) *
    { {1, 2 dZe1 - 2 dSW1/SW + dMHsq1/MH^2 + EL/(2 SW MW MH^2) dTH1 -
            dMWsq1/MW^2 + dZH1 + dZG01} },

  C[ S[1], S[1], S[3], -S[3] ] == -I EL^2 MH^2/(4 SW^2 MW^2) *
    { {1, 2 dZe1 - 2 dSW1/SW + dMHsq1/MH^2 + EL/(2 SW MW MH^2) dTH1 -
            dMWsq1/MW^2 + dZH1 + dZGp1} },

  C[ S[2], S[2], S[2], S[2] ] == -3 I EL^2 MH^2/(4 SW^2 MW^2) *
    { {1, 2 dZe1 - 2 dSW1/SW + dMHsq1/MH^2 + EL/(2 SW MW MH^2) dTH1 -
            dMWsq1/MW^2 + 2 dZG01} },

  C[ S[2], S[2], S[3], -S[3] ] == -I EL^2 MH^2/(4 SW^2 MW^2) *
    { {1, 2 dZe1 - 2 dSW1/SW + dMHsq1/MH^2 + EL/(2 SW MW MH^2) dTH1 -
            dMWsq1/MW^2 + dZG01 + dZGp1} },

  C[ S[3], S[3], -S[3], -S[3] ] == -I EL^2 MH^2/(2 SW^2 MW^2) *
    { {1, 2 dZe1 - 2 dSW1/SW + dMHsq1/MH^2 + EL/(2 SW MW MH^2) dTH1 -
            dMWsq1/MW^2 + 2 dZGp1} },

	(* S-S-S:  G(+) . 1 *)

  C[ S[1], S[1], S[1] ] == -3 I EL MH^2/(2 SW MW) *
    { {1, dZe1 - dSW1/SW + dMHsq1/MH^2 + EL/(2 SW MW MH^2) dTH1 -
            dMWsq1/(2 MW^2) + 3/2 dZH1} },
 
  C[ S[1], S[2], S[2] ] == -I EL MH^2/(2 SW MW) *
    { {1, dZe1 - dSW1/SW + dMHsq1/MH^2 + EL/(2 SW MW MH^2) dTH1 -
            dMWsq1/(2 MW^2) + dZH1/2 + dZG01} },

  C[ S[3], S[1], -S[3] ] == -I EL MH^2/(2 SW MW) *
    { {1, dZe1 - dSW1/SW + dMHsq1/MH^2 + EL/(2 SW MW MH^2) dTH1 -
            dMWsq1/(2 MW^2) + dZH1/2 + dZGp1} },

	(* S-S-V-V:  G(+) . g[mu3, mu4] *)

  C[ S[1], S[1], V[3], -V[3] ] == I EL^2/(2 SW^2) *
    { {1, 2 dZe1 - 2 dSW1/SW + dZW1 + dZH1} },

  C[ S[2], S[2], V[3], -V[3] ] == I EL^2/(2 SW^2) *
    { {1, 2 dZe1 - 2 dSW1/SW + dZW1 + dZG01} },

  C[ S[3], -S[3], V[3], -V[3] ] == I EL^2/(2 SW^2) *
    { {1, 2 dZe1 - 2 dSW1/SW + dZW1 + dZGp1} },

  C[ S[3], -S[3], V[2], V[2] ] == I EL^2 (SW^2 - CW^2)^2/(2 CW^2 SW^2) *
    { {1, 2 dZe1 + 2/(SW CW^2 (SW^2 - CW^2)) dSW1 + dZZZ1 + dZGp1 +
            2 SW CW/(SW^2 - CW^2) dZAZ1} },

  C[ S[3], -S[3], V[1], V[2] ] == I EL^2 (SW^2 - CW^2)/(CW SW) *
    { {1, 2 dZe1 + dSW1/(SW CW^2 (SW^2 - CW^2)) + dZZZ1/2 + dZAA1/2 +
            dZGp1 + (SW^2 - CW^2)/(4 SW CW) dZZA1 +
            SW CW/(SW^2 - CW^2) dZAZ1} },

  C[ S[3], -S[3], V[1], V[1] ] == 2 I EL^2 *
    { {1, 2 dZe1 + dZAA1 + dZGp1 + (SW^2 - CW^2)/(2 SW CW) dZZA1} },

  C[ S[1], S[1], V[2], V[2] ] == I EL^2/(2 CW^2 SW^2) *
    { {1, 2 dZe1 + 2 (SW^2 - CW^2)/(SW CW^2) dSW1 + dZZZ1 + dZH1} },

  C[ S[2], S[2], V[2],  V[2] ] == I EL^2/(2 CW^2 SW^2) *
    { {1, 2 dZe1 + 2 (SW^2 - CW^2)/(SW CW^2) dSW1 + dZZZ1 + dZG01} },

  C[ S[1], S[1], V[1], V[2] ] == I EL^2/(2 CW^2 SW^2) *
    { {0, dZZA1/2} },

  C[ S[2], S[2], V[1], V[2] ] == I EL^2/(2 CW^2 SW^2) *
    { {0, dZZA1/2} },

  C[ S[1], -S[3], V[3], V[2] ] == -I EL^2/(2 CW) *
    { {1, 2 dZe1 - dCW1/CW + dZW1/2 + dZH1/2 + dZGp1/2 +
            dZZZ1/2 + CW/SW dZAZ1/2} },

  C[ S[1], S[3], -V[3], V[2] ] == -I EL^2/(2 CW) *
    { {1, 2 dZe1 - dCW1/CW + dZW1/2 + dZH1/2 + dZGp1/2 +
            dZZZ1/2 + CW/SW dZAZ1/2} },

  C[ S[1], S[3], -V[3], V[1] ] == -I EL^2/(2 SW) *
    { {1, 2 dZe1 - dSW1/SW + dZW1/2 + dZH1/2 + dZGp1/2 +
            dZAA1/2 + SW/CW dZZA1/2} },

  C[ S[1], -S[3], V[3], V[1] ] == -I EL^2/(2 SW) *
    { {1, 2 dZe1 - dSW1/SW + dZW1/2 + dZH1/2 + dZGp1/2 +
            dZAA1/2 + SW/CW dZZA1/2} },

  C[ S[3], S[2], V[2], -V[3] ] == EL^2/(2 CW) *
    { {1, 2 dZe1 - dCW1/CW + dZW1/2 + dZZZ1/2 + dZGp1/2 + dZG01/2 +
            CW/SW dZAZ1/2} },

  C[ -S[3], S[2], V[2], V[3] ] == -EL^2/(2 CW) *
    { {1, 2 dZe1 - dCW1/CW + dZW1/2 + dZZZ1/2 + dZGp1/2 + dZG01/2 +
            CW/SW dZAZ1/2} },

  C[ S[3], S[2], V[1], -V[3] ] == EL^2/(2 SW) *
    { {1, 2 dZe1 - dSW1/SW + dZW1/2 + dZAA1/2 + dZGp1/2 + dZG01/2 +
            SW/CW dZZA1/2} },

  C[ -S[3], S[2], V[1], V[3] ] == -EL^2/(2 SW) *
    { {1, 2 dZe1 - dSW1/SW + dZW1/2 + dZAA1/2 + dZGp1/2 + dZG01/2 +
            SW/CW dZZA1/2} },

	(* S-S-V:  G(-) . (p1 - p2)_mu3 *)

  C[ S[2], S[1], V[1] ] == EL/(2 CW SW) *
    { {0, dZZA1/2} },

  C[ S[2], S[1], V[2] ] == EL/(2 CW SW) *
    { {1, dZe1 + (SW^2 - CW^2)/(CW^2 SW) dSW1 + dZH1/2 + dZZZ1/2 +
            dZG01/2} },

  C[ -S[3], S[3], V[1] ] == -I EL *
    { {1, dZe1 + dZAA1/2 + dZGp1 + (SW^2 - CW^2)/(2 SW CW) dZZA1/2} },

  C[ -S[3], S[3], V[2] ] == -I EL (SW^2 - CW^2)/(2 CW SW) *
    { {1, dZe1 + dSW1/((SW^2 - CW^2) CW^2 SW) + dZZZ1/2 + dZGp1 +
            2 SW CW/(SW^2 - CW^2) dZAZ1/2} },

  C[ S[3], S[1], -V[3] ] == -I EL/(2 SW) *
    { {1, dZe1 - dSW1/SW + dZW1/2 + dZH1/2 + dZGp1/2} },

  C[ -S[3], S[1], V[3] ] == I EL/(2 SW) *
    { {1, dZe1 - dSW1/SW + dZW1/2 + dZH1/2 + dZGp1/2} },

  C[ S[3], S[2], -V[3] ] == EL/(2 SW) *
    { {1, dZe1 - dSW1/SW + dZW1/2 + dZGp1/2 + dZG01/2} },

  C[ -S[3], S[2], V[3] ] == EL/(2 SW) *
    { {1, dZe1 - dSW1/SW + dZW1/2 + dZGp1/2 + dZG01/2} },

	(* S-V-V:  G(+) . g[mu2, mu3] *)

  C[ S[1], -V[3], V[3] ] == I EL MW/SW *
    { {1, dZe1 - dSW1/SW + dMWsq1/(2 MW^2) + dZH1/2 + dZW1} },

  C[ S[1], V[2], V[2] ] == I EL MW/(SW CW^2) *
    { {1, dZe1 + (2 SW^2 - CW^2)/(CW^2 SW) dSW1 + dMWsq1/(2 MW^2) +
            dZH1/2 + dZZZ1} },

  C[ S[1], V[2], V[1] ] == I EL MW/(SW CW^2) *
    { {0, dZZA1/2} },

  C[ -S[3], V[3], V[2] ] == -I EL MW SW/CW *
    { {1, dZe1 + dSW1/(CW^2 SW) + dMWsq1/(2 MW^2) + dZW1/2 + dZZZ1/2 +
            dZGp1/2 + CW/SW dZAZ1/2} },

  C[ S[3], -V[3], V[2] ] == -I EL MW SW/CW *
    { {1, dZe1 + dSW1/(CW^2 SW) + dMWsq1/(2 MW^2) + dZW1/2 + dZZZ1/2 +
            dZGp1/2 + CW/SW dZAZ1/2} },

  C[ -S[3], V[3], V[1] ] == -I EL MW *
    { {1, dZe1 + dMWsq1/(2 MW^2) + dZW1/2 + dZAA1/2 + dZGp1/2 +
            SW/CW dZZA1/2} },

  C[ S[3], -V[3], V[1] ] == -I EL MW *
    { {1, dZe1 + dMWsq1/(2 MW^2) + dZW1/2 + dZAA1/2 + dZGp1/2 +
            SW/CW dZZA1/2} },

	(* F-F-V:  G(-) . { gamma[mu3] omega[-], gamma[mu3] omega[+] } *)

  C[ -F[1, {j1}], F[1, {j2}], V[1] ] == I EL *
    { {0, gL[1] IndexDelta[j1, j2] dZZA1/2},
      {0, 0} },

  C[ -F[2, {j1}], F[2, {j2}], V[1] ] == I EL *
    { {-FermionCharge[2] IndexDelta[j1, j2],
        -FermionCharge[2] *
          (IndexDelta[j1, j2] (dZe1 + dZAA1/2) + AddHC[dZfL1[2, j1, j2]]) +
          gL[2] IndexDelta[j1, j2] dZZA1/2},
      {-FermionCharge[2] IndexDelta[j1, j2],
        -FermionCharge[2] *
          (IndexDelta[j1, j2] (dZe1 + dZAA1/2) + AddHC[dZfR1[2, j1, j2]]) +
          gR[2] IndexDelta[j1, j2] dZZA1/2} },

  C[ -F[3, {j1, o1}], F[3, {j2, o2}], V[1] ] == I EL IndexDelta[o1, o2] *
    { {-FermionCharge[3] IndexDelta[j1, j2],
        -FermionCharge[3] *
          (IndexDelta[j1, j2] (dZe1 + dZAA1/2) + AddHC[dZfL1[3, j1, j2]]) +
          gL[3] IndexDelta[j1, j2] dZZA1/2},
      {-FermionCharge[3] IndexDelta[j1, j2],
        -FermionCharge[3] *
          (IndexDelta[j1, j2] (dZe1 + dZAA1/2) + AddHC[dZfR1[3, j1, j2]]) +
          gR[3] IndexDelta[j1, j2] dZZA1/2} },

  C[ -F[4, {j1, o1}], F[4, {j2, o2}], V[1] ] == I EL IndexDelta[o1, o2] *
    { {-FermionCharge[4] IndexDelta[j1, j2],
        -FermionCharge[4] *
          (IndexDelta[j1, j2] (dZe1 + dZAA1/2) + AddHC[dZfL1[4, j1, j2]]) +
          gL[4] IndexDelta[j1, j2] dZZA1/2},
      {-FermionCharge[4] IndexDelta[j1, j2],
        -FermionCharge[4] *
          (IndexDelta[j1, j2] (dZe1 + dZAA1/2) + AddHC[dZfR1[4, j1, j2]]) +
          gR[4] IndexDelta[j1, j2] dZZA1/2} },

  C[ -F[1, {j1}], F[1, {j2}], V[2] ] == I EL *
    { {gL[1] IndexDelta[j1, j2],
        IndexDelta[j1, j2] (gL[1] dZZZ1/2 + dgL[1]) +
        gL[1] AddHC[dZfL1[1, j1, j2]]},
      {0, 0} },

  C[ -F[2, {j1}], F[2, {j2}], V[2] ] == I EL *
    { {gL[2] IndexDelta[j1, j2],
        IndexDelta[j1, j2] *
          (gL[2] dZZZ1/2 + dgL[2] - FermionCharge[2] dZAZ1/2) +
          gL[2] AddHC[dZfL1[2, j1, j2]]},
      {gR[2] IndexDelta[j1, j2],
        IndexDelta[j1, j2] *
          (gR[2] dZZZ1/2 + dgR[2] - FermionCharge[2] dZAZ1/2) +
          gR[2] AddHC[dZfR1[2, j1, j2]]} },

  C[ -F[3, {j1, o1}], F[3, {j2, o2}], V[2] ] == I EL IndexDelta[o1, o2] *
    { {gL[3] IndexDelta[j1, j2],
        IndexDelta[j1, j2] *
          (gL[3] dZZZ1/2 + dgL[3] - FermionCharge[3] dZAZ1/2) +
          gL[3] AddHC[dZfL1[3, j1, j2]]},
      {gR[3] IndexDelta[j1, j2],
        IndexDelta[j1, j2] *
          (gR[3] dZZZ1/2 + dgR[3] - FermionCharge[3] dZAZ1/2) +
          gR[3] AddHC[dZfR1[3, j1, j2]]} },

  C[ -F[4, {j1, o1}], F[4, {j2, o2}], V[2] ] == I EL IndexDelta[o1, o2] *
    { {gL[4] IndexDelta[j1, j2],
        IndexDelta[j1, j2] *
          (gL[4] dZZZ1/2 + dgL[4] - FermionCharge[4] dZAZ1/2) +
          gL[4] AddHC[dZfL1[4, j1, j2]]},
      {gR[4] IndexDelta[j1, j2],
        IndexDelta[j1, j2] *
          (gR[4] dZZZ1/2 + dgR[4] - FermionCharge[4] dZAZ1/2) +
          gR[4] AddHC[dZfR1[4, j1, j2]]} },

  C[ -F[1, {j1}], F[2, {j2}], -V[3] ] ==
    I EL/(Sqrt[2] SW) IndexDelta[j1, j2] *
    { {1, dZe1 - dSW1/SW + dZW1/2 +
            Conjugate[dZfL1[1, j1, j1]]/2 + dZfL1[2, j1, j1]/2,
          dZe2 - dSW2/SW +
            1/2 (dZW2 + Conjugate[dZfL2[1, j1, j1]] + dZfL2[2, j1, j1]) +
            (dSW1/SW)^2 - dSW1/SW dZe1 -
            1/8 (dZW1^2 + Conjugate[dZfL1[1, j1, j1]]^2 +
              dZfL1[2, j1, j1]^2) +
            (dZe1 - dSW1/SW) *
              1/2 (dZW1 + Conjugate[dZfL1[1, j1, j1]] + dZfL1[2, j1, j1]) +
            1/4 (dZW1 dZfL1[2, j1, j1] + dZW1 Conjugate[dZfL1[1, j1, j1]] +
                   Conjugate[dZfL1[1, j1, j1]] dZfL1[2, j1, j1]) },
      {0, 0, 0} },

  C[ -F[2, {j1}], F[1, {j2}], V[3] ] ==
    I EL/(Sqrt[2] SW) IndexDelta[j1, j2] *
    { {1, dZe1 - dSW1/SW + dZW1/2 +
            dZfL1[1, j1, j1]/2 + Conjugate[dZfL1[2, j1, j1]]/2,
          dZe2 - dSW2/SW +
            1/2 (dZW2 + dZfL2[1, j1, j1] + Conjugate[dZfL2[2, j1, j1]]) +
            (dSW1/SW)^2 - dSW1/SW dZe1 -
            1/8 (dZW1^2 + dZfL1[1, j1, j1]^2 +
              Conjugate[dZfL1[2, j1, j1]]^2) +
            (dZe1 - dSW1/SW) *
              1/2 (dZW1 + dZfL1[1, j1, j1] + Conjugate[dZfL1[2, j1, j1]]) +
	    1/4 (dZW1 Conjugate[dZfL1[2, j1, j1]] + dZW1 dZfL1[1, j1, j1] +
              dZfL1[1, j1, j1] Conjugate[dZfL1[2, j1, j1]]) },
      {0, 0, 0} },

  C[ -F[3, {j1, o1}], F[4, {j2, o2}], -V[3] ] ==
    I EL/(Sqrt[2] SW) IndexDelta[o1, o2] *
    { {CKM[j1, j2],
        CKM[j1, j2] (dZe1 - dSW1/SW + dZW1/2) + dCKM1[j1, j2] +
        1/2 IndexSum[
          Conjugate[dZfL1[3, gn, j1]] CKM[gn, j2] +
          CKM[j1, gn] dZfL1[4, gn, j2],
        {gn, MaxGenerationIndex}]},
      {0, 0} },

  C[ -F[4, {j2, o2}], F[3, {j1, o1}], V[3] ] ==
    I EL/(Sqrt[2] SW) IndexDelta[o1, o2] *
    { {Conjugate[CKM[j1, j2]],
        Conjugate[CKM[j1, j2]] (dZe1 - dSW1/SW + dZW1/2) +
          Conjugate[dCKM1[j1, j2]] +
          1/2 IndexSum[
            Conjugate[dZfL1[4, gn, j2]] Conjugate[CKM[j1, gn]] +
            Conjugate[CKM[gn, j2]] dZfL1[3, gn, j1],
          {gn, MaxGenerationIndex}]},
      {0, 0} },

	(* F-F-S:  G(+) . { omega[-], omega[+] } *)

  C[ -F[2, {j1}], F[2, {j2}], S[1] ] == -I EL/(2 SW MW) *
    { {Mass[F[2, {j1}]] IndexDelta[j1, j2],
        Mass[F[2, {j1}]] IndexDelta[j1, j2] (dZe1 - dSW1/SW +
          dMf1[2, j1]/Mass[F[2, {j1}]] - dMWsq1/(2 MW^2) + dZH1/2) +
          mdZfLR1[2, j1, j2]},
      {Mass[F[2, {j1}]] IndexDelta[j1, j2],
        Mass[F[2, {j1}]] IndexDelta[j1, j2] (dZe1 - dSW1/SW +
          dMf1[2, j1]/Mass[F[2, {j1}]] - dMWsq1/(2 MW^2) + dZH1/2) +
          mdZfRL1[2, j1, j2]} },

  C[ -F[3, {j1, o1}], F[3, {j2, o2}], S[1] ] ==
    -I EL/(2 SW MW) IndexDelta[o1, o2] *
    { {Mass[F[3, {j1}]] IndexDelta[j1, j2],
        Mass[F[3, {j1}]] IndexDelta[j1, j2] (dZe1 - dSW1/SW +
          dMf1[3, j1]/Mass[F[3, {j1}]] - dMWsq1/(2 MW^2) + dZH1/2) +
          mdZfLR1[3, j1, j2]},
      {Mass[F[3, {j1}]] IndexDelta[j1, j2],
        Mass[F[3, {j1}]] IndexDelta[j1, j2] (dZe1 - dSW1/SW +
          dMf1[3, j1]/Mass[F[3, {j1}]] - dMWsq1/(2 MW^2) + dZH1/2) +
          mdZfRL1[3, j1, j2]} },

  C[ -F[4, {j1, o1}], F[4, {j2, o2}], S[1] ] ==
    -I EL/(2 SW MW) IndexDelta[o1, o2] *
    { {Mass[F[4, {j1}]] IndexDelta[j1, j2],
        Mass[F[4, {j1}]] IndexDelta[j1, j2] (dZe1 - dSW1/SW +
          dMf1[4, j1]/Mass[F[4, {j1}]] - dMWsq1/(2 MW^2) + dZH1/2) +
          mdZfLR1[4, j1, j2]},
      {Mass[F[4, {j1}]] IndexDelta[j1, j2],
        Mass[F[4, {j1}]] IndexDelta[j1, j2] (dZe1 - dSW1/SW +
          dMf1[4, j1]/Mass[F[4, {j1}]] - dMWsq1/(2 MW^2) + dZH1/2) +
          mdZfRL1[4, j1, j2]} },

  C[ -F[2, {j1}], F[2, {j2}], S[2] ] == -EL/(2 SW MW) *
    { {Mass[F[2, {j1}]] IndexDelta[j1, j2],
        Mass[F[2, {j1}]] IndexDelta[j1, j2] (dZe1 - dSW1/SW +
          dMf1[2, j1]/Mass[F[2, {j1}]] - dMWsq1/(2 MW^2) + dZG01/2) +
          mdZfLR1[2, j1, j2]},
      {-Mass[F[2, {j1}]] IndexDelta[j1, j2],
        -Mass[F[2, {j1}]] IndexDelta[j1, j2] (dZe1 - dSW1/SW +
          dMf1[2, j1]/Mass[F[2, {j1}]] - dMWsq1/(2 MW^2) + dZG01/2) -
          mdZfRL1[2, j1, j2]} },

  C[ -F[3, {j1, o1}], F[3, {j2, o2}], S[2] ] ==
    EL/(2 SW MW) IndexDelta[o1, o2] *
    { {Mass[F[3, {j1}]] IndexDelta[j1, j2],
        Mass[F[3, {j1}]] IndexDelta[j1, j2] (dZe1 - dSW1/SW +
          dMf1[3, j1]/Mass[F[3, {j1}]] - dMWsq1/(2 MW^2) + dZG01/2) +
          mdZfLR1[3, j1, j2]},
      {-Mass[F[3, {j1}]] IndexDelta[j1, j2],
        -Mass[F[3, {j1}]] IndexDelta[j1, j2] (dZe1 - dSW1/SW +
          dMf1[3, j1]/Mass[F[3, {j1}]] - dMWsq1/(2 MW^2) + dZG01/2) -
          mdZfRL1[3, j1, j2]} },

  C[ -F[4, {j1, o1}], F[4, {j2, o2}], S[2] ] ==
    -EL/(2 SW MW) IndexDelta[o1, o2] *
    { {Mass[F[4, {j1}]] IndexDelta[j1, j2],
        Mass[F[4, {j1}]] IndexDelta[j1, j2] (dZe1 - dSW1/SW +
          dMf1[4, j1]/Mass[F[4, {j1}]] - dMWsq1/(2 MW^2) + dZG01/2) +
          mdZfLR1[4, j1, j2]},
      {-Mass[F[4, {j1}]] IndexDelta[j1, j2],
        -Mass[F[4, {j1}]] IndexDelta[j1, j2] (dZe1 - dSW1/SW +
          dMf1[4, j1]/Mass[F[4, {j1}]] - dMWsq1/(2 MW^2) + dZG01/2) -
          mdZfRL1[4, j1, j2]} },

  C[ -F[3, {j1, o1}], F[4, {j2, o2}], -S[3] ] ==
    I EL/(Sqrt[2] SW MW) IndexDelta[o1, o2] *
    { {Mass[F[3, {j1}]] CKM[j1, j2],
        Mass[F[3, {j1}]] *
          (CKM[j1, j2] (dZe1 - dSW1/SW + dMf1[3, j1]/Mass[F[3, {j1}]] -
            dMWsq1/(2 MW^2) + dZGp1/2) + dCKM1[j1, j2]) +
          1/2 IndexSum[
            Mass[F[3, {gn}]] Conjugate[dZfR1[3, gn, j1]] CKM[gn, j2] +
            Mass[F[3, {j1}]] CKM[j1, gn] dZfL1[4, gn, j2],
          {gn, MaxGenerationIndex}]},
      {-Mass[F[4, {j2}]] CKM[j1, j2],
        -Mass[F[4, {j2}]] *
          (CKM[j1, j2] (dZe1 - dSW1/SW + dMf1[4, j2]/Mass[F[4, {j2}]] -
            dMWsq1/(2 MW^2) + dZGp1/2) + dCKM1[j1, j2]) -
          1/2 IndexSum[
            Mass[F[4, {j2}]] Conjugate[dZfL1[3, gn, j1]] CKM[gn, j2] +
            Mass[F[4, {gn}]] CKM[j1, gn] dZfR1[4, gn, j2],
          {gn, MaxGenerationIndex}]} },

  C[ -F[4, {j2, o2}], F[3, {j1, o1}], S[3] ] ==
    -I EL/(Sqrt[2] SW MW) IndexDelta[o1, o2] *
    { {Mass[F[4, {j2}]] Conjugate[CKM[j1, j2]],
        Mass[F[4, {j2}]] (
            Conjugate[CKM[j1, j2]] (dZe1 - dSW1/SW +
              dMf1[4, j2]/Mass[F[4, {j2}]] - dMWsq1/(2 MW^2) + dZGp1/2) +
            Conjugate[dCKM1[j1, j2]] ) +
          1/2 IndexSum[
            Mass[F[4, {gn}]] Conjugate[dZfR1[4, gn, j2]] *
              Conjugate[CKM[j1, gn]] +
            Mass[F[4, {j2}]] Conjugate[CKM[gn, j2]] dZfL1[3, gn, j1],
          {gn, MaxGenerationIndex}]},
      {-Mass[F[3, {j1}]] Conjugate[CKM[j1, j2]],
        -Mass[F[3, {j1}]] (
            Conjugate[CKM[j1, j2]] (dZe1 - dSW1/SW +
              dMf1[3, j2]/Mass[F[3, {j2}]] - dMWsq1/(2 MW^2) + dZGp1/2) +
            Conjugate[dCKM1[j1, j2]] ) -
          1/2 IndexSum[
            Mass[F[3, {j1}]] Conjugate[dZfL1[4, gn, j2]] *
              Conjugate[CKM[j1, gn]] +
            Mass[F[3, {gn}]] Conjugate[CKM[gn, j2]] *
              dZfR1[3, gn, j1],
          {gn, MaxGenerationIndex}]} },

  C[ -F[1, {j1}], F[2, {j2}], -S[3] ] ==
    -I EL Mass[F[2, {j1}]]/(Sqrt[2] SW MW) IndexDelta[j1, j2] *
    { {0, 0},
      {1, dZe1 - dSW1/SW + dMf1[2, j1]/Mass[F[2, {j1}]] -
            dMWsq1/(2 MW^2) + dZGp1/2 +
            Conjugate[dZfL1[1, j1, j1]]/2 + dZfR1[2, j1, j1]/2} },

  C[ -F[2, {j1}], F[1, {j2}], S[3] ] ==
    -I EL Mass[F[2, {j1}]]/(Sqrt[2] SW MW) IndexDelta[j1, j2] *
    { {1, dZe1 - dSW1/SW + dMf1[2, j1]/Mass[F[2, {j1}]] -
            dMWsq1/(2 MW^2) + dZGp1/2 +
            dZfL1[1, j1, j1]/2 + Conjugate[dZfR1[2, j1, j1]]/2},
      {0, 0} },

	(* U-U-V:  G(+) . { p1_mu3, p2_mu3 } *)

  C[ -U[3], U[3], V[1] ] == -I EL *
    { {1, dZe1 + dZAA1/2 - dZW1/2 + dUW1 - CW/SW dZZA1/2},
      {0, 0} },

  C[ -U[4], U[4], V[1] ] == I EL *
    { {1, dZe1 + dZAA1/2 - dZW1/2 + dUW1 - CW/SW dZZA1/2},
      {0, 0} },

  C[ -U[3], U[3], V[2] ] == I EL CW/SW *
    { {1, dZe1 - 1/(CW^2 SW) dSW1 + dZZZ1/2 - dZW1/2 + dUW1 - SW/CW dZAZ1/2},
      {0, 0} },

  C[ -U[4], U[4], V[2] ] == -I EL CW/SW *
    { {1, dZe1 - 1/(CW^2 SW) dSW1 + dZZZ1/2 - dZW1/2 + dUW1 - SW/CW dZAZ1/2},
      {0, 0} },

  C[ -U[3], U[2], V[3] ] == -I EL CW/SW *
    { {1, dZe1 - 1/(CW^2 SW) dSW1 + dUZZ1 - SW/CW dUAZ1},
      {0, 0} },

  C[ -U[2], U[3], -V[3] ] == -I EL *
    { {CW/SW,
       CW/SW (dZe1 - 1/(CW^2 SW) dSW1 + dZW1/2 - dZZZ1/2 + dUW1) + dZZA1/2},
      {0, 0} }, 

  C[ -U[4], U[2], -V[3] ] == I EL CW/SW *
    { {1, dZe1 - 1/(CW^2 SW) dSW1 + dUZZ1 - SW/CW dUAZ1},
      {0, 0} },

  C[ -U[2], U[4], V[3] ] == I EL *
    { {CW/SW,
       CW/SW (dZe1 - 1/(CW^2 SW) dSW1 + dZW1/2 - dZZZ1/2 + dUW1) + dZZA1/2},
      {0, 0} },

  C[ -U[3], U[1], V[3] ] == I EL *
    { {1, dZe1 + dUAA1 - CW/SW dUZA1},
      {0, 0} },

  C[ -U[1], U[3], -V[3] ] == I EL *
    { {1, dZe1 + dZW1/2 - dZAA1/2 + dUW1 + CW/SW dZAZ1/2},
      {0, 0} },

  C[ -U[4], U[1], -V[3] ] == -I EL *
    { {1, dZe1 + dUAA1 - CW/SW dUZA1},
      {0, 0} },

  C[ -U[1], U[4], V[3] ] == -I EL *
    { {1, dZe1 + dZW1/2 - dZAA1/2 + dUW1 + CW/SW dZAZ1/2},
      {0, 0} },

	(* S-U-U:  G(+) . 1 *)

  C[ S[1], -U[2], U[2] ] == -I EL MZ GaugeXi[Z]/(2 SW CW) *
    { {1, dZe1 + (SW^2 - CW^2)/(CW^2 SW) dSW1 + dZH1/2 - dZG01/2 + dUZZ1} },

  C[ S[1], -U[3], U[3] ] == -I EL MW GaugeXi[W]/(2 SW) *
    { {1, dZe1 - dSW1/SW + dZH1/2 - dZGp1/2 + dUW1} },

  C[ S[1], -U[4], U[4] ] == -I EL MW GaugeXi[W]/(2 SW) *
    { {1, dZe1 - dSW1/SW + dZH1/2 - dZGp1/2 + dUW1} },

  C[ S[2], -U[4], U[4] ] == EL MW GaugeXi[W]/(2 SW) *
    { {1, dZe1 - dSW1/SW + dZG01/2 - dZGp1/2 + dUW1} },

  C[ S[2], -U[3], U[3] ] == -EL MW GaugeXi[W]/(2 SW) *
    { {1, dZe1 - dSW1/SW + dZG01/2 - dZGp1/2 + dUW1} },

  C[ -S[3], -U[2], U[3] ] == I EL MZ GaugeXi[Z]/(2 SW) *
    { {1, dZe1 - dSW1/SW + dZGp1/2 - dZG01/2 + dUW1} },

  C[ S[3], -U[2], U[4] ] == I EL MZ GaugeXi[Z]/(2 SW) *
    { {1, dZe1 + dSW1/SW + dZGp1/2 - dZG01/2 + dUW1} },

  C[ -S[3], -U[4], U[2] ] == I EL (SW^2 - CW^2) MW GaugeXi[W]/(2 CW SW) *
    { {1, dZe1 + dSW1/((SW^2 - CW^2) CW^2 SW) + dUZZ1 +
        (2 SW CW)/(SW^2 - CW^2) dUAZ1} },

  C[ S[3], -U[3], U[2] ] == I EL (SW^2 - CW^2) MW GaugeXi[W]/(2 CW SW) *
    { {1, dZe1 + dSW1/((SW^2 - CW^2) CW^2 SW) + dUZZ1 +
        (2 SW CW)/(SW^2 - CW^2) dUAZ1} },

  C[ -S[3], -U[4], U[1] ] == I EL MW GaugeXi[W] *
    { {1, dZe1 + dUAA1 + (SW^2 - CW^2)/(2 SW CW) dUZA1} },

  C[ S[3], -U[3], U[1] ] == I EL MW GaugeXi[W] *
    { {1, dZe1 + dUAA1 + (SW^2 - CW^2)/(2 SW CW) dUZA1} }
}


M$LastModelRules = {}


(* some short-hands for excluding classes of particles *)

QEDOnly = ExcludeParticles -> {F[1], V[2], V[3], S, SV, U[2], U[3], U[4]}

NoGeneration1 = ExcludeParticles -> F[_, {1, ___}]

NoGeneration2 = ExcludeParticles -> F[_, {2, ___}]

NoGeneration3 = ExcludeParticles -> F[_, {3, ___}]

NoElectronHCoupling =
  ExcludeFieldPoints -> {
    FieldPoint[_][-F[2, {1}], F[2, {1}], S],
    FieldPoint[_][-F[2, {1}], F[1, {1}], S] }

NoLightFHCoupling =
  ExcludeFieldPoints -> {
    FieldPoint[_][-F[2], F[2], S],
    FieldPoint[_][-F[2], F[1], S],
    FieldPoint[_][-F[3, {1, ___}], F[3, {1, ___}], S],
    FieldPoint[_][-F[3, {2, ___}], F[3, {2, ___}], S],
    FieldPoint[_][-F[4], F[4], S],
    FieldPoint[_][-F[4], F[3, {1, ___}], S],
    FieldPoint[_][-F[4], F[3, {2, ___}], S] }

NoQuarkMixing =
  ExcludeFieldPoints -> {
    FieldPoint[_][-F[4, {1, ___}], F[3, {2, ___}], S[3]],
    FieldPoint[_][-F[4, {1, ___}], F[3, {2, ___}], V[3]],
    FieldPoint[_][-F[4, {1, ___}], F[3, {3, ___}], S[3]],
    FieldPoint[_][-F[4, {1, ___}], F[3, {3, ___}], V[3]],
    FieldPoint[_][-F[4, {2, ___}], F[3, {1, ___}], S[3]],
    FieldPoint[_][-F[4, {2, ___}], F[3, {1, ___}], V[3]],
    FieldPoint[_][-F[4, {2, ___}], F[3, {3, ___}], S[3]],
    FieldPoint[_][-F[4, {2, ___}], F[3, {3, ___}], V[3]],
    FieldPoint[_][-F[4, {3, ___}], F[3, {1, ___}], S[3]],
    FieldPoint[_][-F[4, {3, ___}], F[3, {1, ___}], V[3]],
    FieldPoint[_][-F[4, {3, ___}], F[3, {2, ___}], S[3]],
    FieldPoint[_][-F[4, {3, ___}], F[3, {2, ___}], V[3]] }


(* The following definitions of renormalization constants
   are for the on-shell renormalization of the Standard Model in
   the scheme of A. Denner, Fortschr. d. Physik, 41 (1993) 4.

   The renormalization constants are not directly used by
   FeynArts, and hence do not restrict the generation of diagrams
   and amplitudes in any way. *)

Clear[RenConst]

RenConst[ dMf1[type_, j1_] ] := MassRC[F[type, {j1}]]

RenConst[ dZfL1[type_, j1_, j2_] ] :=
  FieldRC[F[type, {j1}], F[type, {j2}]][[1]]

RenConst[ dZfR1[type_, j1_, j2_] ] :=
  FieldRC[F[type, {j1}], F[type, {j2}]][[2]]

RenConst[ dCKM1[j1_, j2_] ] := 1/4 IndexSum[
  (dZfL1[3, j1, gn] - Conjugate[dZfL1[3, gn, j1]]) CKM[gn, j2] -
  CKM[j1, gn] (dZfL1[4, gn, j2] - Conjugate[dZfL1[4, j2, gn]]),
  {gn, MaxGenerationIndex} ]

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

RenConst[ dSW1 ] := CW^2/SW/2 (dMZsq1/MZ^2 - dMWsq1/MW^2)

RenConst[ dZe1 ] := -1/2 (dZAA1 + SW/CW dZZA1)

