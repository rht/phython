(*
	QED.mod
		Classes model file for leptons-only QED
		by Hagen Eck and Sepp Kueblbeck
		last modified 6 Mar 07 by Thomas Hahn


This file introduces the following symbols:

	coupling constants and masses:
	------------------------------
	EL:		electron charge (Thomson limit)

	MLE:		lepton class mass
	ME, MM, ML:	lepton masses (e, mu, tau)

	GaugeXi[A]:	photon gauge parameter


	one-loop renormalization constants (RCs):
	-----------------------------------------
	dZe1:		electromagnetic charge RC

	dZAA1:		photon field RC

	dMf1:		fermion mass RCs
	dZfL1, dZfR1:	fermion field RCs
*)


IndexRange[ Index[Generation] ] = {1, 2, 3}

IndexStyle[ Index[Generation, i_Integer] ] := Alph[i + 8]


ViolatesQ[ q__ ] := Plus[q] =!= 0


M$ClassesDescription = {

	(* Leptons (e, mu, tau)
	   note that in SM.mod the leptons live in class 2 (F[2]) *)
  F[1] == {
	SelfConjugate -> False,
	Indices -> {Index[Generation]},
	Mass -> MLE,
	QuantumNumbers -> -Charge,
	PropagatorLabel -> ComposedChar["e", Index[Generation]],
	PropagatorType -> Straight,
	PropagatorArrow -> Forward },

	(* Photon *)
  V[1] == {
	SelfConjugate -> True,
	Mass -> 0,
	PropagatorLabel -> "\\gamma",
	PropagatorType -> Sine,
	PropagatorArrow -> None }
}

MLE[1] = ME;
MLE[2] = MM;
MLE[3] = ML

TheLabel[ F[1, {1}] ] = "e";
TheLabel[ F[1, {2}] ] = "\\mu";
TheLabel[ F[1, {3}] ] = "\\tau"

GaugeXi[ V[1] ] = GaugeXi[A]


	(* the nomenclature has been kept compatible with SM.mod even
	   though some indices are unnecessary since everything is
	   flavour-diagonal here *)

mdZfLR1[ type_, j1_ ] :=
  Mass[F[type, j1]]/2 *
    (dZfL1[type, j1, j1] + Conjugate[dZfR1[type, j1, j1]])

mdZfRL1[ type_, j1_ ] :=
  Mass[F[type, j1]]/2 *
    (dZfR1[type, j1, j1] + Conjugate[dZfL1[type, j1, j1]])

dZfL1cc[ type_, j1_ ] :=
  dZfL1[type, j1, j1]/2 + Conjugate[dZfL1[type, j1, j1]]/2

dZfR1cc[ type_, j1_ ] :=
  dZfR1[type, j1, j1]/2 + Conjugate[dZfR1[type, j1, j1]]/2


M$CouplingMatrices = {

	(* F-F:  G(+) . { slash[mom1] omega[-], slash[mom2] omega[+],
	                  omega[-], omega[+] } *)

  C[ -F[1, {j1}], F[1, {j2}] ] == I IndexDelta[j1, j2] *
    { {0, -dZfL1cc[1, j1]},
      {0, dZfR1cc[1, j1]},
      {0, -mdZfLR1[1, j1] - dMf1[1, j1]},
      {0, -mdZfRL1[1, j1] - dMf1[1, j1]} },

	(* V-V:  G(+) . { -g[mu, nu] mom^2, g[mu, nu], -mom[mu] mom[nu] } *)

  C[ V[1], V[1] ] == I * 
    { {0, dZAA1},
      {0, 0},
      {0, -dZAA1} },

	(* F-F-V:  G(-) . { gamma[mu3] omega[-], gamma[mu3] omega[+] } *)

  C[ -F[1, {j1}], F[1, {j2}], V[1] ] == I EL IndexDelta[j1, j2] *
    { {1, dZe1 + dZAA1/2 + dZfL1cc[1, j1]},
      {1, dZe1 + dZAA1/2 + dZfR1cc[1, j1]} }
}


M$LastModelRules = {}


(* some short-hands for excluding classes of particles *)

NoGeneration1 = ExcludeParticles -> F[_, {1}]

NoGeneration2 = ExcludeParticles -> F[_, {2}]

NoGeneration3 = ExcludeParticles -> F[_, {3}]

