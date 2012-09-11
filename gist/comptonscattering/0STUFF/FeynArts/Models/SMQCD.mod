(*
	(SM|THDM)QCD.mod
		Addendum classes model file for SM.mod
		to include the strong interactions
		by Christian Schappacher
		last modified 15 Jun 12 by Thomas Hahn

This file introduces the following symbols in addition to the ones in
SM.mod:

	GS: the strong coupling constant

	SUNT[a, i, j]: the generators of SU(N)
		(half the Gell-Mann matrices)

	SUNF[a, b, c]: the structure constants of SU(N)

	SUNF[a, b, c, d]: a short-hand for the sum
		\sum_i SUNF[a, b, i] SUNF[i, c, d]

	GaugeXi[g]: gluon gauge parameter

	dZGG1: gluon field RC
        dZg1: strong coupling-constant RC
*)


LoadModel[ StringReplace[$Input, "QCD" -> ""], "" ]

If[ $NoElectroweak === True, M$CouplingMatrices = {} ]


IndexRange[ Index[Gluon] ] = NoUnfold[Range[8]]

M$ClassesDescription = Join[ M$ClassesDescription, {

  V[5] == {
	SelfConjugate -> True,
	Indices -> {Index[Gluon]},
	Mass -> 0,
	PropagatorLabel -> "g",
	PropagatorType -> Cycles,
	PropagatorArrow -> None },

  U[5] == {
	SelfConjugate -> False,
	Indices -> {Index[Gluon]},
	Mass -> 0,
	QuantumNumbers -> GhostNumber,
	PropagatorLabel -> ComposedChar["u", "g"],
	PropagatorType -> GhostDash,
	PropagatorArrow -> Forward }
} ]


GaugeXi[ V[5, ___] ] = GaugeXi[g];
GaugeXi[ U[5, ___] ] = GaugeXi[g]


M$CouplingMatrices = Join[ M$CouplingMatrices, {


(*--- gluon-gluon counter term -----------------------------------------*)

  C[ V[5, {g1}], V[5, {g2}] ] == I IndexDelta[g1, g2] *
    { {0, dZGG1},
      {0, 0},
      {0, -dZGG1} },

(*--- gluon-gluon-gluon-gluon ------------------------------------------*)

  C[ V[5, {g1}], V[5, {g2}], V[5, {g3}], V[5, {g4}] ] == -I GS^2 *
    { ( SUNF[g1, g3, g2, g4] - SUNF[g1, g4, g3, g2] ) * 
      {1, 2 dZg1 + 2 dZGG1},
      ( SUNF[g1, g2, g3, g4] + SUNF[g1, g4, g3, g2] ) * 
      {1, 2 dZg1 + 2 dZGG1},
      (-SUNF[g1, g2, g3, g4] - SUNF[g1, g3, g2, g4] ) * 
      {1, 2 dZg1 + 2 dZGG1} },


(*--- gluon-gluon-gluon ------------------------------------------------*)

  C[ V[5, {g1}], V[5, {g2}], V[5, {g3}] ] == GS SUNF[g1, g2, g3] *
    { {1, dZg1 + dZGG1 + dZGG1/2} },


(*--- ghost-ghost-gluon ------------------------------------------------*)

  C[ -U[5, {g1}], U[5, {g2}], V[5, {g3}] ] == GS SUNF[g1, g2, g3] *
    { {1, dZg1 + dZGG1/2}, {0, 0} },


(*--- quark-quark-gluon ------------------------------------------------*)

  C[ -F[3, {j1, o1}], F[3, {j2, o2}], V[5, {g1}] ] == -I GS *
    SUNT[g1, o1, o2] *
    { {IndexDelta[j1, j2],
        (dZg1 + dZGG1/2) IndexDelta[j1, j2] + AddHC[dZfL1[3, j1, j2]]}, 
      {IndexDelta[j1, j2],
        (dZg1 + dZGG1/2) IndexDelta[j1, j2] + AddHC[dZfR1[3, j1, j2]]} },

  C[ -F[4, {j1, o1}], F[4, {j2, o2}], V[5, {g1}] ] == -I GS *
    SUNT[g1, o1, o2] *
    { {IndexDelta[j1, j2],
         (dZg1 + dZGG1/2) IndexDelta[j1, j2] + AddHC[dZfL1[4, j1, j2]]}, 
      {IndexDelta[j1, j2],
         (dZg1 + dZGG1/2) IndexDelta[j1, j2] + AddHC[dZfR1[4, j1, j2]]} }

} ]


RenConst[ dZGG1 ] := UVDivergentPart[FieldRC[V[5]]]

RenConst[ dZg1 ] := -1/2 dZGG1

