(*
	Toy.mod
		a toy model with a fermion, a vector boson, and
		a scalar, with fictitious f-f-v and f-f-s couplings
		last modified 26 Oct 07 th
*)


M$ClassesDescription = {
  F[1] == {
	SelfConjugate -> False,
	Mass -> mf,
	PropagatorLabel -> "f",
	PropagatorType -> Straight,
	PropagatorArrow -> Forward },
  V[1] == {
	SelfConjugate -> True,
	Mass -> mv,
	PropagatorLabel -> "v",
	PropagatorType -> Sine,
	PropagatorArrow -> None },
  S[1] == {
	SelfConjugate -> True,
	Mass -> ms,
	PropagatorLabel -> "s",
	PropagatorType -> ScalarDash,
	PropagatorArrow -> None }
}

M$CouplingMatrices = {
  C[ -F[1], F[1], V[1] ] == {{ffvL}, {ffvR}},
  C[ -F[1], F[1], S[1] ] == {{ffsL}, {ffsR}}
}

M$LastModelRules = {}

