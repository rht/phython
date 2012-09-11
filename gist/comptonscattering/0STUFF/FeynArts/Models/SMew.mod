(*
	SMew.mod
		The Standard Model without colour indices
		last modified 3 Oct 07 th
*)


LoadModel["SM"]

M$ClassesDescription = M$ClassesDescription /.
  (Indices -> {g_, Index[Colour]}) ->
  Sequence[Indices -> {g}, MatrixTraceFactor -> 3]

M$CouplingMatrices = M$CouplingMatrices /.
  o2 -> o1 /. o1 -> Sequence[]

