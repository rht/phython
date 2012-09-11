(*
	SMbgf.mod
		Classes model file for the Standard Model in
		the background field formalism
		last modified 13 Aug 09 by Thomas Hahn

Reference:
	A. Denner, S. Dittmaier, and G. Weiglein
	Nucl. Phys. B440 (1995) 95

This file introduces the following symbols in addition to the ones in
SM.mod:

	GaugeXi[Q, bg]:		gauge parameters

V[10, 20, 30] and S[10, 20, 30] are the background fields respectively of
V[1, 2, 3] and S[1, 2, 3].

*)


IndexRange[ Index[Generation] ] = Range[3]

IndexRange[ Index[Colour] ] = NoUnfold[Range[3]]

IndexStyle[ Index[Generation, i_Integer] ] := Alph[i + 8]

MaxGenerationIndex = 3


ViolatesQ[ q__ ] := Plus[q] =!= 0


(* the leptonic field RCs are diagonal: *)

dZfL1[ type:1 | 2, j1_, j2_ ] :=
  IndexDelta[j1, j2] dZfL1[type, j1, j1] /; j1 =!= j2

dZfR1[ type:1 | 2, j1_, j2_ ] :=
  IndexDelta[j1, j2] dZfR1[type, j1, j1] /; j1 =!= j2



M$ClassesDescription =
{F[1] == {SelfConjugate -> False, Indices -> {Index[Generation]}, Mass -> 0, 
   QuantumNumbers -> LeptonNumber, PropagatorLabel -> 
    ComposedChar["\\nu", Index[Generation]], PropagatorType -> Straight, 
   PropagatorArrow -> Forward}, F[2] == {SelfConjugate -> False, 
   Indices -> {Index[Generation]}, Mass -> MLE, 
   QuantumNumbers -> {-Charge, LeptonNumber}, PropagatorLabel -> 
    ComposedChar["e", Index[Generation]], PropagatorType -> Straight, 
   PropagatorArrow -> Forward}, F[3] == {SelfConjugate -> False, 
   Indices -> {Index[Generation], Index[Colour]}, Mass -> MQU, 
   QuantumNumbers -> (2*Charge)/3, PropagatorLabel -> 
    ComposedChar["u", Index[Generation]], PropagatorType -> Straight, 
   PropagatorArrow -> Forward}, F[4] == {SelfConjugate -> False, 
   Indices -> {Index[Generation], Index[Colour]}, Mass -> MQD, 
   QuantumNumbers -> -Charge/3, PropagatorLabel -> 
    ComposedChar["d", Index[Generation]], PropagatorType -> Straight, 
   PropagatorArrow -> Forward}, V[1] == {SelfConjugate -> True, 
   Indices -> {}, Mass -> 0, PropagatorLabel -> "\\gamma", 
   PropagatorType -> Sine, PropagatorArrow -> None, InsertOnly -> Loop}, 
 V[10] == {SelfConjugate -> True, Indices -> {}, Mass -> 0, 
   PropagatorLabel -> ComposedChar["\\gamma", Null, Null, "\\hat"], 
   PropagatorType -> Sine, PropagatorArrow -> None, 
   InsertOnly -> {Internal, External}}, 
 V[2] == {SelfConjugate -> True, Indices -> {}, Mass -> MZ, 
   PropagatorLabel -> "Z", PropagatorType -> Sine, PropagatorArrow -> None, 
   InsertOnly -> Loop}, V[20] == {SelfConjugate -> True, Indices -> {}, 
   Mass -> MZ, PropagatorLabel -> ComposedChar["Z", Null, Null, "\\hat"], 
   PropagatorType -> Sine, PropagatorArrow -> None, 
   InsertOnly -> {Internal, External}}, 
 V[3] == {SelfConjugate -> False, Indices -> {}, Mass -> MW, 
   QuantumNumbers -> -Charge, PropagatorLabel -> "W", PropagatorType -> Sine, 
   PropagatorArrow -> Forward, InsertOnly -> Loop}, 
 V[30] == {SelfConjugate -> False, Indices -> {}, Mass -> MW, 
   QuantumNumbers -> -Charge, PropagatorLabel -> ComposedChar["W", Null, 
     Null, "\\hat"], PropagatorType -> Sine, PropagatorArrow -> Forward, 
   InsertOnly -> {Internal, External}}, 
 SV[2] == {SelfConjugate -> True, Indices -> {}, Mass -> MZ, 
   MixingPartners -> {S[2], V[2]}, PropagatorLabel -> 
    {ComposedChar["G", Null, "0"], "Z"}, PropagatorType -> 
    {ScalarDash, Sine}, PropagatorArrow -> None, InsertOnly -> Loop}, 
 SV[20] == {SelfConjugate -> True, Indices -> {}, Mass -> MZ, 
   MixingPartners -> {S[2], V[2]}, PropagatorLabel -> 
    {ComposedChar["G", Null, "0", "\\hat"], ComposedChar["Z", Null, Null, 
      "\\hat"]}, PropagatorType -> {ScalarDash, Sine}, 
   PropagatorArrow -> None, InsertOnly -> {Internal, External}}, 
 SV[3] == {SelfConjugate -> False, Indices -> {}, Mass -> MW, 
   QuantumNumbers -> -Charge, MixingPartners -> {S[3], V[3]}, 
   PropagatorLabel -> {"G", "W"}, PropagatorType -> {ScalarDash, Sine}, 
   PropagatorArrow -> Forward, InsertOnly -> Loop}, 
 SV[30] == {SelfConjugate -> False, Indices -> {}, Mass -> MW, 
   QuantumNumbers -> -Charge, MixingPartners -> {S[3], V[3]}, 
   PropagatorLabel -> {ComposedChar["G", Null, Null, "\\hat"], 
     ComposedChar["W", Null, Null, "\\hat"]}, 
   PropagatorType -> {ScalarDash, Sine}, PropagatorArrow -> Forward, 
   InsertOnly -> {Internal, External}}, 
 S[1] == {SelfConjugate -> True, Indices -> {}, Mass -> MH, 
   PropagatorLabel -> "H", PropagatorType -> ScalarDash, 
   PropagatorArrow -> None, InsertOnly -> Loop}, 
 S[10] == {SelfConjugate -> True, Indices -> {}, Mass -> MH, 
   PropagatorLabel -> ComposedChar["H", Null, Null, "\\hat"], 
   PropagatorType -> ScalarDash, PropagatorArrow -> None, 
   InsertOnly -> {Internal, External}}, 
 S[2] == {SelfConjugate -> True, Indices -> {}, Mass -> MZ, 
   PropagatorLabel -> ComposedChar["G", Null, "0"], 
   PropagatorType -> ScalarDash, PropagatorArrow -> None, 
   InsertOnly -> Loop}, S[20] == {SelfConjugate -> True, Indices -> {}, 
   Mass -> MZ, PropagatorLabel -> ComposedChar["G", Null, "0", "\\hat"], 
   PropagatorType -> ScalarDash, PropagatorArrow -> None, 
   InsertOnly -> {Internal, External}}, 
 S[3] == {SelfConjugate -> False, Indices -> {}, Mass -> MW, 
   QuantumNumbers -> -Charge, PropagatorLabel -> "G", 
   PropagatorType -> ScalarDash, PropagatorArrow -> Forward, 
   InsertOnly -> Loop}, S[30] == {SelfConjugate -> False, Indices -> {}, 
   Mass -> MW, QuantumNumbers -> -Charge, PropagatorLabel -> 
    ComposedChar["G", Null, Null, "\\hat"], PropagatorType -> ScalarDash, 
   PropagatorArrow -> Forward, InsertOnly -> {Internal, External}}, 
 U[1] == {SelfConjugate -> False, Indices -> {}, Mass -> 0, 
   QuantumNumbers -> GhostNumber, PropagatorLabel -> 
    ComposedChar["u", "\\gamma"], PropagatorType -> GhostDash, 
   PropagatorArrow -> Forward}, U[2] == {SelfConjugate -> False, 
   Indices -> {}, Mass -> MZ, QuantumNumbers -> GhostNumber, 
   PropagatorLabel -> ComposedChar["u", "Z"], PropagatorType -> GhostDash, 
   PropagatorArrow -> Forward}, U[3] == {SelfConjugate -> False, 
   Indices -> {}, Mass -> MW, QuantumNumbers -> {-Charge, GhostNumber}, 
   PropagatorLabel -> ComposedChar["u", "-"], PropagatorType -> GhostDash, 
   PropagatorArrow -> Forward}, U[4] == {SelfConjugate -> False, 
   Indices -> {}, Mass -> MW, QuantumNumbers -> {Charge, GhostNumber}, 
   PropagatorLabel -> ComposedChar["u", "+"], PropagatorType -> GhostDash, 
   PropagatorArrow -> Forward}}

M$CouplingMatrices =
{C[-V[3], -V[3], V[3], V[3]] == {{((2*I)*EL^2)/SW^2}, {((-I)*EL^2)/SW^2}, 
   {((-I)*EL^2)/SW^2}}, C[-V[3], V[3], V[2], V[2]] == 
  {{((-2*I)*CW^2*EL^2)/SW^2}, {(I*CW^2*EL^2)/SW^2}, {(I*CW^2*EL^2)/SW^2}}, 
 C[-V[3], V[3], V[1], V[2]] == {{((2*I)*CW*EL^2)/SW}, {((-I)*CW*EL^2)/SW}, 
   {((-I)*CW*EL^2)/SW}}, C[-V[3], V[3], V[1], V[1]] == 
  {{(-2*I)*EL^2}, {I*EL^2}, {I*EL^2}}, C[V[1], -V[3], V[3]] == 
  {{(-I)*EL}, {0}, {0}, {0}}, C[V[2], -V[3], V[3]] == 
  {{(I*CW*EL)/SW}, {0}, {0}, {0}}, C[S[1], S[1], S[1], S[1]] == 
  {{(((-3*I)/4)*EL^2*MH^2)/(MW^2*SW^2)}}, C[S[1], S[1], S[2], S[2]] == 
  {{((-I/4)*EL^2*MH^2)/(MW^2*SW^2)}}, C[S[1], S[1], S[3], -S[3]] == 
  {{((-I/4)*EL^2*MH^2)/(MW^2*SW^2)}}, C[S[2], S[2], S[2], S[2]] == 
  {{(((-3*I)/4)*EL^2*MH^2)/(MW^2*SW^2)}}, C[S[2], S[2], S[3], -S[3]] == 
  {{((-I/4)*EL^2*MH^2)/(MW^2*SW^2)}}, C[S[3], S[3], -S[3], -S[3]] == 
  {{((-I/2)*EL^2*MH^2)/(MW^2*SW^2)}}, C[S[1], S[1], S[1]] == 
  {{(((-3*I)/2)*EL*MH^2)/(MW*SW)}}, C[S[1], S[2], S[2]] == 
  {{((-I/2)*EL*MH^2)/(MW*SW)}}, C[S[3], S[1], -S[3]] == 
  {{((-I/2)*EL*MH^2)/(MW*SW)}}, C[S[1], S[1], V[3], -V[3]] == 
  {{((I/2)*EL^2)/SW^2}}, C[S[2], S[2], V[3], -V[3]] == {{((I/2)*EL^2)/SW^2}}, 
 C[S[3], -S[3], V[3], -V[3]] == {{((I/2)*EL^2)/SW^2}}, 
 C[S[3], -S[3], V[2], V[2]] == {{((I/2)*EL^2*(-CW^2 + SW^2)^2)/(CW^2*SW^2)}}, 
 C[S[3], -S[3], V[1], V[2]] == {{(I*EL^2*(-CW^2 + SW^2))/(CW*SW)}}, 
 C[S[3], -S[3], V[1], V[1]] == {{(2*I)*EL^2}}, 
 C[S[1], S[1], V[2], V[2]] == {{((I/2)*EL^2)/(CW^2*SW^2)}}, 
 C[S[2], S[2], V[2], V[2]] == {{((I/2)*EL^2)/(CW^2*SW^2)}}, 
 C[S[1], -S[3], V[3], V[2]] == {{((-I/2)*EL^2)/CW}}, 
 C[S[1], S[3], -V[3], V[2]] == {{((-I/2)*EL^2)/CW}}, 
 C[S[1], S[3], -V[3], V[1]] == {{((-I/2)*EL^2)/SW}}, 
 C[S[1], -S[3], V[3], V[1]] == {{((-I/2)*EL^2)/SW}}, 
 C[S[3], S[2], V[2], -V[3]] == {{EL^2/(2*CW)}}, 
 C[-S[3], S[2], V[2], V[3]] == {{-EL^2/(2*CW)}}, 
 C[S[3], S[2], V[1], -V[3]] == {{EL^2/(2*SW)}}, 
 C[-S[3], S[2], V[1], V[3]] == {{-EL^2/(2*SW)}}, 
 C[S[2], S[1], V[2]] == {{EL/(2*CW*SW)}, {-EL/(2*CW*SW)}}, 
 C[-S[3], S[3], V[1]] == {{(-I)*EL}, {I*EL}}, 
 C[-S[3], S[3], V[2]] == {{((-I/2)*EL*(-CW^2 + SW^2))/(CW*SW)}, 
   {((I/2)*EL*(-CW^2 + SW^2))/(CW*SW)}}, C[S[3], S[1], -V[3]] == 
  {{((-I/2)*EL)/SW}, {((I/2)*EL)/SW}}, C[-S[3], S[1], V[3]] == 
  {{((I/2)*EL)/SW}, {((-I/2)*EL)/SW}}, C[S[3], S[2], -V[3]] == 
  {{EL/(2*SW)}, {-EL/(2*SW)}}, C[-S[3], S[2], V[3]] == 
  {{EL/(2*SW)}, {-EL/(2*SW)}}, C[S[1], -V[3], V[3]] == {{(I*EL*MW)/SW}}, 
 C[S[1], V[2], V[2]] == {{(I*EL*MW)/(CW^2*SW)}}, 
 C[-S[3], V[3], V[2]] == {{((-I)*EL*MW*SW)/CW}}, 
 C[S[3], -V[3], V[2]] == {{((-I)*EL*MW*SW)/CW}}, 
 C[-S[3], V[3], V[1]] == {{(-I)*EL*MW}}, C[S[3], -V[3], V[1]] == 
  {{(-I)*EL*MW}}, C[-F[2, {j1}], F[2, {j2}], V[1]] == 
  {{I*EL*IndexDelta[j1, j2]}, {I*EL*IndexDelta[j1, j2]}}, 
 C[-F[3, {j1, o1}], F[3, {j2, o2}], V[1]] == 
  {{((-2*I)/3)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]}, 
   {((-2*I)/3)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]}}, 
 C[-F[4, {j1, o1}], F[4, {j2, o2}], V[1]] == 
  {{(I/3)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]}, 
   {(I/3)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]}}, 
 C[-F[1, {j1}], F[1, {j2}], V[2]] == 
  {{((I/2)*EL*IndexDelta[j1, j2])/(CW*SW)}, {0}}, 
 C[-F[2, {j1}], F[2, {j2}], V[2]] == 
  {{(I*EL*(-1/2 + SW^2)*IndexDelta[j1, j2])/(CW*SW)}, 
   {(I*EL*SW*IndexDelta[j1, j2])/CW}}, 
 C[-F[3, {j1, o1}], F[3, {j2, o2}], V[2]] == 
  {{(I*EL*(1/2 - (2*SW^2)/3)*IndexDelta[j1, j2]*IndexDelta[o1, o2])/(CW*SW)}, 
   {(((-2*I)/3)*EL*SW*IndexDelta[j1, j2]*IndexDelta[o1, o2])/CW}}, 
 C[-F[4, {j1, o1}], F[4, {j2, o2}], V[2]] == 
  {{(I*EL*(-1/2 + SW^2/3)*IndexDelta[j1, j2]*IndexDelta[o1, o2])/(CW*SW)}, 
   {((I/3)*EL*SW*IndexDelta[j1, j2]*IndexDelta[o1, o2])/CW}}, 
 C[-F[1, {j1}], F[2, {j2}], -V[3]] == 
  {{(I*EL*IndexDelta[j1, j2])/(Sqrt[2]*SW)}, {0}}, 
 C[-F[2, {j1}], F[1, {j2}], V[3]] == 
  {{(I*EL*IndexDelta[j1, j2])/(Sqrt[2]*SW)}, {0}}, 
 C[-F[3, {j1, o1}], F[4, {j2, o2}], -V[3]] == 
  {{(I*EL*CKM[j1, j2]*IndexDelta[o1, o2])/(Sqrt[2]*SW)}, {0}}, 
 C[-F[4, {j2, o2}], F[3, {j1, o1}], V[3]] == 
  {{(I*EL*Conjugate[CKM[j1, j2]]*IndexDelta[o1, o2])/(Sqrt[2]*SW)}, {0}}, 
 C[-F[2, {j1}], F[2, {j2}], S[1]] == 
  {{((-I/2)*EL*IndexDelta[j1, j2]*Mass[F[2, {j1}]])/(MW*SW)}, 
   {((-I/2)*EL*IndexDelta[j1, j2]*Mass[F[2, {j1}]])/(MW*SW)}}, 
 C[-F[3, {j1, o1}], F[3, {j2, o2}], S[1]] == 
  {{((-I/2)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[3, {j1}]])/
     (MW*SW)}, {((-I/2)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      Mass[F[3, {j1}]])/(MW*SW)}}, 
 C[-F[4, {j1, o1}], F[4, {j2, o2}], S[1]] == 
  {{((-I/2)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[4, {j1}]])/
     (MW*SW)}, {((-I/2)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*
      Mass[F[4, {j1}]])/(MW*SW)}}, C[-F[2, {j1}], F[2, {j2}], S[2]] == 
  {{-(EL*IndexDelta[j1, j2]*Mass[F[2, {j1}]])/(2*MW*SW)}, 
   {(EL*IndexDelta[j1, j2]*Mass[F[2, {j1}]])/(2*MW*SW)}}, 
 C[-F[3, {j1, o1}], F[3, {j2, o2}], S[2]] == 
  {{(EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[3, {j1}]])/(2*MW*SW)}, 
   {-(EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[3, {j1}]])/(2*MW*SW)}}, 
 C[-F[4, {j1, o1}], F[4, {j2, o2}], S[2]] == 
  {{-(EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[4, {j1}]])/(2*MW*SW)}, 
   {(EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[4, {j1}]])/(2*MW*SW)}}, 
 C[-F[3, {j1, o1}], F[4, {j2, o2}], -S[3]] == 
  {{(I*EL*CKM[j1, j2]*IndexDelta[o1, o2]*Mass[F[3, {j1}]])/(Sqrt[2]*MW*SW)}, 
   {((-I)*EL*CKM[j1, j2]*IndexDelta[o1, o2]*Mass[F[4, {j2}]])/
     (Sqrt[2]*MW*SW)}}, C[-F[4, {j2, o2}], F[3, {j1, o1}], S[3]] == 
  {{((-I)*EL*Conjugate[CKM[j1, j2]]*IndexDelta[o1, o2]*Mass[F[4, {j2}]])/
     (Sqrt[2]*MW*SW)}, {(I*EL*Conjugate[CKM[j1, j2]]*IndexDelta[o1, o2]*
      Mass[F[3, {j1}]])/(Sqrt[2]*MW*SW)}}, 
 C[-F[1, {j1}], F[2, {j2}], -S[3]] == 
  {{0}, {((-I)*EL*IndexDelta[j1, j2]*Mass[F[2, {j1}]])/(Sqrt[2]*MW*SW)}}, 
 C[-F[2, {j1}], F[1, {j2}], S[3]] == 
  {{((-I)*EL*IndexDelta[j1, j2]*Mass[F[2, {j1}]])/(Sqrt[2]*MW*SW)}, {0}}, 
 C[-U[3], U[3], V[1]] == {{(-I)*EL}, {0}}, C[-U[4], U[4], V[1]] == 
  {{I*EL}, {0}}, C[-U[3], U[3], V[2]] == {{(I*CW*EL)/SW}, {0}}, 
 C[-U[4], U[4], V[2]] == {{((-I)*CW*EL)/SW}, {0}}, 
 C[-U[3], U[2], V[3]] == {{((-I)*CW*EL)/SW}, {0}}, 
 C[-U[2], U[3], -V[3]] == {{((-I)*CW*EL)/SW}, {0}}, 
 C[-U[4], U[2], -V[3]] == {{(I*CW*EL)/SW}, {0}}, 
 C[-U[2], U[4], V[3]] == {{(I*CW*EL)/SW}, {0}}, 
 C[-U[3], U[1], V[3]] == {{I*EL}, {0}}, C[-U[1], U[3], -V[3]] == 
  {{I*EL}, {0}}, C[-U[4], U[1], -V[3]] == {{(-I)*EL}, {0}}, 
 C[-U[1], U[4], V[3]] == {{(-I)*EL}, {0}}, C[S[1], -U[2], U[2]] == 
  {{((-I/2)*EL*MZ*GaugeXi[Z])/(CW*SW)}}, C[S[1], -U[3], U[3]] == 
  {{((-I/2)*EL*MW*GaugeXi[W])/SW}}, C[S[1], -U[4], U[4]] == 
  {{((-I/2)*EL*MW*GaugeXi[W])/SW}}, C[S[2], -U[4], U[4]] == 
  {{(EL*MW*GaugeXi[W])/(2*SW)}}, C[S[2], -U[3], U[3]] == 
  {{-(EL*MW*GaugeXi[W])/(2*SW)}}, C[-S[3], -U[2], U[3]] == 
  {{((I/2)*EL*MZ*GaugeXi[Z])/SW}}, C[S[3], -U[2], U[4]] == 
  {{((I/2)*EL*MZ*GaugeXi[Z])/SW}}, C[-S[3], -U[4], U[2]] == 
  {{((I/2)*EL*MW*(-CW^2 + SW^2)*GaugeXi[W])/(CW*SW)}}, 
 C[S[3], -U[3], U[2]] == {{((I/2)*EL*MW*(-CW^2 + SW^2)*GaugeXi[W])/(CW*SW)}}, 
 C[-S[3], -U[4], U[1]] == {{I*EL*MW*GaugeXi[W]}}, 
 C[S[3], -U[3], U[1]] == {{I*EL*MW*GaugeXi[W]}}, 
 C[-F[1, {j1}], F[1, {j2}]] == 
  {{0, (-I/2)*(Conjugate[dZfL1[1, j1, j1]] + dZfL1[1, j1, j1])*
     IndexDelta[j1, j2]}, {0, (I/2)*(Conjugate[dZfR1[1, j1, j1]] + 
      dZfR1[1, j1, j1])*IndexDelta[j1, j2]}, {0, 0}, {0, 0}}, 
 C[-F[2, {j1}], F[2, {j2}]] == 
  {{0, (-I/2)*(Conjugate[dZfL1[2, j1, j1]] + dZfL1[2, j1, j1])*
     IndexDelta[j1, j2]}, {0, (I/2)*(Conjugate[dZfR1[2, j1, j1]] + 
      dZfR1[2, j1, j1])*IndexDelta[j1, j2]}, 
   {0, (-I/2)*IndexDelta[j1, j2]*(2*dMf1[2, j1] + dZfL1[2, j1, j1]*
       Mass[F[2, {j1}]] + Conjugate[dZfR1[2, j2, j2]]*Mass[F[2, {j2}]])}, 
   {0, (-I/2)*IndexDelta[j1, j2]*(2*dMf1[2, j1] + dZfR1[2, j1, j1]*
       Mass[F[2, {j1}]] + Conjugate[dZfL1[2, j2, j2]]*Mass[F[2, {j2}]])}}, 
 C[-F[3, {j1, o1}], F[3, {j2, o2}]] == 
  {{0, (-I/2)*(Conjugate[dZfL1[3, j2, j1]] + dZfL1[3, j1, j2])*
     IndexDelta[o1, o2]}, {0, (I/2)*(Conjugate[dZfR1[3, j2, j1]] + 
      dZfR1[3, j1, j2])*IndexDelta[o1, o2]}, 
   {0, (-I/2)*IndexDelta[o1, o2]*(2*dMf1[3, j1]*IndexDelta[j1, j2] + 
      dZfL1[3, j1, j2]*Mass[F[3, {j1}]] + Conjugate[dZfR1[3, j2, j1]]*
       Mass[F[3, {j2}]])}, {0, (-I/2)*IndexDelta[o1, o2]*
     (2*dMf1[3, j1]*IndexDelta[j1, j2] + dZfR1[3, j1, j2]*Mass[F[3, {j1}]] + 
      Conjugate[dZfL1[3, j2, j1]]*Mass[F[3, {j2}]])}}, 
 C[-F[4, {j1, o1}], F[4, {j2, o2}]] == 
  {{0, (-I/2)*(Conjugate[dZfL1[4, j2, j1]] + dZfL1[4, j1, j2])*
     IndexDelta[o1, o2]}, {0, (I/2)*(Conjugate[dZfR1[4, j2, j1]] + 
      dZfR1[4, j1, j2])*IndexDelta[o1, o2]}, 
   {0, (-I/2)*IndexDelta[o1, o2]*(2*dMf1[4, j1]*IndexDelta[j1, j2] + 
      dZfL1[4, j1, j2]*Mass[F[4, {j1}]] + Conjugate[dZfR1[4, j2, j1]]*
       Mass[F[4, {j2}]])}, {0, (-I/2)*IndexDelta[o1, o2]*
     (2*dMf1[4, j1]*IndexDelta[j1, j2] + dZfR1[4, j1, j2]*Mass[F[4, {j1}]] + 
      Conjugate[dZfL1[4, j2, j1]]*Mass[F[4, {j2}]])}}, 
 C[-V[30], -V[30], V[30], V[30]] == 
  {{((2*I)*EL^2)/SW^2, ((2*I)*EL^2*(-dMWsq1 + dZH1*MW^2))/(MW^2*SW^2)}, 
   {((-I)*EL^2)/SW^2, (I*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW^2)}, 
   {((-I)*EL^2)/SW^2, (I*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW^2)}}, 
 C[-V[30], -V[30], V[30], V[3]] == 
  {{((2*I)*EL^2)/SW^2, ((2*I)*EL^2*(-dMWsq1 + dZH1*MW^2))/(MW^2*SW^2)}, 
   {((-I)*EL^2)/SW^2, (I*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW^2)}, 
   {((-I)*EL^2)/SW^2, (I*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW^2)}}, 
 C[-V[30], -V[3], V[30], V[30]] == 
  {{((2*I)*EL^2)/SW^2, ((2*I)*EL^2*(-dMWsq1 + dZH1*MW^2))/(MW^2*SW^2)}, 
   {((-I)*EL^2)/SW^2, (I*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW^2)}, 
   {((-I)*EL^2)/SW^2, (I*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW^2)}}, 
 C[-V[30], -V[3], V[3], V[3]] == 
  {{((2*I)*EL^2)/SW^2, ((2*I)*EL^2*(-dMWsq1 + dZH1*MW^2))/(MW^2*SW^2)}, 
   {((-I)*EL^2)/SW^2, (I*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW^2)}, 
   {((-I)*EL^2)/SW^2, (I*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW^2)}}, 
 C[-V[3], -V[3], V[30], V[3]] == 
  {{((2*I)*EL^2)/SW^2, ((2*I)*EL^2*(-dMWsq1 + dZH1*MW^2))/(MW^2*SW^2)}, 
   {((-I)*EL^2)/SW^2, (I*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW^2)}, 
   {((-I)*EL^2)/SW^2, (I*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW^2)}}, 
 C[-V[30], V[30], V[20], V[20]] == 
  {{((-2*I)*CW^2*EL^2)/SW^2, ((2*I)*CW^2*EL^2*(dMWsq1 - dZH1*MW^2))/
     (MW^2*SW^2)}, {(I*CW^2*EL^2)/SW^2, (I*CW^2*EL^2*(-dMWsq1 + dZH1*MW^2))/
     (MW^2*SW^2)}, {(I*CW^2*EL^2)/SW^2, (I*CW^2*EL^2*(-dMWsq1 + dZH1*MW^2))/
     (MW^2*SW^2)}}, C[-V[30], V[30], V[20], V[2]] == 
  {{((-2*I)*CW^2*EL^2)/SW^2, ((2*I)*CW^2*EL^2*(dMWsq1 - dZH1*MW^2))/
     (MW^2*SW^2)}, {(I*CW^2*EL^2)/SW^2, (I*CW^2*EL^2*(-dMWsq1 + dZH1*MW^2))/
     (MW^2*SW^2)}, {(I*CW^2*EL^2)/SW^2, (I*CW^2*EL^2*(-dMWsq1 + dZH1*MW^2))/
     (MW^2*SW^2)}}, C[-V[30], V[3], V[20], V[20]] == 
  {{((-2*I)*CW^2*EL^2)/SW^2, ((2*I)*CW^2*EL^2*(dMWsq1 - dZH1*MW^2))/
     (MW^2*SW^2)}, {(I*CW^2*EL^2)/SW^2, (I*CW^2*EL^2*(-dMWsq1 + dZH1*MW^2))/
     (MW^2*SW^2)}, {(I*CW^2*EL^2)/SW^2, (I*CW^2*EL^2*(-dMWsq1 + dZH1*MW^2))/
     (MW^2*SW^2)}}, C[-V[3], V[30], V[20], V[20]] == 
  {{((-2*I)*CW^2*EL^2)/SW^2, ((2*I)*CW^2*EL^2*(dMWsq1 - dZH1*MW^2))/
     (MW^2*SW^2)}, {(I*CW^2*EL^2)/SW^2, (I*CW^2*EL^2*(-dMWsq1 + dZH1*MW^2))/
     (MW^2*SW^2)}, {(I*CW^2*EL^2)/SW^2, (I*CW^2*EL^2*(-dMWsq1 + dZH1*MW^2))/
     (MW^2*SW^2)}}, C[-V[30], V[3], V[2], V[2]] == 
  {{((-2*I)*CW^2*EL^2)/SW^2, ((2*I)*CW^2*EL^2*(dMWsq1 - dZH1*MW^2))/
     (MW^2*SW^2)}, {(I*CW^2*EL^2)/SW^2, (I*CW^2*EL^2*(-dMWsq1 + dZH1*MW^2))/
     (MW^2*SW^2)}, {(I*CW^2*EL^2)/SW^2, (I*CW^2*EL^2*(-dMWsq1 + dZH1*MW^2))/
     (MW^2*SW^2)}}, C[-V[3], V[30], V[2], V[2]] == 
  {{((-2*I)*CW^2*EL^2)/SW^2, ((2*I)*CW^2*EL^2*(dMWsq1 - dZH1*MW^2))/
     (MW^2*SW^2)}, {(I*CW^2*EL^2)/SW^2, (I*CW^2*EL^2*(-dMWsq1 + dZH1*MW^2))/
     (MW^2*SW^2)}, {(I*CW^2*EL^2)/SW^2, (I*CW^2*EL^2*(-dMWsq1 + dZH1*MW^2))/
     (MW^2*SW^2)}}, C[-V[3], V[3], V[20], V[2]] == 
  {{((-2*I)*CW^2*EL^2)/SW^2, ((2*I)*CW^2*EL^2*(dMWsq1 - dZH1*MW^2))/
     (MW^2*SW^2)}, {(I*CW^2*EL^2)/SW^2, (I*CW^2*EL^2*(-dMWsq1 + dZH1*MW^2))/
     (MW^2*SW^2)}, {(I*CW^2*EL^2)/SW^2, (I*CW^2*EL^2*(-dMWsq1 + dZH1*MW^2))/
     (MW^2*SW^2)}}, C[-V[30], V[30], V[10], V[20]] == 
  {{((2*I)*CW*EL^2)/SW, ((2*I)*CW*EL^2*(-dMWsq1 + dZH1*MW^2))/(MW^2*SW)}, 
   {((-I)*CW*EL^2)/SW, (I*CW*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW)}, 
   {((-I)*CW*EL^2)/SW, (I*CW*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW)}}, 
 C[-V[30], V[30], V[10], V[2]] == 
  {{((2*I)*CW*EL^2)/SW, ((2*I)*CW*EL^2*(-dMWsq1 + dZH1*MW^2))/(MW^2*SW)}, 
   {((-I)*CW*EL^2)/SW, (I*CW*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW)}, 
   {((-I)*CW*EL^2)/SW, (I*CW*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW)}}, 
 C[-V[30], V[30], V[1], V[20]] == 
  {{((2*I)*CW*EL^2)/SW, ((2*I)*CW*EL^2*(-dMWsq1 + dZH1*MW^2))/(MW^2*SW)}, 
   {((-I)*CW*EL^2)/SW, (I*CW*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW)}, 
   {((-I)*CW*EL^2)/SW, (I*CW*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW)}}, 
 C[-V[30], V[3], V[10], V[20]] == 
  {{((2*I)*CW*EL^2)/SW, ((2*I)*CW*EL^2*(-dMWsq1 + dZH1*MW^2))/(MW^2*SW)}, 
   {((-I)*CW*EL^2)/SW, (I*CW*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW)}, 
   {((-I)*CW*EL^2)/SW, (I*CW*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW)}}, 
 C[-V[3], V[30], V[10], V[20]] == 
  {{((2*I)*CW*EL^2)/SW, ((2*I)*CW*EL^2*(-dMWsq1 + dZH1*MW^2))/(MW^2*SW)}, 
   {((-I)*CW*EL^2)/SW, (I*CW*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW)}, 
   {((-I)*CW*EL^2)/SW, (I*CW*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW)}}, 
 C[-V[30], V[3], V[1], V[2]] == 
  {{((2*I)*CW*EL^2)/SW, ((2*I)*CW*EL^2*(-dMWsq1 + dZH1*MW^2))/(MW^2*SW)}, 
   {((-I)*CW*EL^2)/SW, (I*CW*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW)}, 
   {((-I)*CW*EL^2)/SW, (I*CW*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW)}}, 
 C[-V[3], V[30], V[1], V[2]] == 
  {{((2*I)*CW*EL^2)/SW, ((2*I)*CW*EL^2*(-dMWsq1 + dZH1*MW^2))/(MW^2*SW)}, 
   {((-I)*CW*EL^2)/SW, (I*CW*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW)}, 
   {((-I)*CW*EL^2)/SW, (I*CW*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW)}}, 
 C[-V[3], V[3], V[10], V[2]] == 
  {{((2*I)*CW*EL^2)/SW, ((2*I)*CW*EL^2*(-dMWsq1 + dZH1*MW^2))/(MW^2*SW)}, 
   {((-I)*CW*EL^2)/SW, (I*CW*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW)}, 
   {((-I)*CW*EL^2)/SW, (I*CW*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW)}}, 
 C[-V[3], V[3], V[1], V[20]] == 
  {{((2*I)*CW*EL^2)/SW, ((2*I)*CW*EL^2*(-dMWsq1 + dZH1*MW^2))/(MW^2*SW)}, 
   {((-I)*CW*EL^2)/SW, (I*CW*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW)}, 
   {((-I)*CW*EL^2)/SW, (I*CW*EL^2*(dMWsq1 - dZH1*MW^2))/(MW^2*SW)}}, 
 C[-V[30], V[30], V[10], V[10]] == 
  {{(-2*I)*EL^2, (-2*I)*EL^2*(dZH1 - dMWsq1/MW^2)}, 
   {I*EL^2, (-I)*EL^2*(-dZH1 + dMWsq1/MW^2)}, 
   {I*EL^2, (-I)*EL^2*(-dZH1 + dMWsq1/MW^2)}}, 
 C[-V[30], V[30], V[10], V[1]] == 
  {{(-2*I)*EL^2, (-2*I)*EL^2*(dZH1 - dMWsq1/MW^2)}, 
   {I*EL^2, (-I)*EL^2*(-dZH1 + dMWsq1/MW^2)}, 
   {I*EL^2, (-I)*EL^2*(-dZH1 + dMWsq1/MW^2)}}, 
 C[-V[30], V[3], V[10], V[10]] == 
  {{(-2*I)*EL^2, (-2*I)*EL^2*(dZH1 - dMWsq1/MW^2)}, 
   {I*EL^2, (-I)*EL^2*(-dZH1 + dMWsq1/MW^2)}, 
   {I*EL^2, (-I)*EL^2*(-dZH1 + dMWsq1/MW^2)}}, 
 C[-V[3], V[30], V[10], V[10]] == 
  {{(-2*I)*EL^2, (-2*I)*EL^2*(dZH1 - dMWsq1/MW^2)}, 
   {I*EL^2, (-I)*EL^2*(-dZH1 + dMWsq1/MW^2)}, 
   {I*EL^2, (-I)*EL^2*(-dZH1 + dMWsq1/MW^2)}}, 
 C[-V[30], V[3], V[1], V[1]] == 
  {{(-2*I)*EL^2, (-2*I)*EL^2*(dZH1 - dMWsq1/MW^2)}, 
   {I*EL^2, (-I)*EL^2*(-dZH1 + dMWsq1/MW^2)}, 
   {I*EL^2, (-I)*EL^2*(-dZH1 + dMWsq1/MW^2)}}, 
 C[-V[3], V[30], V[1], V[1]] == 
  {{(-2*I)*EL^2, (-2*I)*EL^2*(dZH1 - dMWsq1/MW^2)}, 
   {I*EL^2, (-I)*EL^2*(-dZH1 + dMWsq1/MW^2)}, 
   {I*EL^2, (-I)*EL^2*(-dZH1 + dMWsq1/MW^2)}}, 
 C[-V[3], V[3], V[10], V[1]] == 
  {{(-2*I)*EL^2, (-2*I)*EL^2*(dZH1 - dMWsq1/MW^2)}, 
   {I*EL^2, (-I)*EL^2*(-dZH1 + dMWsq1/MW^2)}, 
   {I*EL^2, (-I)*EL^2*(-dZH1 + dMWsq1/MW^2)}}, 
 C[V[10], -V[30], V[30]] == {{(-I)*EL, (-I)*EL*(dZH1 - dMWsq1/MW^2)}, {0, 0}, 
   {0, 0}, {0, 0}}, C[V[10], -V[30], V[3]] == 
  {{(-I)*EL, (-I)*EL*(dZH1 - dMWsq1/MW^2)}, {0, 0}, {0, 0}, {0, 0}}, 
 C[V[10], -V[3], V[30]] == {{(-I)*EL, (-I)*EL*(dZH1 - dMWsq1/MW^2)}, {0, 0}, 
   {0, 0}, {0, 0}}, C[V[1], -V[30], V[30]] == 
  {{(-I)*EL, (-I)*EL*(dZH1 - dMWsq1/MW^2)}, {0, 0}, {0, 0}, {0, 0}}, 
 C[V[20], -V[30], V[30]] == 
  {{(I*CW*EL)/SW, (I*CW*EL*(-dMWsq1 + dZH1*MW^2))/(MW^2*SW)}, {0, 0}, {0, 0}, 
   {0, 0}}, C[V[20], -V[30], V[3]] == 
  {{(I*CW*EL)/SW, (I*CW*EL*(-dMWsq1 + dZH1*MW^2))/(MW^2*SW)}, {0, 0}, {0, 0}, 
   {0, 0}}, C[V[20], -V[3], V[30]] == 
  {{(I*CW*EL)/SW, (I*CW*EL*(-dMWsq1 + dZH1*MW^2))/(MW^2*SW)}, {0, 0}, {0, 0}, 
   {0, 0}}, C[V[2], -V[30], V[30]] == 
  {{(I*CW*EL)/SW, (I*CW*EL*(-dMWsq1 + dZH1*MW^2))/(MW^2*SW)}, {0, 0}, {0, 0}, 
   {0, 0}}, C[S[10], S[10], S[10], S[10]] == 
  {{(((-3*I)/4)*EL^2*MH^2)/(MW^2*SW^2), 
    (((-3*I)/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[10], S[10], S[10], S[1]] == {{(((-3*I)/4)*EL^2*MH^2)/(MW^2*SW^2), 
    (((-3*I)/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[10], S[1], S[1], S[1]] == {{(((-3*I)/4)*EL^2*MH^2)/(MW^2*SW^2), 
    (((-3*I)/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[10], S[10], S[20], S[20]] == {{((-I/4)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[10], S[10], S[20], S[2]] == {{((-I/4)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[10], S[1], S[20], S[20]] == {{((-I/4)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[10], S[1], S[2], S[2]] == {{((-I/4)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[1], S[1], S[20], S[2]] == {{((-I/4)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[10], S[10], S[30], -S[30]] == {{((-I/4)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[10], S[10], S[30], -S[3]] == {{((-I/4)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[10], S[10], S[3], -S[30]] == {{((-I/4)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[10], S[1], S[30], -S[30]] == {{((-I/4)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[10], S[1], S[3], -S[3]] == {{((-I/4)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[1], S[1], S[30], -S[3]] == {{((-I/4)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[1], S[1], S[3], -S[30]] == {{((-I/4)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[20], S[20], S[20], S[20]] == {{(((-3*I)/4)*EL^2*MH^2)/(MW^2*SW^2), 
    (((-3*I)/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[20], S[20], S[20], S[2]] == {{(((-3*I)/4)*EL^2*MH^2)/(MW^2*SW^2), 
    (((-3*I)/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[20], S[2], S[2], S[2]] == {{(((-3*I)/4)*EL^2*MH^2)/(MW^2*SW^2), 
    (((-3*I)/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[20], S[20], S[30], -S[30]] == {{((-I/4)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[20], S[20], S[30], -S[3]] == {{((-I/4)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[20], S[20], S[3], -S[30]] == {{((-I/4)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[20], S[2], S[30], -S[30]] == {{((-I/4)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[20], S[2], S[3], -S[3]] == {{((-I/4)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[2], S[2], S[30], -S[3]] == {{((-I/4)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[2], S[2], S[3], -S[30]] == {{((-I/4)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/8)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[30], S[30], -S[30], -S[30]] == 
  {{((-I/2)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/4)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[30], S[30], -S[30], -S[3]] == {{((-I/2)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/4)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[30], S[3], -S[30], -S[30]] == {{((-I/2)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/4)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[30], S[3], -S[3], -S[3]] == {{((-I/2)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/4)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[3], S[3], -S[30], -S[3]] == {{((-I/2)*EL^2*MH^2)/(MW^2*SW^2), 
    ((-I/4)*EL^2*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^3*SW^3)}}, 
 C[S[10], S[10], S[10]] == {{(((-3*I)/2)*EL*MH^2)/(MW*SW), 
    (((-3*I)/4)*EL*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^2*SW^2)}}, 
 C[S[10], S[10], S[1]] == {{(((-3*I)/2)*EL*MH^2)/(MW*SW), 
    (((-3*I)/4)*EL*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^2*SW^2)}}, 
 C[S[10], S[20], S[20]] == {{((-I/2)*EL*MH^2)/(MW*SW), 
    ((-I/4)*EL*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^2*SW^2)}}, 
 C[S[10], S[20], S[2]] == {{((-I/2)*EL*MH^2)/(MW*SW), 
    ((-I/4)*EL*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^2*SW^2)}}, 
 C[S[1], S[20], S[20]] == {{((-I/2)*EL*MH^2)/(MW*SW), 
    ((-I/4)*EL*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^2*SW^2)}}, 
 C[S[30], S[10], -S[30]] == {{((-I/2)*EL*MH^2)/(MW*SW), 
    ((-I/4)*EL*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^2*SW^2)}}, 
 C[S[30], S[10], -S[3]] == {{((-I/2)*EL*MH^2)/(MW*SW), 
    ((-I/4)*EL*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^2*SW^2)}}, 
 C[S[30], S[1], -S[30]] == {{((-I/2)*EL*MH^2)/(MW*SW), 
    ((-I/4)*EL*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^2*SW^2)}}, 
 C[S[3], S[10], -S[30]] == {{((-I/2)*EL*MH^2)/(MW*SW), 
    ((-I/4)*EL*(dTH1*EL + 2*(dMHsq1 + dZH1*MH^2)*MW*SW))/(MW^2*SW^2)}}, 
 C[S[10], S[10], V[30], -V[30]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[10], S[10], V[30], -V[3]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[10], S[10], V[3], -V[30]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[10], S[1], V[30], -V[30]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[10], S[1], V[3], -V[3]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[1], S[1], V[30], -V[3]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[1], S[1], V[3], -V[30]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[1], S[1], V[30], -V[30]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[10], S[10], V[3], -V[3]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[20], S[20], V[30], -V[30]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[20], S[20], V[30], -V[3]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[20], S[20], V[3], -V[30]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[20], S[2], V[30], -V[30]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[20], S[2], V[3], -V[3]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[2], S[2], V[30], -V[3]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[2], S[2], V[3], -V[30]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[2], S[2], V[30], -V[30]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[20], S[20], V[3], -V[3]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[30], -S[30], V[30], -V[30]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[30], -S[30], V[30], -V[3]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[30], -S[30], V[3], -V[30]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[30], -S[3], V[30], -V[30]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[3], -S[30], V[30], -V[30]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[30], -S[3], V[3], -V[3]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[3], -S[30], V[3], -V[3]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[3], -S[3], V[30], -V[3]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[3], -S[3], V[3], -V[30]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[3], -S[3], V[30], -V[30]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[30], -S[30], V[3], -V[3]] == 
  {{((I/2)*EL^2)/SW^2, ((I/2)*dZH1*EL^2)/SW^2}}, 
 C[S[30], -S[30], V[20], V[20]] == 
  {{((I/2)*EL^2*(-CW^2 + SW^2)^2)/(CW^2*SW^2), 
    ((I/2)*dZH1*EL^2*(-CW^2 + SW^2)^2)/(CW^2*SW^2)}}, 
 C[S[30], -S[30], V[20], V[2]] == 
  {{((I/2)*EL^2*(-CW^2 + SW^2)^2)/(CW^2*SW^2), 
    ((I/2)*dZH1*EL^2*(-CW^2 + SW^2)^2)/(CW^2*SW^2)}}, 
 C[S[30], -S[3], V[20], V[20]] == 
  {{((I/2)*EL^2*(-CW^2 + SW^2)^2)/(CW^2*SW^2), 
    ((I/2)*dZH1*EL^2*(-CW^2 + SW^2)^2)/(CW^2*SW^2)}}, 
 C[S[3], -S[30], V[20], V[20]] == 
  {{((I/2)*EL^2*(-CW^2 + SW^2)^2)/(CW^2*SW^2), 
    ((I/2)*dZH1*EL^2*(-CW^2 + SW^2)^2)/(CW^2*SW^2)}}, 
 C[S[30], -S[3], V[2], V[2]] == {{((I/2)*EL^2*(-CW^2 + SW^2)^2)/(CW^2*SW^2), 
    ((I/2)*dZH1*EL^2*(-CW^2 + SW^2)^2)/(CW^2*SW^2)}}, 
 C[S[3], -S[30], V[2], V[2]] == {{((I/2)*EL^2*(-CW^2 + SW^2)^2)/(CW^2*SW^2), 
    ((I/2)*dZH1*EL^2*(-CW^2 + SW^2)^2)/(CW^2*SW^2)}}, 
 C[S[3], -S[3], V[20], V[2]] == {{((I/2)*EL^2*(-CW^2 + SW^2)^2)/(CW^2*SW^2), 
    ((I/2)*dZH1*EL^2*(-CW^2 + SW^2)^2)/(CW^2*SW^2)}}, 
 C[S[3], -S[3], V[20], V[20]] == {{((I/2)*EL^2*(-CW^2 + SW^2)^2)/(CW^2*SW^2), 
    ((I/2)*dZH1*EL^2*(-CW^2 + SW^2)^2)/(CW^2*SW^2)}}, 
 C[S[30], -S[30], V[2], V[2]] == {{((I/2)*EL^2*(-CW^2 + SW^2)^2)/(CW^2*SW^2), 
    ((I/2)*dZH1*EL^2*(-CW^2 + SW^2)^2)/(CW^2*SW^2)}}, 
 C[S[30], -S[30], V[10], V[20]] == {{(I*EL^2*(-CW^2 + SW^2))/(CW*SW), 
    ((-I)*dZH1*EL^2*(CW^2 - SW^2))/(CW*SW)}}, 
 C[S[30], -S[30], V[10], V[2]] == {{(I*EL^2*(-CW^2 + SW^2))/(CW*SW), 
    ((-I)*dZH1*EL^2*(CW^2 - SW^2))/(CW*SW)}}, 
 C[S[30], -S[30], V[1], V[20]] == {{(I*EL^2*(-CW^2 + SW^2))/(CW*SW), 
    ((-I)*dZH1*EL^2*(CW^2 - SW^2))/(CW*SW)}}, 
 C[S[30], -S[3], V[10], V[20]] == {{(I*EL^2*(-CW^2 + SW^2))/(CW*SW), 
    ((-I)*dZH1*EL^2*(CW^2 - SW^2))/(CW*SW)}}, 
 C[S[3], -S[30], V[10], V[20]] == {{(I*EL^2*(-CW^2 + SW^2))/(CW*SW), 
    ((-I)*dZH1*EL^2*(CW^2 - SW^2))/(CW*SW)}}, C[S[30], -S[3], V[1], V[2]] == 
  {{(I*EL^2*(-CW^2 + SW^2))/(CW*SW), ((-I)*dZH1*EL^2*(CW^2 - SW^2))/
     (CW*SW)}}, C[S[3], -S[30], V[1], V[2]] == 
  {{(I*EL^2*(-CW^2 + SW^2))/(CW*SW), ((-I)*dZH1*EL^2*(CW^2 - SW^2))/
     (CW*SW)}}, C[S[3], -S[3], V[10], V[2]] == 
  {{(I*EL^2*(-CW^2 + SW^2))/(CW*SW), ((-I)*dZH1*EL^2*(CW^2 - SW^2))/
     (CW*SW)}}, C[S[3], -S[3], V[1], V[20]] == 
  {{(I*EL^2*(-CW^2 + SW^2))/(CW*SW), ((-I)*dZH1*EL^2*(CW^2 - SW^2))/
     (CW*SW)}}, C[S[3], -S[3], V[10], V[20]] == 
  {{(I*EL^2*(-CW^2 + SW^2))/(CW*SW), ((-I)*dZH1*EL^2*(CW^2 - SW^2))/
     (CW*SW)}}, C[S[30], -S[30], V[1], V[2]] == 
  {{(I*EL^2*(-CW^2 + SW^2))/(CW*SW), ((-I)*dZH1*EL^2*(CW^2 - SW^2))/
     (CW*SW)}}, C[S[30], -S[30], V[10], V[10]] == 
  {{(2*I)*EL^2, (2*I)*dZH1*EL^2}}, C[S[30], -S[30], V[10], V[1]] == 
  {{(2*I)*EL^2, (2*I)*dZH1*EL^2}}, C[S[30], -S[3], V[10], V[10]] == 
  {{(2*I)*EL^2, (2*I)*dZH1*EL^2}}, C[S[3], -S[30], V[10], V[10]] == 
  {{(2*I)*EL^2, (2*I)*dZH1*EL^2}}, C[S[30], -S[3], V[1], V[1]] == 
  {{(2*I)*EL^2, (2*I)*dZH1*EL^2}}, C[S[3], -S[30], V[1], V[1]] == 
  {{(2*I)*EL^2, (2*I)*dZH1*EL^2}}, C[S[3], -S[3], V[10], V[1]] == 
  {{(2*I)*EL^2, (2*I)*dZH1*EL^2}}, C[S[3], -S[3], V[10], V[10]] == 
  {{(2*I)*EL^2, (2*I)*dZH1*EL^2}}, C[S[30], -S[30], V[1], V[1]] == 
  {{(2*I)*EL^2, (2*I)*dZH1*EL^2}}, C[S[10], S[10], V[20], V[20]] == 
  {{((I/2)*EL^2)/(CW^2*SW^2), ((I/2)*dZH1*EL^2)/(CW^2*SW^2)}}, 
 C[S[10], S[10], V[20], V[2]] == 
  {{((I/2)*EL^2)/(CW^2*SW^2), ((I/2)*dZH1*EL^2)/(CW^2*SW^2)}}, 
 C[S[10], S[1], V[20], V[20]] == 
  {{((I/2)*EL^2)/(CW^2*SW^2), ((I/2)*dZH1*EL^2)/(CW^2*SW^2)}}, 
 C[S[10], S[1], V[2], V[2]] == {{((I/2)*EL^2)/(CW^2*SW^2), 
    ((I/2)*dZH1*EL^2)/(CW^2*SW^2)}}, C[S[1], S[1], V[20], V[2]] == 
  {{((I/2)*EL^2)/(CW^2*SW^2), ((I/2)*dZH1*EL^2)/(CW^2*SW^2)}}, 
 C[S[1], S[1], V[20], V[20]] == {{((I/2)*EL^2)/(CW^2*SW^2), 
    ((I/2)*dZH1*EL^2)/(CW^2*SW^2)}}, C[S[10], S[10], V[2], V[2]] == 
  {{((I/2)*EL^2)/(CW^2*SW^2), ((I/2)*dZH1*EL^2)/(CW^2*SW^2)}}, 
 C[S[20], S[20], V[20], V[20]] == 
  {{((I/2)*EL^2)/(CW^2*SW^2), ((I/2)*dZH1*EL^2)/(CW^2*SW^2)}}, 
 C[S[20], S[20], V[20], V[2]] == 
  {{((I/2)*EL^2)/(CW^2*SW^2), ((I/2)*dZH1*EL^2)/(CW^2*SW^2)}}, 
 C[S[20], S[2], V[20], V[20]] == 
  {{((I/2)*EL^2)/(CW^2*SW^2), ((I/2)*dZH1*EL^2)/(CW^2*SW^2)}}, 
 C[S[20], S[2], V[2], V[2]] == {{((I/2)*EL^2)/(CW^2*SW^2), 
    ((I/2)*dZH1*EL^2)/(CW^2*SW^2)}}, C[S[2], S[2], V[20], V[2]] == 
  {{((I/2)*EL^2)/(CW^2*SW^2), ((I/2)*dZH1*EL^2)/(CW^2*SW^2)}}, 
 C[S[2], S[2], V[20], V[20]] == {{((I/2)*EL^2)/(CW^2*SW^2), 
    ((I/2)*dZH1*EL^2)/(CW^2*SW^2)}}, C[S[20], S[20], V[2], V[2]] == 
  {{((I/2)*EL^2)/(CW^2*SW^2), ((I/2)*dZH1*EL^2)/(CW^2*SW^2)}}, 
 C[S[10], -S[30], V[30], V[20]] == 
  {{((-I/2)*EL^2)/CW, ((-I/2)*dZH1*EL^2)/CW}}, 
 C[S[10], -S[30], V[30], V[2]] == 
  {{((-I/2)*EL^2)/CW, ((-I/2)*dZH1*EL^2)/CW}}, 
 C[S[10], -S[30], V[3], V[20]] == 
  {{((-I/2)*EL^2)/CW, ((-I/2)*dZH1*EL^2)/CW}}, 
 C[S[10], -S[3], V[30], V[20]] == 
  {{((-I/2)*EL^2)/CW, ((-I/2)*dZH1*EL^2)/CW}}, 
 C[S[1], -S[30], V[30], V[20]] == 
  {{((-I/2)*EL^2)/CW, ((-I/2)*dZH1*EL^2)/CW}}, 
 C[S[10], -S[3], V[3], V[2]] == {{((-I/2)*EL^2)/CW, ((-I/2)*dZH1*EL^2)/CW}}, 
 C[S[1], -S[30], V[3], V[2]] == {{((-I/2)*EL^2)/CW, ((-I/2)*dZH1*EL^2)/CW}}, 
 C[S[1], -S[3], V[30], V[2]] == {{((-I/2)*EL^2)/CW, ((-I/2)*dZH1*EL^2)/CW}}, 
 C[S[1], -S[3], V[3], V[20]] == {{((-I/2)*EL^2)/CW, ((-I/2)*dZH1*EL^2)/CW}}, 
 C[S[1], -S[3], V[30], V[20]] == {{((-I/2)*EL^2)/CW, ((-I/2)*dZH1*EL^2)/CW}}, 
 C[S[10], -S[30], V[3], V[2]] == {{((-I/2)*EL^2)/CW, ((-I/2)*dZH1*EL^2)/CW}}, 
 C[S[10], S[30], -V[30], V[20]] == 
  {{((-I/2)*EL^2)/CW, ((-I/2)*dZH1*EL^2)/CW}}, 
 C[S[10], S[30], -V[30], V[2]] == 
  {{((-I/2)*EL^2)/CW, ((-I/2)*dZH1*EL^2)/CW}}, 
 C[S[10], S[30], -V[3], V[20]] == 
  {{((-I/2)*EL^2)/CW, ((-I/2)*dZH1*EL^2)/CW}}, 
 C[S[10], S[3], -V[30], V[20]] == 
  {{((-I/2)*EL^2)/CW, ((-I/2)*dZH1*EL^2)/CW}}, 
 C[S[1], S[30], -V[30], V[20]] == 
  {{((-I/2)*EL^2)/CW, ((-I/2)*dZH1*EL^2)/CW}}, 
 C[S[10], S[3], -V[3], V[2]] == {{((-I/2)*EL^2)/CW, ((-I/2)*dZH1*EL^2)/CW}}, 
 C[S[1], S[30], -V[3], V[2]] == {{((-I/2)*EL^2)/CW, ((-I/2)*dZH1*EL^2)/CW}}, 
 C[S[1], S[3], -V[30], V[2]] == {{((-I/2)*EL^2)/CW, ((-I/2)*dZH1*EL^2)/CW}}, 
 C[S[1], S[3], -V[3], V[20]] == {{((-I/2)*EL^2)/CW, ((-I/2)*dZH1*EL^2)/CW}}, 
 C[S[1], S[3], -V[30], V[20]] == {{((-I/2)*EL^2)/CW, ((-I/2)*dZH1*EL^2)/CW}}, 
 C[S[10], S[30], -V[3], V[2]] == {{((-I/2)*EL^2)/CW, ((-I/2)*dZH1*EL^2)/CW}}, 
 C[S[10], S[30], -V[30], V[10]] == 
  {{((-I/2)*EL^2)/SW, ((-I/2)*dZH1*EL^2)/SW}}, 
 C[S[10], S[30], -V[30], V[1]] == 
  {{((-I/2)*EL^2)/SW, ((-I/2)*dZH1*EL^2)/SW}}, 
 C[S[10], S[30], -V[3], V[10]] == 
  {{((-I/2)*EL^2)/SW, ((-I/2)*dZH1*EL^2)/SW}}, 
 C[S[10], S[3], -V[30], V[10]] == 
  {{((-I/2)*EL^2)/SW, ((-I/2)*dZH1*EL^2)/SW}}, 
 C[S[1], S[30], -V[30], V[10]] == 
  {{((-I/2)*EL^2)/SW, ((-I/2)*dZH1*EL^2)/SW}}, 
 C[S[10], S[3], -V[3], V[1]] == {{((-I/2)*EL^2)/SW, ((-I/2)*dZH1*EL^2)/SW}}, 
 C[S[1], S[30], -V[3], V[1]] == {{((-I/2)*EL^2)/SW, ((-I/2)*dZH1*EL^2)/SW}}, 
 C[S[1], S[3], -V[30], V[1]] == {{((-I/2)*EL^2)/SW, ((-I/2)*dZH1*EL^2)/SW}}, 
 C[S[1], S[3], -V[3], V[10]] == {{((-I/2)*EL^2)/SW, ((-I/2)*dZH1*EL^2)/SW}}, 
 C[S[1], S[3], -V[30], V[10]] == {{((-I/2)*EL^2)/SW, ((-I/2)*dZH1*EL^2)/SW}}, 
 C[S[10], S[30], -V[3], V[1]] == {{((-I/2)*EL^2)/SW, ((-I/2)*dZH1*EL^2)/SW}}, 
 C[S[10], -S[30], V[30], V[10]] == 
  {{((-I/2)*EL^2)/SW, ((-I/2)*dZH1*EL^2)/SW}}, 
 C[S[10], -S[30], V[30], V[1]] == 
  {{((-I/2)*EL^2)/SW, ((-I/2)*dZH1*EL^2)/SW}}, 
 C[S[10], -S[30], V[3], V[10]] == 
  {{((-I/2)*EL^2)/SW, ((-I/2)*dZH1*EL^2)/SW}}, 
 C[S[10], -S[3], V[30], V[10]] == 
  {{((-I/2)*EL^2)/SW, ((-I/2)*dZH1*EL^2)/SW}}, 
 C[S[1], -S[30], V[30], V[10]] == 
  {{((-I/2)*EL^2)/SW, ((-I/2)*dZH1*EL^2)/SW}}, 
 C[S[10], -S[3], V[3], V[1]] == {{((-I/2)*EL^2)/SW, ((-I/2)*dZH1*EL^2)/SW}}, 
 C[S[1], -S[30], V[3], V[1]] == {{((-I/2)*EL^2)/SW, ((-I/2)*dZH1*EL^2)/SW}}, 
 C[S[1], -S[3], V[30], V[1]] == {{((-I/2)*EL^2)/SW, ((-I/2)*dZH1*EL^2)/SW}}, 
 C[S[1], -S[3], V[3], V[10]] == {{((-I/2)*EL^2)/SW, ((-I/2)*dZH1*EL^2)/SW}}, 
 C[S[1], -S[3], V[30], V[10]] == {{((-I/2)*EL^2)/SW, ((-I/2)*dZH1*EL^2)/SW}}, 
 C[S[10], -S[30], V[3], V[1]] == {{((-I/2)*EL^2)/SW, ((-I/2)*dZH1*EL^2)/SW}}, 
 C[S[30], S[20], V[20], -V[30]] == {{EL^2/(2*CW), (dZH1*EL^2)/(2*CW)}}, 
 C[S[30], S[20], V[20], -V[3]] == {{EL^2/(2*CW), (dZH1*EL^2)/(2*CW)}}, 
 C[S[30], S[20], V[2], -V[30]] == {{EL^2/(2*CW), (dZH1*EL^2)/(2*CW)}}, 
 C[S[30], S[2], V[20], -V[30]] == {{EL^2/(2*CW), (dZH1*EL^2)/(2*CW)}}, 
 C[S[3], S[20], V[20], -V[30]] == {{EL^2/(2*CW), (dZH1*EL^2)/(2*CW)}}, 
 C[S[30], S[2], V[2], -V[3]] == {{EL^2/(2*CW), (dZH1*EL^2)/(2*CW)}}, 
 C[S[3], S[20], V[2], -V[3]] == {{EL^2/(2*CW), (dZH1*EL^2)/(2*CW)}}, 
 C[S[3], S[2], V[20], -V[3]] == {{EL^2/(2*CW), (dZH1*EL^2)/(2*CW)}}, 
 C[S[3], S[2], V[2], -V[30]] == {{EL^2/(2*CW), (dZH1*EL^2)/(2*CW)}}, 
 C[S[3], S[2], V[20], -V[30]] == {{EL^2/(2*CW), (dZH1*EL^2)/(2*CW)}}, 
 C[S[30], S[20], V[2], -V[3]] == {{EL^2/(2*CW), (dZH1*EL^2)/(2*CW)}}, 
 C[-S[30], S[20], V[20], V[30]] == {{-EL^2/(2*CW), -(dZH1*EL^2)/(2*CW)}}, 
 C[-S[30], S[20], V[20], V[3]] == {{-EL^2/(2*CW), -(dZH1*EL^2)/(2*CW)}}, 
 C[-S[30], S[20], V[2], V[30]] == {{-EL^2/(2*CW), -(dZH1*EL^2)/(2*CW)}}, 
 C[-S[30], S[2], V[20], V[30]] == {{-EL^2/(2*CW), -(dZH1*EL^2)/(2*CW)}}, 
 C[-S[3], S[20], V[20], V[30]] == {{-EL^2/(2*CW), -(dZH1*EL^2)/(2*CW)}}, 
 C[-S[30], S[2], V[2], V[3]] == {{-EL^2/(2*CW), -(dZH1*EL^2)/(2*CW)}}, 
 C[-S[3], S[20], V[2], V[3]] == {{-EL^2/(2*CW), -(dZH1*EL^2)/(2*CW)}}, 
 C[-S[3], S[2], V[20], V[3]] == {{-EL^2/(2*CW), -(dZH1*EL^2)/(2*CW)}}, 
 C[-S[3], S[2], V[2], V[30]] == {{-EL^2/(2*CW), -(dZH1*EL^2)/(2*CW)}}, 
 C[-S[3], S[2], V[20], V[30]] == {{-EL^2/(2*CW), -(dZH1*EL^2)/(2*CW)}}, 
 C[-S[30], S[20], V[2], V[3]] == {{-EL^2/(2*CW), -(dZH1*EL^2)/(2*CW)}}, 
 C[S[30], S[20], V[10], -V[30]] == {{EL^2/(2*SW), (dZH1*EL^2)/(2*SW)}}, 
 C[S[30], S[20], V[10], -V[3]] == {{EL^2/(2*SW), (dZH1*EL^2)/(2*SW)}}, 
 C[S[30], S[20], V[1], -V[30]] == {{EL^2/(2*SW), (dZH1*EL^2)/(2*SW)}}, 
 C[S[30], S[2], V[10], -V[30]] == {{EL^2/(2*SW), (dZH1*EL^2)/(2*SW)}}, 
 C[S[3], S[20], V[10], -V[30]] == {{EL^2/(2*SW), (dZH1*EL^2)/(2*SW)}}, 
 C[S[30], S[2], V[1], -V[3]] == {{EL^2/(2*SW), (dZH1*EL^2)/(2*SW)}}, 
 C[S[3], S[20], V[1], -V[3]] == {{EL^2/(2*SW), (dZH1*EL^2)/(2*SW)}}, 
 C[S[3], S[2], V[10], -V[3]] == {{EL^2/(2*SW), (dZH1*EL^2)/(2*SW)}}, 
 C[S[3], S[2], V[1], -V[30]] == {{EL^2/(2*SW), (dZH1*EL^2)/(2*SW)}}, 
 C[S[3], S[2], V[10], -V[30]] == {{EL^2/(2*SW), (dZH1*EL^2)/(2*SW)}}, 
 C[S[30], S[20], V[1], -V[3]] == {{EL^2/(2*SW), (dZH1*EL^2)/(2*SW)}}, 
 C[-S[30], S[20], V[10], V[30]] == {{-EL^2/(2*SW), -(dZH1*EL^2)/(2*SW)}}, 
 C[-S[30], S[20], V[10], V[3]] == {{-EL^2/(2*SW), -(dZH1*EL^2)/(2*SW)}}, 
 C[-S[30], S[20], V[1], V[30]] == {{-EL^2/(2*SW), -(dZH1*EL^2)/(2*SW)}}, 
 C[-S[30], S[2], V[10], V[30]] == {{-EL^2/(2*SW), -(dZH1*EL^2)/(2*SW)}}, 
 C[-S[3], S[20], V[10], V[30]] == {{-EL^2/(2*SW), -(dZH1*EL^2)/(2*SW)}}, 
 C[-S[30], S[2], V[1], V[3]] == {{-EL^2/(2*SW), -(dZH1*EL^2)/(2*SW)}}, 
 C[-S[3], S[20], V[1], V[3]] == {{-EL^2/(2*SW), -(dZH1*EL^2)/(2*SW)}}, 
 C[-S[3], S[2], V[10], V[3]] == {{-EL^2/(2*SW), -(dZH1*EL^2)/(2*SW)}}, 
 C[-S[3], S[2], V[1], V[30]] == {{-EL^2/(2*SW), -(dZH1*EL^2)/(2*SW)}}, 
 C[-S[3], S[2], V[10], V[30]] == {{-EL^2/(2*SW), -(dZH1*EL^2)/(2*SW)}}, 
 C[-S[30], S[20], V[1], V[3]] == {{-EL^2/(2*SW), -(dZH1*EL^2)/(2*SW)}}, 
 C[S[20], S[10], V[20]] == {{EL/(2*CW*SW), (dZH1*EL)/(2*CW*SW)}, 
   {-EL/(2*CW*SW), -(dZH1*EL)/(2*CW*SW)}}, C[S[20], S[10], V[2]] == 
  {{EL/(2*CW*SW), (dZH1*EL)/(2*CW*SW)}, {-EL/(2*CW*SW), 
    -(dZH1*EL)/(2*CW*SW)}}, C[S[20], S[1], V[20]] == 
  {{EL/(2*CW*SW), (dZH1*EL)/(2*CW*SW)}, {-EL/(2*CW*SW), 
    -(dZH1*EL)/(2*CW*SW)}}, C[S[2], S[10], V[20]] == 
  {{EL/(2*CW*SW), (dZH1*EL)/(2*CW*SW)}, {-EL/(2*CW*SW), 
    -(dZH1*EL)/(2*CW*SW)}}, C[S[2], S[1], V[20]] == 
  {{EL/(2*CW*SW), (dZH1*EL)/(2*CW*SW)}, {-EL/(2*CW*SW), 
    -(dZH1*EL)/(2*CW*SW)}}, C[-S[30], S[30], V[10]] == 
  {{(-I)*EL, (-I)*dZH1*EL}, {I*EL, I*dZH1*EL}}, 
 C[-S[30], S[30], V[1]] == {{(-I)*EL, (-I)*dZH1*EL}, {I*EL, I*dZH1*EL}}, 
 C[-S[30], S[3], V[10]] == {{(-I)*EL, (-I)*dZH1*EL}, {I*EL, I*dZH1*EL}}, 
 C[-S[3], S[30], V[10]] == {{(-I)*EL, (-I)*dZH1*EL}, {I*EL, I*dZH1*EL}}, 
 C[-S[3], S[3], V[10]] == {{(-I)*EL, (-I)*dZH1*EL}, {I*EL, I*dZH1*EL}}, 
 C[-S[30], S[30], V[20]] == {{((I/2)*EL*(CW^2 - SW^2))/(CW*SW), 
    ((I/2)*dZH1*EL*(CW^2 - SW^2))/(CW*SW)}, 
   {((-I/2)*EL*(CW^2 - SW^2))/(CW*SW), ((-I/2)*dZH1*EL*(CW^2 - SW^2))/
     (CW*SW)}}, C[-S[30], S[30], V[2]] == 
  {{((I/2)*EL*(CW^2 - SW^2))/(CW*SW), ((I/2)*dZH1*EL*(CW^2 - SW^2))/(CW*SW)}, 
   {((-I/2)*EL*(CW^2 - SW^2))/(CW*SW), ((-I/2)*dZH1*EL*(CW^2 - SW^2))/
     (CW*SW)}}, C[-S[30], S[3], V[20]] == 
  {{((I/2)*EL*(CW^2 - SW^2))/(CW*SW), ((I/2)*dZH1*EL*(CW^2 - SW^2))/(CW*SW)}, 
   {((-I/2)*EL*(CW^2 - SW^2))/(CW*SW), ((-I/2)*dZH1*EL*(CW^2 - SW^2))/
     (CW*SW)}}, C[-S[3], S[30], V[20]] == 
  {{((I/2)*EL*(CW^2 - SW^2))/(CW*SW), ((I/2)*dZH1*EL*(CW^2 - SW^2))/(CW*SW)}, 
   {((-I/2)*EL*(CW^2 - SW^2))/(CW*SW), ((-I/2)*dZH1*EL*(CW^2 - SW^2))/
     (CW*SW)}}, C[-S[3], S[3], V[20]] == 
  {{((I/2)*EL*(CW^2 - SW^2))/(CW*SW), ((I/2)*dZH1*EL*(CW^2 - SW^2))/(CW*SW)}, 
   {((-I/2)*EL*(CW^2 - SW^2))/(CW*SW), ((-I/2)*dZH1*EL*(CW^2 - SW^2))/
     (CW*SW)}}, C[S[30], S[10], -V[30]] == 
  {{((-I/2)*EL)/SW, ((-I/2)*dZH1*EL)/SW}, {((I/2)*EL)/SW, 
    ((I/2)*dZH1*EL)/SW}}, C[S[30], S[10], -V[3]] == 
  {{((-I/2)*EL)/SW, ((-I/2)*dZH1*EL)/SW}, {((I/2)*EL)/SW, 
    ((I/2)*dZH1*EL)/SW}}, C[S[30], S[1], -V[30]] == 
  {{((-I/2)*EL)/SW, ((-I/2)*dZH1*EL)/SW}, {((I/2)*EL)/SW, 
    ((I/2)*dZH1*EL)/SW}}, C[S[3], S[10], -V[30]] == 
  {{((-I/2)*EL)/SW, ((-I/2)*dZH1*EL)/SW}, {((I/2)*EL)/SW, 
    ((I/2)*dZH1*EL)/SW}}, C[S[3], S[1], -V[30]] == 
  {{((-I/2)*EL)/SW, ((-I/2)*dZH1*EL)/SW}, {((I/2)*EL)/SW, 
    ((I/2)*dZH1*EL)/SW}}, C[-S[30], S[10], V[30]] == 
  {{((I/2)*EL)/SW, ((I/2)*dZH1*EL)/SW}, {((-I/2)*EL)/SW, 
    ((-I/2)*dZH1*EL)/SW}}, C[-S[30], S[10], V[3]] == 
  {{((I/2)*EL)/SW, ((I/2)*dZH1*EL)/SW}, {((-I/2)*EL)/SW, 
    ((-I/2)*dZH1*EL)/SW}}, C[-S[30], S[1], V[30]] == 
  {{((I/2)*EL)/SW, ((I/2)*dZH1*EL)/SW}, {((-I/2)*EL)/SW, 
    ((-I/2)*dZH1*EL)/SW}}, C[-S[3], S[10], V[30]] == 
  {{((I/2)*EL)/SW, ((I/2)*dZH1*EL)/SW}, {((-I/2)*EL)/SW, 
    ((-I/2)*dZH1*EL)/SW}}, C[-S[3], S[1], V[30]] == 
  {{((I/2)*EL)/SW, ((I/2)*dZH1*EL)/SW}, {((-I/2)*EL)/SW, 
    ((-I/2)*dZH1*EL)/SW}}, C[S[30], S[20], -V[30]] == 
  {{EL/(2*SW), (dZH1*EL)/(2*SW)}, {-EL/(2*SW), -(dZH1*EL)/(2*SW)}}, 
 C[S[30], S[20], -V[3]] == {{EL/(2*SW), (dZH1*EL)/(2*SW)}, 
   {-EL/(2*SW), -(dZH1*EL)/(2*SW)}}, C[S[30], S[2], -V[30]] == 
  {{EL/(2*SW), (dZH1*EL)/(2*SW)}, {-EL/(2*SW), -(dZH1*EL)/(2*SW)}}, 
 C[S[3], S[20], -V[30]] == {{EL/(2*SW), (dZH1*EL)/(2*SW)}, 
   {-EL/(2*SW), -(dZH1*EL)/(2*SW)}}, C[S[3], S[2], -V[30]] == 
  {{EL/(2*SW), (dZH1*EL)/(2*SW)}, {-EL/(2*SW), -(dZH1*EL)/(2*SW)}}, 
 C[-S[30], S[20], V[30]] == {{EL/(2*SW), (dZH1*EL)/(2*SW)}, 
   {-EL/(2*SW), -(dZH1*EL)/(2*SW)}}, C[-S[30], S[20], V[3]] == 
  {{EL/(2*SW), (dZH1*EL)/(2*SW)}, {-EL/(2*SW), -(dZH1*EL)/(2*SW)}}, 
 C[-S[30], S[2], V[30]] == {{EL/(2*SW), (dZH1*EL)/(2*SW)}, 
   {-EL/(2*SW), -(dZH1*EL)/(2*SW)}}, C[-S[3], S[20], V[30]] == 
  {{EL/(2*SW), (dZH1*EL)/(2*SW)}, {-EL/(2*SW), -(dZH1*EL)/(2*SW)}}, 
 C[-S[3], S[2], V[30]] == {{EL/(2*SW), (dZH1*EL)/(2*SW)}, 
   {-EL/(2*SW), -(dZH1*EL)/(2*SW)}}, C[S[10], -V[30], V[30]] == 
  {{(I*EL*MW)/SW, (I*dZH1*EL*MW)/SW}}, C[S[10], -V[30], V[3]] == 
  {{(I*EL*MW)/SW, (I*dZH1*EL*MW)/SW}}, C[S[10], -V[3], V[30]] == 
  {{(I*EL*MW)/SW, (I*dZH1*EL*MW)/SW}}, C[S[1], -V[30], V[30]] == 
  {{(I*EL*MW)/SW, (I*dZH1*EL*MW)/SW}}, C[S[10], -V[3], V[3]] == 
  {{(I*EL*MW)/SW, (I*dZH1*EL*MW)/SW}}, C[S[10], V[20], V[20]] == 
  {{(I*EL*MW)/(CW^2*SW), (I*dZH1*EL*MW)/(CW^2*SW)}}, 
 C[S[10], V[20], V[2]] == {{(I*EL*MW)/(CW^2*SW), (I*dZH1*EL*MW)/(CW^2*SW)}}, 
 C[S[1], V[20], V[20]] == {{(I*EL*MW)/(CW^2*SW), (I*dZH1*EL*MW)/(CW^2*SW)}}, 
 C[S[10], V[2], V[2]] == {{(I*EL*MW)/(CW^2*SW), (I*dZH1*EL*MW)/(CW^2*SW)}}, 
 C[-S[30], V[30], V[20]] == {{((-I)*EL*MW*SW)/CW, ((-I)*dZH1*EL*MW*SW)/CW}}, 
 C[-S[30], V[30], V[2]] == {{((-I)*EL*MW*SW)/CW, ((-I)*dZH1*EL*MW*SW)/CW}}, 
 C[-S[30], V[3], V[20]] == {{((-I)*EL*MW*SW)/CW, ((-I)*dZH1*EL*MW*SW)/CW}}, 
 C[-S[3], V[30], V[20]] == {{((-I)*EL*MW*SW)/CW, ((-I)*dZH1*EL*MW*SW)/CW}}, 
 C[-S[30], V[3], V[2]] == {{((-I)*EL*MW*SW)/CW, ((-I)*dZH1*EL*MW*SW)/CW}}, 
 C[S[30], -V[30], V[20]] == {{((-I)*EL*MW*SW)/CW, ((-I)*dZH1*EL*MW*SW)/CW}}, 
 C[S[30], -V[30], V[2]] == {{((-I)*EL*MW*SW)/CW, ((-I)*dZH1*EL*MW*SW)/CW}}, 
 C[S[30], -V[3], V[20]] == {{((-I)*EL*MW*SW)/CW, ((-I)*dZH1*EL*MW*SW)/CW}}, 
 C[S[3], -V[30], V[20]] == {{((-I)*EL*MW*SW)/CW, ((-I)*dZH1*EL*MW*SW)/CW}}, 
 C[S[30], -V[3], V[2]] == {{((-I)*EL*MW*SW)/CW, ((-I)*dZH1*EL*MW*SW)/CW}}, 
 C[-S[30], V[30], V[10]] == {{(-I)*EL*MW, (-I)*dZH1*EL*MW}}, 
 C[-S[30], V[30], V[1]] == {{(-I)*EL*MW, (-I)*dZH1*EL*MW}}, 
 C[-S[30], V[3], V[10]] == {{(-I)*EL*MW, (-I)*dZH1*EL*MW}}, 
 C[-S[3], V[30], V[10]] == {{(-I)*EL*MW, (-I)*dZH1*EL*MW}}, 
 C[-S[30], V[3], V[1]] == {{(-I)*EL*MW, (-I)*dZH1*EL*MW}}, 
 C[S[30], -V[30], V[10]] == {{(-I)*EL*MW, (-I)*dZH1*EL*MW}}, 
 C[S[30], -V[30], V[1]] == {{(-I)*EL*MW, (-I)*dZH1*EL*MW}}, 
 C[S[30], -V[3], V[10]] == {{(-I)*EL*MW, (-I)*dZH1*EL*MW}}, 
 C[S[3], -V[30], V[10]] == {{(-I)*EL*MW, (-I)*dZH1*EL*MW}}, 
 C[S[30], -V[3], V[1]] == {{(-I)*EL*MW, (-I)*dZH1*EL*MW}}, 
 C[-F[2, {j1}], F[2, {j2}], V[10]] == 
  {{I*EL*IndexDelta[j1, j2], (I/2)*EL*(Conjugate[dZfL1[2, j1, j1]] + 
      dZfL1[2, j1, j1])*IndexDelta[j1, j2]}, {I*EL*IndexDelta[j1, j2], 
    (I/2)*EL*(Conjugate[dZfR1[2, j1, j1]] + dZfR1[2, j1, j1])*
     IndexDelta[j1, j2]}}, C[-F[3, {j1, o1}], F[3, {j2, o2}], V[10]] == 
  {{((-2*I)/3)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2], 
    (-I/3)*EL*(Conjugate[dZfL1[3, j2, j1]] + dZfL1[3, j1, j2])*
     IndexDelta[o1, o2]}, {((-2*I)/3)*EL*IndexDelta[j1, j2]*
     IndexDelta[o1, o2], (-I/3)*EL*(Conjugate[dZfR1[3, j2, j1]] + 
      dZfR1[3, j1, j2])*IndexDelta[o1, o2]}}, 
 C[-F[4, {j1, o1}], F[4, {j2, o2}], V[10]] == 
  {{(I/3)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2], 
    (I/6)*EL*(Conjugate[dZfL1[4, j2, j1]] + dZfL1[4, j1, j2])*
     IndexDelta[o1, o2]}, {(I/3)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2], 
    (I/6)*EL*(Conjugate[dZfR1[4, j2, j1]] + dZfR1[4, j1, j2])*
     IndexDelta[o1, o2]}}, C[-F[1, {j1}], F[1, {j2}], V[20]] == 
  {{((I/2)*EL*IndexDelta[j1, j2])/(CW*SW), 
    ((I/4)*EL*(Conjugate[dZfL1[1, j1, j1]] + dZfL1[1, j1, j1])*
      IndexDelta[j1, j2])/(CW*SW)}, {0, 0}}, 
 C[-F[2, {j1}], F[2, {j2}], V[20]] == 
  {{(I*(1/2 - CW^2)*EL*IndexDelta[j1, j2])/(CW*SW), 
    ((I/4)*EL*(CW^2*MW^2*(-CW^2 + SW^2)*Conjugate[dZfL1[2, j1, j1]] + 
       CW^2*MW^2*(-CW^2 + SW^2)*dZfL1[2, j1, j1])*IndexDelta[j1, j2])/
     (CW^3*MW^2*SW)}, {(I*EL*SW*IndexDelta[j1, j2])/CW, 
    ((I/2)*EL*SW*(Conjugate[dZfR1[2, j1, j1]] + dZfR1[2, j1, j1])*
      IndexDelta[j1, j2])/CW}}, C[-F[3, {j1, o1}], F[3, {j2, o2}], V[20]] == 
  {{((I/6)*(-1 + 4*CW^2)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2])/(CW*SW), 
    ((I/12)*EL*(CW^2*(-1 + 4*CW^2)*MW^2*Conjugate[dZfL1[3, j2, j1]] + 
       CW^2*(-1 + 4*CW^2)*MW^2*dZfL1[3, j1, j2])*IndexDelta[o1, o2])/
     (CW^3*MW^2*SW)}, {(((-2*I)/3)*EL*SW*IndexDelta[j1, j2]*
      IndexDelta[o1, o2])/CW, ((-I/3)*EL*SW*(Conjugate[dZfR1[3, j2, j1]] + 
       dZfR1[3, j1, j2])*IndexDelta[o1, o2])/CW}}, 
 C[-F[4, {j1, o1}], F[4, {j2, o2}], V[20]] == 
  {{((-I/6)*(1 + 2*CW^2)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2])/(CW*SW), 
    ((-I/12)*EL*(CW^2*(1 + 2*CW^2)*MW^2*Conjugate[dZfL1[4, j2, j1]] + 
       CW^2*(1 + 2*CW^2)*MW^2*dZfL1[4, j1, j2])*IndexDelta[o1, o2])/
     (CW^3*MW^2*SW)}, {((I/3)*EL*SW*IndexDelta[j1, j2]*IndexDelta[o1, o2])/
     CW, ((I/6)*EL*SW*(Conjugate[dZfR1[4, j2, j1]] + dZfR1[4, j1, j2])*
      IndexDelta[o1, o2])/CW}}, C[-F[1, {j1}], F[2, {j2}], -V[30]] == 
  {{(I*EL*IndexDelta[j1, j2])/(Sqrt[2]*SW), 
    ((I/2)*EL*(Conjugate[dZfL1[1, j1, j1]] + dZfL1[2, j1, j1])*
      IndexDelta[j1, j2])/(Sqrt[2]*SW), 
    ((-I/8)*EL*(8*dSW2*MW^4 + dMWsq1^2*SW + 4*dMWsq1*dZe1*MW^2*SW - 
       2*dMWsq1*dZH1*MW^2*SW - 8*dZe2*MW^4*SW - 4*dZe1*dZH1*MW^4*SW + 
       dZH1^2*MW^4*SW - 4*dZW2*MW^4*SW + MW^4*SW*Conjugate[dZfL1[1, j1, j1]]^
         2 - 4*MW^4*SW*Conjugate[dZfL2[1, j1, j1]] - 
       2*MW^4*SW*Conjugate[dZfL1[1, j1, j1]]*dZfL1[2, j1, j1] + 
       MW^4*SW*dZfL1[2, j1, j1]^2 - 4*MW^4*SW*dZfL2[2, j1, j1])*
      IndexDelta[j1, j2])/(Sqrt[2]*MW^4*SW^2)}, {0, 0, 0}}, 
 C[-F[2, {j1}], F[1, {j2}], V[30]] == 
  {{(I*EL*IndexDelta[j1, j2])/(Sqrt[2]*SW), 
    ((I/2)*EL*(Conjugate[dZfL1[2, j1, j1]] + dZfL1[1, j1, j1])*
      IndexDelta[j1, j2])/(Sqrt[2]*SW), 
    ((-I/8)*EL*(8*dSW2*MW^4 + dMWsq1^2*SW + 4*dMWsq1*dZe1*MW^2*SW - 
       2*dMWsq1*dZH1*MW^2*SW - 8*dZe2*MW^4*SW - 4*dZe1*dZH1*MW^4*SW + 
       dZH1^2*MW^4*SW - 4*dZW2*MW^4*SW + MW^4*SW*Conjugate[dZfL1[2, j1, j1]]^
         2 - 4*MW^4*SW*Conjugate[dZfL2[2, j1, j1]] - 
       2*MW^4*SW*Conjugate[dZfL1[2, j1, j1]]*dZfL1[1, j1, j1] + 
       MW^4*SW*dZfL1[1, j1, j1]^2 - 4*MW^4*SW*dZfL2[1, j1, j1])*
      IndexDelta[j1, j2])/(Sqrt[2]*MW^4*SW^2)}, {0, 0, 0}}, 
 C[-F[3, {j1, o1}], F[4, {j2, o2}], -V[30]] == 
  {{(I*EL*CKM[j1, j2]*IndexDelta[o1, o2])/(Sqrt[2]*SW), 
    ((I/2)*EL*IndexDelta[o1, o2]*(2*dCKM1[j1, j2] + 
       IndexSum[CKM[gn, j2]*Conjugate[dZfL1[3, gn, j1]] + 
         CKM[j1, gn]*dZfL1[4, gn, j2], {gn, 3}]))/(Sqrt[2]*SW)}, {0, 0}}, 
 C[-F[4, {j2, o2}], F[3, {j1, o1}], V[30]] == 
  {{(I*EL*Conjugate[CKM[j1, j2]]*IndexDelta[o1, o2])/(Sqrt[2]*SW), 
    ((I/2)*EL*IndexDelta[o1, o2]*(2*Conjugate[dCKM1[j1, j2]] + 
       IndexSum[Conjugate[CKM[j1, gn]]*Conjugate[dZfL1[4, gn, j2]] + 
         Conjugate[CKM[gn, j2]]*dZfL1[3, gn, j1], {gn, 3}]))/(Sqrt[2]*SW)}, 
   {0, 0}}, C[-F[2, {j1}], F[2, {j2}], S[10]] == 
  {{((-I/2)*EL*IndexDelta[j1, j2]*Mass[F[2, {j1}]])/(MW*SW), 
    ((-I/4)*EL*IndexDelta[j1, j2]*(2*dMf1[2, j1] + dZfL1[2, j1, j1]*
        Mass[F[2, {j1}]] + Conjugate[dZfR1[2, j2, j2]]*Mass[F[2, {j2}]]))/
     (MW*SW)}, {((-I/2)*EL*IndexDelta[j1, j2]*Mass[F[2, {j1}]])/(MW*SW), 
    ((-I/4)*EL*IndexDelta[j1, j2]*(2*dMf1[2, j1] + dZfR1[2, j1, j1]*
        Mass[F[2, {j1}]] + Conjugate[dZfL1[2, j2, j2]]*Mass[F[2, {j2}]]))/
     (MW*SW)}}, C[-F[3, {j1, o1}], F[3, {j2, o2}], S[10]] == 
  {{((-I/2)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[3, {j1}]])/
     (MW*SW), ((-I/4)*EL*IndexDelta[o1, o2]*
      (2*dMf1[3, j1]*IndexDelta[j1, j2] + dZfL1[3, j1, j2]*Mass[F[3, {j1}]] + 
       Conjugate[dZfR1[3, j2, j1]]*Mass[F[3, {j2}]]))/(MW*SW)}, 
   {((-I/2)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[3, {j1}]])/
     (MW*SW), ((-I/4)*EL*IndexDelta[o1, o2]*
      (2*dMf1[3, j1]*IndexDelta[j1, j2] + dZfR1[3, j1, j2]*Mass[F[3, {j1}]] + 
       Conjugate[dZfL1[3, j2, j1]]*Mass[F[3, {j2}]]))/(MW*SW)}}, 
 C[-F[4, {j1, o1}], F[4, {j2, o2}], S[10]] == 
  {{((-I/2)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[4, {j1}]])/
     (MW*SW), ((-I/4)*EL*IndexDelta[o1, o2]*
      (2*dMf1[4, j1]*IndexDelta[j1, j2] + dZfL1[4, j1, j2]*Mass[F[4, {j1}]] + 
       Conjugate[dZfR1[4, j2, j1]]*Mass[F[4, {j2}]]))/(MW*SW)}, 
   {((-I/2)*EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[4, {j1}]])/
     (MW*SW), ((-I/4)*EL*IndexDelta[o1, o2]*
      (2*dMf1[4, j1]*IndexDelta[j1, j2] + dZfR1[4, j1, j2]*Mass[F[4, {j1}]] + 
       Conjugate[dZfL1[4, j2, j1]]*Mass[F[4, {j2}]]))/(MW*SW)}}, 
 C[-F[2, {j1}], F[2, {j2}], S[20]] == 
  {{-(EL*IndexDelta[j1, j2]*Mass[F[2, {j1}]])/(2*MW*SW), 
    -(EL*IndexDelta[j1, j2]*(2*dMf1[2, j1] + dZfL1[2, j1, j1]*
         Mass[F[2, {j1}]] + Conjugate[dZfR1[2, j2, j2]]*Mass[F[2, {j2}]]))/
     (4*MW*SW)}, {(EL*IndexDelta[j1, j2]*Mass[F[2, {j1}]])/(2*MW*SW), 
    (EL*IndexDelta[j1, j2]*(2*dMf1[2, j1] + dZfR1[2, j1, j1]*
        Mass[F[2, {j1}]] + Conjugate[dZfL1[2, j2, j2]]*Mass[F[2, {j2}]]))/
     (4*MW*SW)}}, C[-F[3, {j1, o1}], F[3, {j2, o2}], S[20]] == 
  {{(EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[3, {j1}]])/(2*MW*SW), 
    (EL*IndexDelta[o1, o2]*(2*dMf1[3, j1]*IndexDelta[j1, j2] + 
       dZfL1[3, j1, j2]*Mass[F[3, {j1}]] + Conjugate[dZfR1[3, j2, j1]]*
        Mass[F[3, {j2}]]))/(4*MW*SW)}, 
   {-(EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[3, {j1}]])/(2*MW*SW), 
    -(EL*IndexDelta[o1, o2]*(2*dMf1[3, j1]*IndexDelta[j1, j2] + 
        dZfR1[3, j1, j2]*Mass[F[3, {j1}]] + Conjugate[dZfL1[3, j2, j1]]*
         Mass[F[3, {j2}]]))/(4*MW*SW)}}, 
 C[-F[4, {j1, o1}], F[4, {j2, o2}], S[20]] == 
  {{-(EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[4, {j1}]])/(2*MW*SW), 
    -(EL*IndexDelta[o1, o2]*(2*dMf1[4, j1]*IndexDelta[j1, j2] + 
        dZfL1[4, j1, j2]*Mass[F[4, {j1}]] + Conjugate[dZfR1[4, j2, j1]]*
         Mass[F[4, {j2}]]))/(4*MW*SW)}, 
   {(EL*IndexDelta[j1, j2]*IndexDelta[o1, o2]*Mass[F[4, {j1}]])/(2*MW*SW), 
    (EL*IndexDelta[o1, o2]*(2*dMf1[4, j1]*IndexDelta[j1, j2] + 
       dZfR1[4, j1, j2]*Mass[F[4, {j1}]] + Conjugate[dZfL1[4, j2, j1]]*
        Mass[F[4, {j2}]]))/(4*MW*SW)}}, 
 C[-F[3, {j1, o1}], F[4, {j2, o2}], -S[30]] == 
  {{(I*EL*CKM[j1, j2]*IndexDelta[o1, o2]*Mass[F[3, {j1}]])/(Sqrt[2]*MW*SW), 
    ((I/2)*EL*IndexDelta[o1, o2]*(2*CKM[j1, j2]*dMf1[3, j1] + 
       IndexSum[CKM[gn, j2]*Conjugate[dZfR1[3, gn, j1]]*Mass[F[3, {gn}]] + 
         CKM[j1, gn]*dZfL1[4, gn, j2]*Mass[F[3, {j1}]], {gn, 3}] + 
       2*dCKM1[j1, j2]*Mass[F[3, {j1}]]))/(Sqrt[2]*MW*SW)}, 
   {((-I)*EL*CKM[j1, j2]*IndexDelta[o1, o2]*Mass[F[4, {j2}]])/
     (Sqrt[2]*MW*SW), ((-I/2)*EL*IndexDelta[o1, o2]*
      (2*CKM[j1, j2]*dMf1[4, j2] + IndexSum[CKM[j1, gn]*dZfR1[4, gn, j2]*
          Mass[F[4, {gn}]] + CKM[gn, j2]*Conjugate[dZfL1[3, gn, j1]]*
          Mass[F[4, {j2}]], {gn, 3}] + 2*dCKM1[j1, j2]*Mass[F[4, {j2}]]))/
     (Sqrt[2]*MW*SW)}}, C[-F[4, {j2, o2}], F[3, {j1, o1}], S[30]] == 
  {{((-I)*EL*Conjugate[CKM[j1, j2]]*IndexDelta[o1, o2]*Mass[F[4, {j2}]])/
     (Sqrt[2]*MW*SW), ((-I/2)*EL*IndexDelta[o1, o2]*
      (2*Conjugate[CKM[j1, j2]]*dMf1[4, j2] + 
       IndexSum[Conjugate[CKM[j1, gn]]*Conjugate[dZfR1[4, gn, j2]]*
          Mass[F[4, {gn}]] + Conjugate[CKM[gn, j2]]*dZfL1[3, gn, j1]*
          Mass[F[4, {j2}]], {gn, 3}] + 2*Conjugate[dCKM1[j1, j2]]*
        Mass[F[4, {j2}]]))/(Sqrt[2]*MW*SW)}, 
   {(I*EL*Conjugate[CKM[j1, j2]]*IndexDelta[o1, o2]*Mass[F[3, {j1}]])/
     (Sqrt[2]*MW*SW), ((I/2)*EL*IndexDelta[o1, o2]*
      (2*Conjugate[CKM[j1, j2]]*dMf1[3, j2]*Mass[F[3, {j1}]] + 
       (IndexSum[Conjugate[CKM[gn, j2]]*dZfR1[3, gn, j1]*Mass[F[3, {gn}]] + 
           Conjugate[CKM[j1, gn]]*Conjugate[dZfL1[4, gn, j2]]*
            Mass[F[3, {j1}]], {gn, 3}] + 2*Conjugate[dCKM1[j1, j2]]*
          Mass[F[3, {j1}]])*Mass[F[3, {j2}]]))/(Sqrt[2]*MW*SW*
      Mass[F[3, {j2}]])}}, C[-F[1, {j1}], F[2, {j2}], -S[30]] == 
  {{0, 0}, {((-I)*EL*IndexDelta[j1, j2]*Mass[F[2, {j1}]])/(Sqrt[2]*MW*SW), 
    ((-I/2)*EL*IndexDelta[j1, j2]*(2*dMf1[2, j1] + 
       (Conjugate[dZfL1[1, j1, j1]] + dZfR1[2, j1, j1])*Mass[F[2, {j1}]]))/
     (Sqrt[2]*MW*SW)}}, C[-F[2, {j1}], F[1, {j2}], S[30]] == 
  {{((-I)*EL*IndexDelta[j1, j2]*Mass[F[2, {j1}]])/(Sqrt[2]*MW*SW), 
    ((-I/2)*EL*IndexDelta[j1, j2]*(2*dMf1[2, j1] + 
       (Conjugate[dZfR1[2, j1, j1]] + dZfL1[1, j1, j1])*Mass[F[2, {j1}]]))/
     (Sqrt[2]*MW*SW)}, {0, 0}}, C[-V[30], V[30]] == 
  {{0, I*dZW1}, {0, I*(dMWsq1 + dZW1*MW^2)}, {0, (-I)*dZW1}}, 
 C[V[20], V[20]] == {{0, I*dZZZ1}, {0, I*(dMZsq1 + dZZZ1*MZ^2)}, 
   {0, (-I)*dZZZ1}}, C[V[10], V[10]] == {{0, I*dZAA1}, {0, 0}, 
   {0, (-I)*dZAA1}}, C[V[10], V[20]] == {{0, (I/2)*dZAZ1}, {0, 0}, 
   {0, (-I/2)*dZAZ1}}, C[S[30], -V[30]] == {{0, 0}, {0, I*dZH1*MW}}, 
 C[-S[30], V[30]] == {{0, 0}, {0, (-I)*dZH1*MW}}, 
 C[S[20], V[20]] == {{0, 0}, {0, -(dZH1*MZ)}}, 
 C[S[10], S[10]] == {{0, (-I)*dZH1}, {0, I*(-dMHsq1 - dZH1*MH^2)}}, 
 C[S[20], S[20]] == {{0, (-I)*dZH1}, {0, ((I/2)*dTH1*EL)/(MW*SW)}}, 
 C[S[30], -S[30]] == {{0, (-I)*dZH1}, {0, ((I/2)*dTH1*EL)/(MW*SW)}}, 
 C[V[10], -V[3], V[3]] == {{(-I)*EL}, {((-I)*EL)/GaugeXi[Q]}, 
   {(I*EL)/GaugeXi[Q]}, {0}}, C[-V[30], V[3], V[1]] == 
  {{(-I)*EL}, {((-I)*EL)/GaugeXi[Q]}, {(I*EL)/GaugeXi[Q]}, {0}}, 
 C[V[30], V[1], -V[3]] == {{(-I)*EL}, {((-I)*EL)/GaugeXi[Q]}, 
   {(I*EL)/GaugeXi[Q]}, {0}}, C[V[20], -V[3], V[3]] == 
  {{(I*CW*EL)/SW}, {(I*CW*EL)/(SW*GaugeXi[Q])}, 
   {((-I)*CW*EL)/(SW*GaugeXi[Q])}, {0}}, C[-V[30], V[3], V[2]] == 
  {{(I*CW*EL)/SW}, {(I*CW*EL)/(SW*GaugeXi[Q])}, 
   {((-I)*CW*EL)/(SW*GaugeXi[Q])}, {0}}, C[V[30], V[2], -V[3]] == 
  {{(I*CW*EL)/SW}, {(I*CW*EL)/(SW*GaugeXi[Q])}, 
   {((-I)*CW*EL)/(SW*GaugeXi[Q])}, {0}}, C[S[10], S[1], S[1]] == 
  {{(((-3*I)/2)*EL*MH^2)/(MW*SW)}}, C[S[10], S[2], S[2]] == 
  {{I*EL*(-MH^2/(2*MW*SW) - (MW*GaugeXi[Q])/(CW^2*SW))}}, 
 C[S[1], S[20], S[2]] == 
  {{I*EL*(-MH^2/(2*MW*SW) + (MW*GaugeXi[Q])/(CW^2*SW))}}, 
 C[S[10], -S[3], S[3]] == {{I*EL*(-MH^2/(2*MW*SW) - (MW*GaugeXi[Q])/SW)}}, 
 C[S[1], -S[30], S[3]] == 
  {{I*EL*(-MH^2/(2*MW*SW) + (MW*GaugeXi[Q])/(2*SW))}}, 
 C[S[1], -S[3], S[30]] == 
  {{I*EL*(-MH^2/(2*MW*SW) + (MW*GaugeXi[Q])/(2*SW))}}, 
 C[S[2], -S[30], S[3]] == {{(EL*MW*SW*GaugeXi[Q])/(2*CW^2)}}, 
 C[S[2], -S[3], S[30]] == {{-(EL*MW*SW*GaugeXi[Q])/(2*CW^2)}}, 
 C[S[30], S[1], -V[3]] == {{((-I)*EL)/SW}, {0}}, 
 C[-S[30], S[1], V[3]] == {{(I*EL)/SW}, {0}}, 
 C[S[10], S[3], -V[3]] == {{(I*EL)/SW}, {0}}, 
 C[S[10], -S[3], V[3]] == {{((-I)*EL)/SW}, {0}}, 
 C[S[30], S[2], -V[3]] == {{EL/SW}, {0}}, C[-S[30], S[2], V[3]] == 
  {{EL/SW}, {0}}, C[S[20], S[3], -V[3]] == {{-(EL/SW)}, {0}}, 
 C[S[20], -S[3], V[3]] == {{-(EL/SW)}, {0}}, C[S[30], -S[3], V[2]] == 
  {{((-I)*EL*(CW^2 - SW^2))/(CW*SW)}, {0}}, C[-S[30], S[3], V[2]] == 
  {{(I*EL*(CW^2 - SW^2))/(CW*SW)}, {0}}, C[S[30], -S[3], V[1]] == 
  {{(2*I)*EL}, {0}}, C[-S[30], S[3], V[1]] == {{(-2*I)*EL}, {0}}, 
 C[S[20], S[1], V[2]] == {{EL/(CW*SW)}, {0}}, 
 C[S[10], S[2], V[2]] == {{-(EL/(CW*SW))}, {0}}, 
 C[S[1], -V[30], V[3]] == {{(I*EL*MW)/SW}}, C[S[1], V[30], -V[3]] == 
  {{(I*EL*MW)/SW}}, C[S[1], V[20], V[2]] == {{(I*EL*MW)/(CW^2*SW)}}, 
 C[-S[3], V[20], V[3]] == {{((-I)*EL*MW)/(CW*SW)}}, 
 C[S[3], V[20], -V[3]] == {{((-I)*EL*MW)/(CW*SW)}}, 
 C[-S[3], V[30], V[2]] == {{(I*EL*MW*(CW^2 - SW^2))/(CW*SW)}}, 
 C[S[3], -V[30], V[2]] == {{(I*EL*MW*(CW^2 - SW^2))/(CW*SW)}}, 
 C[-S[3], V[30], V[1]] == {{(-2*I)*EL*MW}}, C[S[3], -V[30], V[1]] == 
  {{(-2*I)*EL*MW}}, C[S[2], -V[30], V[3]] == {{(EL*MW)/SW}}, 
 C[S[2], V[30], -V[3]] == {{-((EL*MW)/SW)}}, C[-U[4], U[2], -V[30]] == 
  {{(I*CW*EL)/SW}, {((-I)*CW*EL)/SW}}, C[-U[3], U[2], V[30]] == 
  {{((-I)*CW*EL)/SW}, {(I*CW*EL)/SW}}, C[-U[4], U[1], -V[30]] == 
  {{(-I)*EL}, {I*EL}}, C[-U[3], U[1], V[30]] == {{I*EL}, {(-I)*EL}}, 
 C[-U[2], U[3], -V[30]] == {{((-I)*CW*EL)/SW}, {(I*CW*EL)/SW}}, 
 C[-U[2], U[4], V[30]] == {{(I*CW*EL)/SW}, {((-I)*CW*EL)/SW}}, 
 C[-U[1], U[3], -V[30]] == {{I*EL}, {(-I)*EL}}, 
 C[-U[1], U[4], V[30]] == {{(-I)*EL}, {I*EL}}, 
 C[-U[4], U[4], V[10]] == {{I*EL}, {(-I)*EL}}, 
 C[-U[3], U[3], V[10]] == {{(-I)*EL}, {I*EL}}, 
 C[-U[4], U[4], V[20]] == {{((-I)*CW*EL)/SW}, {(I*CW*EL)/SW}}, 
 C[-U[3], U[3], V[20]] == {{(I*CW*EL)/SW}, {((-I)*CW*EL)/SW}}, 
 C[-S[30], -U[4], U[2]] == {{(I*EL*MW*SW*GaugeXi[Q])/CW}}, 
 C[S[30], -U[3], U[2]] == {{(I*EL*MW*SW*GaugeXi[Q])/CW}}, 
 C[-S[30], -U[4], U[1]] == {{I*EL*MW*GaugeXi[Q]}}, 
 C[S[30], -U[3], U[1]] == {{I*EL*MW*GaugeXi[Q]}}, 
 C[-S[30], -U[2], U[3]] == {{(I*EL*MW*SW*GaugeXi[Q])/CW}}, 
 C[S[30], -U[2], U[4]] == {{(I*EL*MW*SW*GaugeXi[Q])/CW}}, 
 C[-S[30], -U[1], U[3]] == {{I*EL*MW*GaugeXi[Q]}}, 
 C[S[30], -U[1], U[4]] == {{I*EL*MW*GaugeXi[Q]}}, 
 C[S[10], -U[4], U[4]] == {{((-I)*EL*MW*GaugeXi[Q])/SW}}, 
 C[S[10], -U[3], U[3]] == {{((-I)*EL*MW*GaugeXi[Q])/SW}}, 
 C[S[10], -U[2], U[2]] == {{((-I)*EL*MW*GaugeXi[Q])/(CW^2*SW)}}, 
 C[V[10], V[10], -V[3], V[3]] == {{(-2*I)*EL^2}, 
   {(-I)*EL^2*(-1 + GaugeXi[Q]^(-1))}, {(-I)*EL^2*(-1 + GaugeXi[Q]^(-1))}}, 
 C[-V[30], V[30], V[1], V[1]] == {{(-2*I)*EL^2}, 
   {(-I)*EL^2*(-1 + GaugeXi[Q]^(-1))}, {(-I)*EL^2*(-1 + GaugeXi[Q]^(-1))}}, 
 C[V[20], V[20], -V[3], V[3]] == {{((-2*I)*CW^2*EL^2)/SW^2}, 
   {((-I)*CW^2*EL^2*(-1 + GaugeXi[Q]^(-1)))/SW^2}, 
   {((-I)*CW^2*EL^2*(-1 + GaugeXi[Q]^(-1)))/SW^2}}, 
 C[-V[30], V[30], V[2], V[2]] == {{((-2*I)*CW^2*EL^2)/SW^2}, 
   {((-I)*CW^2*EL^2*(-1 + GaugeXi[Q]^(-1)))/SW^2}, 
   {((-I)*CW^2*EL^2*(-1 + GaugeXi[Q]^(-1)))/SW^2}}, 
 C[V[10], V[20], -V[3], V[3]] == {{((2*I)*CW*EL^2)/SW}, 
   {(I*CW*EL^2*(-1 + GaugeXi[Q]^(-1)))/SW}, 
   {(I*CW*EL^2*(-1 + GaugeXi[Q]^(-1)))/SW}}, C[-V[30], V[30], V[1], V[2]] == 
  {{((2*I)*CW*EL^2)/SW}, {(I*CW*EL^2*(-1 + GaugeXi[Q]^(-1)))/SW}, 
   {(I*CW*EL^2*(-1 + GaugeXi[Q]^(-1)))/SW}}, C[V[30], V[30], -V[3], -V[3]] == 
  {{((2*I)*EL^2)/SW^2}, {(I*EL^2*(-1 + GaugeXi[Q]^(-1)))/SW^2}, 
   {(I*EL^2*(-1 + GaugeXi[Q]^(-1)))/SW^2}}, C[-V[30], -V[30], V[3], V[3]] == 
  {{((2*I)*EL^2)/SW^2}, {(I*EL^2*(-1 + GaugeXi[Q]^(-1)))/SW^2}, 
   {(I*EL^2*(-1 + GaugeXi[Q]^(-1)))/SW^2}}, C[-V[30], V[10], V[3], V[1]] == 
  {{I*EL^2}, {(-2*I)*EL^2}, {(-I)*EL^2*(-1 - GaugeXi[Q]^(-1))}}, 
 C[V[30], V[10], -V[3], V[1]] == {{I*EL^2}, {(-2*I)*EL^2}, 
   {(-I)*EL^2*(-1 - GaugeXi[Q]^(-1))}}, C[-V[30], V[20], V[3], V[2]] == 
  {{(I*CW^2*EL^2)/SW^2}, {((-2*I)*CW^2*EL^2)/SW^2}, 
   {((-I)*CW^2*EL^2*(-1 - GaugeXi[Q]^(-1)))/SW^2}}, 
 C[V[30], V[20], -V[3], V[2]] == {{(I*CW^2*EL^2)/SW^2}, 
   {((-2*I)*CW^2*EL^2)/SW^2}, {((-I)*CW^2*EL^2*(-1 - GaugeXi[Q]^(-1)))/
     SW^2}}, C[-V[30], V[10], V[3], V[2]] == {{((-I)*CW*EL^2)/SW}, 
   {((2*I)*CW*EL^2)/SW}, {(I*CW*EL^2*(-1 - GaugeXi[Q]^(-1)))/SW}}, 
 C[V[30], V[10], -V[3], V[2]] == {{((-I)*CW*EL^2)/SW}, {((2*I)*CW*EL^2)/SW}, 
   {(I*CW*EL^2*(-1 - GaugeXi[Q]^(-1)))/SW}}, C[-V[30], V[20], V[3], V[1]] == 
  {{((-I)*CW*EL^2)/SW}, {((2*I)*CW*EL^2)/SW}, 
   {(I*CW*EL^2*(-1 - GaugeXi[Q]^(-1)))/SW}}, C[V[30], V[20], -V[3], V[1]] == 
  {{((-I)*CW*EL^2)/SW}, {((2*I)*CW*EL^2)/SW}, 
   {(I*CW*EL^2*(-1 - GaugeXi[Q]^(-1)))/SW}}, C[-V[30], V[30], -V[3], V[3]] == 
  {{((-I)*EL^2)/SW^2}, {((2*I)*EL^2)/SW^2}, 
   {(I*EL^2*(-1 - GaugeXi[Q]^(-1)))/SW^2}}, C[S[1], S[1], -S[30], S[30]] == 
  {{I*EL^2*(-MH^2/(4*MW^2*SW^2) - GaugeXi[Q]/(2*SW^2))}}, 
 C[S[10], S[10], -S[3], S[3]] == 
  {{I*EL^2*(-MH^2/(4*MW^2*SW^2) - GaugeXi[Q]/(2*SW^2))}}, 
 C[S[2], S[2], -S[30], S[30]] == 
  {{I*EL^2*(-MH^2/(4*MW^2*SW^2) - GaugeXi[Q]/(2*SW^2))}}, 
 C[S[20], S[20], -S[3], S[3]] == 
  {{I*EL^2*(-MH^2/(4*MW^2*SW^2) - GaugeXi[Q]/(2*SW^2))}}, 
 C[S[10], S[1], -S[30], S[3]] == 
  {{I*EL^2*(-MH^2/(4*MW^2*SW^2) + GaugeXi[Q]/(4*SW^2))}}, 
 C[S[10], S[1], -S[3], S[30]] == 
  {{I*EL^2*(-MH^2/(4*MW^2*SW^2) + GaugeXi[Q]/(4*SW^2))}}, 
 C[S[20], S[2], -S[30], S[3]] == 
  {{I*EL^2*(-MH^2/(4*MW^2*SW^2) + GaugeXi[Q]/(4*SW^2))}}, 
 C[S[20], S[2], -S[3], S[30]] == 
  {{I*EL^2*(-MH^2/(4*MW^2*SW^2) + GaugeXi[Q]/(4*SW^2))}}, 
 C[S[10], S[2], -S[30], S[3]] == {{(EL^2*GaugeXi[Q])/(4*CW^2)}}, 
 C[S[1], S[20], -S[30], S[3]] == {{-(EL^2*GaugeXi[Q])/(4*CW^2)}}, 
 C[S[10], S[2], -S[3], S[30]] == {{-(EL^2*GaugeXi[Q])/(4*CW^2)}}, 
 C[S[1], S[20], -S[3], S[30]] == {{(EL^2*GaugeXi[Q])/(4*CW^2)}}, 
 C[S[1], S[1], S[20], S[20]] == 
  {{I*EL^2*(-MH^2/(4*MW^2*SW^2) - GaugeXi[Q]/(2*CW^2*SW^2))}}, 
 C[S[10], S[10], S[2], S[2]] == 
  {{I*EL^2*(-MH^2/(4*MW^2*SW^2) - GaugeXi[Q]/(2*CW^2*SW^2))}}, 
 C[S[10], S[1], S[20], S[2]] == 
  {{I*EL^2*(-MH^2/(4*MW^2*SW^2) + GaugeXi[Q]/(4*CW^2*SW^2))}}, 
 C[-S[30], -S[30], S[3], S[3]] == 
  {{I*EL^2*(-MH^2/(2*MW^2*SW^2) + GaugeXi[Q]/(2*CW^2*SW^2))}}, 
 C[S[30], S[30], -S[3], -S[3]] == 
  {{I*EL^2*(-MH^2/(2*MW^2*SW^2) + GaugeXi[Q]/(2*CW^2*SW^2))}}, 
 C[S[30], -S[30], S[3], -S[3]] == 
  {{I*EL^2*(-MH^2/(2*MW^2*SW^2) - GaugeXi[Q]/(4*CW^2*SW^2))}}, 
 C[S[10], S[10], S[1], S[1]] == {{(((-3*I)/4)*EL^2*MH^2)/(MW^2*SW^2)}}, 
 C[S[20], S[20], S[2], S[2]] == {{(((-3*I)/4)*EL^2*MH^2)/(MW^2*SW^2)}}, 
 C[-S[30], S[3], V[10], V[1]] == {{(2*I)*EL^2}}, 
 C[-S[3], S[30], V[10], V[1]] == {{(2*I)*EL^2}}, 
 C[-S[30], S[3], V[10], V[2]] == {{((-I)*EL^2*(CW^2 - SW^2))/(CW*SW)}}, 
 C[-S[3], S[30], V[10], V[2]] == {{((-I)*EL^2*(CW^2 - SW^2))/(CW*SW)}}, 
 C[-S[30], S[3], V[1], V[20]] == {{((-I)*EL^2*(CW^2 - SW^2))/(CW*SW)}}, 
 C[-S[3], S[30], V[1], V[20]] == {{((-I)*EL^2*(CW^2 - SW^2))/(CW*SW)}}, 
 C[-S[30], S[3], V[20], V[2]] == 
  {{((I/2)*EL^2*(CW^2 - SW^2)^2)/(CW^2*SW^2)}}, 
 C[-S[3], S[30], V[20], V[2]] == 
  {{((I/2)*EL^2*(CW^2 - SW^2)^2)/(CW^2*SW^2)}}, 
 C[S[10], S[1], V[20], V[2]] == {{((I/2)*EL^2)/(CW^2*SW^2)}}, 
 C[S[20], S[2], V[20], V[2]] == {{((I/2)*EL^2)/(CW^2*SW^2)}}, 
 C[S[10], S[1], -V[30], V[3]] == {{((I/2)*EL^2)/SW^2}}, 
 C[S[10], S[1], -V[3], V[30]] == {{((I/2)*EL^2)/SW^2}}, 
 C[S[20], S[2], -V[30], V[3]] == {{((I/2)*EL^2)/SW^2}}, 
 C[S[20], S[2], -V[3], V[30]] == {{((I/2)*EL^2)/SW^2}}, 
 C[S[1], S[30], V[10], -V[3]] == {{((-I)*EL^2)/SW}}, 
 C[S[1], -S[30], V[10], V[3]] == {{((-I)*EL^2)/SW}}, 
 C[S[10], S[3], V[1], -V[30]] == {{((-I)*EL^2)/SW}}, 
 C[S[10], -S[3], V[1], V[30]] == {{((-I)*EL^2)/SW}}, 
 C[S[2], S[30], V[10], -V[3]] == {{EL^2/SW}}, C[S[2], -S[30], V[10], V[3]] == 
  {{-(EL^2/SW)}}, C[S[20], S[3], V[1], -V[30]] == {{EL^2/SW}}, 
 C[S[20], -S[3], V[1], V[30]] == {{-(EL^2/SW)}}, 
 C[S[10], S[3], V[20], -V[3]] == {{((-I/2)*EL^2)/(CW*SW^2)}}, 
 C[S[10], -S[3], V[20], V[3]] == {{((-I/2)*EL^2)/(CW*SW^2)}}, 
 C[S[1], S[30], V[2], -V[30]] == {{((-I/2)*EL^2)/(CW*SW^2)}}, 
 C[S[1], -S[30], V[2], V[30]] == {{((-I/2)*EL^2)/(CW*SW^2)}}, 
 C[S[1], S[30], V[20], -V[3]] == {{((I/2)*EL^2*(CW^2 - SW^2))/(CW*SW^2)}}, 
 C[S[1], -S[30], V[20], V[3]] == {{((I/2)*EL^2*(CW^2 - SW^2))/(CW*SW^2)}}, 
 C[S[10], S[3], V[2], -V[30]] == {{((I/2)*EL^2*(CW^2 - SW^2))/(CW*SW^2)}}, 
 C[S[10], -S[3], V[2], V[30]] == {{((I/2)*EL^2*(CW^2 - SW^2))/(CW*SW^2)}}, 
 C[S[20], S[3], V[20], -V[3]] == {{EL^2/(2*CW*SW^2)}}, 
 C[S[20], -S[3], V[20], V[3]] == {{-EL^2/(2*CW*SW^2)}}, 
 C[S[2], S[30], V[2], -V[30]] == {{EL^2/(2*CW*SW^2)}}, 
 C[S[2], -S[30], V[2], V[30]] == {{-EL^2/(2*CW*SW^2)}}, 
 C[S[2], S[30], V[20], -V[3]] == {{-(EL^2*(CW^2 - SW^2))/(2*CW*SW^2)}}, 
 C[S[2], -S[30], V[20], V[3]] == {{(EL^2*(CW^2 - SW^2))/(2*CW*SW^2)}}, 
 C[S[20], S[3], V[2], -V[30]] == {{-(EL^2*(CW^2 - SW^2))/(2*CW*SW^2)}}, 
 C[S[20], -S[3], V[2], V[30]] == {{(EL^2*(CW^2 - SW^2))/(2*CW*SW^2)}}, 
 C[-S[3], S[30], -V[30], V[3]] == {{(I*EL^2)/SW^2}}, 
 C[-S[30], S[3], -V[3], V[30]] == {{(I*EL^2)/SW^2}}, 
 C[S[10], S[2], -V[30], V[3]] == {{EL^2/(2*SW^2)}}, 
 C[S[10], S[2], -V[3], V[30]] == {{-EL^2/(2*SW^2)}}, 
 C[S[1], S[20], -V[30], V[3]] == {{-EL^2/(2*SW^2)}}, 
 C[S[1], S[20], -V[3], V[30]] == {{EL^2/(2*SW^2)}}, 
 C[-U[4], U[4], V[10], V[1]] == {{I*EL^2}}, C[-U[3], U[3], V[10], V[1]] == 
  {{I*EL^2}}, C[-U[4], U[4], V[10], V[2]] == {{((-I)*CW*EL^2)/SW}}, 
 C[-U[3], U[3], V[10], V[2]] == {{((-I)*CW*EL^2)/SW}}, 
 C[-U[4], U[4], V[20], V[1]] == {{((-I)*CW*EL^2)/SW}}, 
 C[-U[3], U[3], V[20], V[1]] == {{((-I)*CW*EL^2)/SW}}, 
 C[-U[4], U[4], V[20], V[2]] == {{(I*CW^2*EL^2)/SW^2}}, 
 C[-U[3], U[3], V[20], V[2]] == {{(I*CW^2*EL^2)/SW^2}}, 
 C[-U[4], U[1], V[10], -V[3]] == {{(-I)*EL^2}}, 
 C[-U[3], U[1], V[10], V[3]] == {{(-I)*EL^2}}, 
 C[-U[4], U[1], V[20], -V[3]] == {{(I*CW*EL^2)/SW}}, 
 C[-U[3], U[1], V[20], V[3]] == {{(I*CW*EL^2)/SW}}, 
 C[-U[4], U[2], V[10], -V[3]] == {{(I*CW*EL^2)/SW}}, 
 C[-U[3], U[2], V[10], V[3]] == {{(I*CW*EL^2)/SW}}, 
 C[-U[4], U[2], V[20], -V[3]] == {{((-I)*CW^2*EL^2)/SW^2}}, 
 C[-U[3], U[2], V[20], V[3]] == {{((-I)*CW^2*EL^2)/SW^2}}, 
 C[-U[1], U[3], -V[30], V[1]] == {{(-I)*EL^2}}, 
 C[-U[1], U[4], V[30], V[1]] == {{(-I)*EL^2}}, 
 C[-U[1], U[3], -V[30], V[2]] == {{(I*CW*EL^2)/SW}}, 
 C[-U[1], U[4], V[30], V[2]] == {{(I*CW*EL^2)/SW}}, 
 C[-U[2], U[3], -V[30], V[1]] == {{(I*CW*EL^2)/SW}}, 
 C[-U[2], U[4], V[30], V[1]] == {{(I*CW*EL^2)/SW}}, 
 C[-U[2], U[3], -V[30], V[2]] == {{((-I)*CW^2*EL^2)/SW^2}}, 
 C[-U[2], U[4], V[30], V[2]] == {{((-I)*CW^2*EL^2)/SW^2}}, 
 C[-U[1], U[1], -V[30], V[3]] == {{I*EL^2}}, C[-U[1], U[1], V[30], -V[3]] == 
  {{I*EL^2}}, C[-U[2], U[1], -V[30], V[3]] == {{((-I)*CW*EL^2)/SW}}, 
 C[-U[2], U[1], V[30], -V[3]] == {{((-I)*CW*EL^2)/SW}}, 
 C[-U[1], U[2], -V[30], V[3]] == {{((-I)*CW*EL^2)/SW}}, 
 C[-U[1], U[2], V[30], -V[3]] == {{((-I)*CW*EL^2)/SW}}, 
 C[-U[2], U[2], -V[30], V[3]] == {{(I*CW^2*EL^2)/SW^2}}, 
 C[-U[2], U[2], V[30], -V[3]] == {{(I*CW^2*EL^2)/SW^2}}, 
 C[-U[4], U[4], -V[30], V[3]] == {{(I*EL^2)/SW^2}}, 
 C[-U[3], U[3], V[30], -V[3]] == {{(I*EL^2)/SW^2}}, 
 C[-U[4], U[3], -V[30], -V[3]] == {{((-I)*EL^2)/SW^2}}, 
 C[-U[3], U[4], V[30], V[3]] == {{((-I)*EL^2)/SW^2}}, 
 C[-U[4], U[4], V[10], V[10]] == {{(2*I)*EL^2}}, 
 C[-U[3], U[3], V[10], V[10]] == {{(2*I)*EL^2}}, 
 C[-U[4], U[4], V[10], V[20]] == {{((-2*I)*CW*EL^2)/SW}}, 
 C[-U[3], U[3], V[10], V[20]] == {{((-2*I)*CW*EL^2)/SW}}, 
 C[-U[4], U[4], V[20], V[10]] == {{((-2*I)*CW*EL^2)/SW}}, 
 C[-U[3], U[3], V[20], V[10]] == {{((-2*I)*CW*EL^2)/SW}}, 
 C[-U[4], U[4], V[20], V[20]] == {{((2*I)*CW^2*EL^2)/SW^2}}, 
 C[-U[3], U[3], V[20], V[20]] == {{((2*I)*CW^2*EL^2)/SW^2}}, 
 C[-U[4], U[1], V[10], -V[30]] == {{(-I)*EL^2}}, 
 C[-U[3], U[1], V[10], V[30]] == {{(-I)*EL^2}}, 
 C[-U[4], U[1], V[20], -V[30]] == {{(I*CW*EL^2)/SW}}, 
 C[-U[3], U[1], V[20], V[30]] == {{(I*CW*EL^2)/SW}}, 
 C[-U[4], U[2], V[10], -V[30]] == {{(I*CW*EL^2)/SW}}, 
 C[-U[3], U[2], V[10], V[30]] == {{(I*CW*EL^2)/SW}}, 
 C[-U[4], U[2], V[20], -V[30]] == {{((-I)*CW^2*EL^2)/SW^2}}, 
 C[-U[3], U[2], V[20], V[30]] == {{((-I)*CW^2*EL^2)/SW^2}}, 
 C[-U[1], U[3], -V[30], V[10]] == {{(-I)*EL^2}}, 
 C[-U[1], U[4], V[30], V[10]] == {{(-I)*EL^2}}, 
 C[-U[1], U[3], -V[30], V[20]] == {{(I*CW*EL^2)/SW}}, 
 C[-U[1], U[4], V[30], V[20]] == {{(I*CW*EL^2)/SW}}, 
 C[-U[2], U[3], -V[30], V[10]] == {{(I*CW*EL^2)/SW}}, 
 C[-U[2], U[4], V[30], V[10]] == {{(I*CW*EL^2)/SW}}, 
 C[-U[2], U[3], -V[30], V[20]] == {{((-I)*CW^2*EL^2)/SW^2}}, 
 C[-U[2], U[4], V[30], V[20]] == {{((-I)*CW^2*EL^2)/SW^2}}, 
 C[-U[1], U[1], V[30], -V[30]] == {{(2*I)*EL^2}}, 
 C[-U[2], U[1], V[30], -V[30]] == {{((-2*I)*CW*EL^2)/SW}}, 
 C[-U[1], U[2], V[30], -V[30]] == {{((-2*I)*CW*EL^2)/SW}}, 
 C[-U[2], U[2], V[30], -V[30]] == {{((2*I)*CW^2*EL^2)/SW^2}}, 
 C[-U[4], U[4], -V[30], V[30]] == {{(I*EL^2)/SW^2}}, 
 C[-U[3], U[3], V[30], -V[30]] == {{(I*EL^2)/SW^2}}, 
 C[-U[4], U[3], -V[30], -V[30]] == {{((-2*I)*EL^2)/SW^2}}, 
 C[-U[3], U[4], V[30], V[30]] == {{((-2*I)*EL^2)/SW^2}}, 
 C[-S[30], S[3], -U[4], U[4]] == {{((-I/2)*EL^2*GaugeXi[Q])/SW^2}}, 
 C[-S[3], S[30], -U[3], U[3]] == {{((-I/2)*EL^2*GaugeXi[Q])/SW^2}}, 
 C[S[10], S[1], -U[4], U[4]] == {{((-I/4)*EL^2*GaugeXi[Q])/SW^2}}, 
 C[S[10], S[1], -U[3], U[3]] == {{((-I/4)*EL^2*GaugeXi[Q])/SW^2}}, 
 C[S[20], S[2], -U[4], U[4]] == {{((-I/4)*EL^2*GaugeXi[Q])/SW^2}}, 
 C[S[20], S[2], -U[3], U[3]] == {{((-I/4)*EL^2*GaugeXi[Q])/SW^2}}, 
 C[S[10], S[1], -U[2], U[2]] == {{((-I/4)*EL^2*GaugeXi[Q])/(CW^2*SW^2)}}, 
 C[S[20], S[2], -U[2], U[2]] == {{((-I/4)*EL^2*GaugeXi[Q])/(CW^2*SW^2)}}, 
 C[S[10], S[2], -U[4], U[4]] == {{(EL^2*GaugeXi[Q])/(4*SW^2)}}, 
 C[S[10], S[2], -U[3], U[3]] == {{-(EL^2*GaugeXi[Q])/(4*SW^2)}}, 
 C[S[1], S[20], -U[4], U[4]] == {{-(EL^2*GaugeXi[Q])/(4*SW^2)}}, 
 C[S[1], S[20], -U[3], U[3]] == {{(EL^2*GaugeXi[Q])/(4*SW^2)}}, 
 C[S[30], -S[3], -U[2], U[2]] == 
  {{((I/4)*EL^2*(CW^2 - SW^2)^2*GaugeXi[Q])/(CW^2*SW^2)}}, 
 C[S[3], -S[30], -U[2], U[2]] == 
  {{((I/4)*EL^2*(CW^2 - SW^2)^2*GaugeXi[Q])/(CW^2*SW^2)}}, 
 C[S[3], -S[30], -U[2], U[1]] == 
  {{((I/2)*EL^2*(CW^2 - SW^2)*GaugeXi[Q])/(CW*SW)}}, 
 C[S[30], -S[3], -U[2], U[1]] == 
  {{((I/2)*EL^2*(CW^2 - SW^2)*GaugeXi[Q])/(CW*SW)}}, 
 C[S[3], -S[30], -U[1], U[2]] == 
  {{((I/2)*EL^2*(CW^2 - SW^2)*GaugeXi[Q])/(CW*SW)}}, 
 C[S[30], -S[3], -U[1], U[2]] == 
  {{((I/2)*EL^2*(CW^2 - SW^2)*GaugeXi[Q])/(CW*SW)}}, 
 C[S[3], -S[30], -U[1], U[1]] == {{(-I)*EL^2*GaugeXi[Q]}}, 
 C[S[30], -S[3], -U[1], U[1]] == {{(-I)*EL^2*GaugeXi[Q]}}, 
 C[-S[30], S[1], -U[4], U[2]] == {{((I/4)*EL^2*GaugeXi[Q])/(CW*SW^2)}}, 
 C[S[30], S[1], -U[3], U[2]] == {{((I/4)*EL^2*GaugeXi[Q])/(CW*SW^2)}}, 
 C[-S[3], S[10], -U[2], U[3]] == {{((I/4)*EL^2*GaugeXi[Q])/(CW*SW^2)}}, 
 C[S[3], S[10], -U[2], U[4]] == {{((I/4)*EL^2*GaugeXi[Q])/(CW*SW^2)}}, 
 C[-S[3], S[10], -U[4], U[2]] == 
  {{((-I/4)*EL^2*(CW^2 - SW^2)*GaugeXi[Q])/(CW*SW^2)}}, 
 C[S[3], S[10], -U[3], U[2]] == 
  {{((-I/4)*EL^2*(CW^2 - SW^2)*GaugeXi[Q])/(CW*SW^2)}}, 
 C[-S[30], S[1], -U[2], U[3]] == 
  {{((-I/4)*EL^2*(CW^2 - SW^2)*GaugeXi[Q])/(CW*SW^2)}}, 
 C[S[30], S[1], -U[2], U[4]] == 
  {{((I/4)*EL^2*(CW^2 - SW^2)*GaugeXi[Q])/(CW*SW^2)}}, 
 C[-S[30], S[1], -U[4], U[1]] == {{((I/2)*EL^2*GaugeXi[Q])/SW}}, 
 C[S[30], S[1], -U[3], U[1]] == {{((I/2)*EL^2*GaugeXi[Q])/SW}}, 
 C[-S[3], S[10], -U[1], U[3]] == {{((I/2)*EL^2*GaugeXi[Q])/SW}}, 
 C[S[3], S[10], -U[1], U[4]] == {{((I/2)*EL^2*GaugeXi[Q])/SW}}, 
 C[-S[30], S[2], -U[4], U[2]] == {{(EL^2*GaugeXi[Q])/(4*CW*SW^2)}}, 
 C[S[30], S[2], -U[3], U[2]] == {{-(EL^2*GaugeXi[Q])/(4*CW*SW^2)}}, 
 C[-S[3], S[20], -U[2], U[3]] == {{(EL^2*GaugeXi[Q])/(4*CW*SW^2)}}, 
 C[S[3], S[20], -U[2], U[4]] == {{-(EL^2*GaugeXi[Q])/(4*CW*SW^2)}}, 
 C[-S[3], S[20], -U[4], U[2]] == 
  {{-(EL^2*(CW^2 - SW^2)*GaugeXi[Q])/(4*CW*SW^2)}}, 
 C[S[3], S[20], -U[3], U[2]] == 
  {{(EL^2*(CW^2 - SW^2)*GaugeXi[Q])/(4*CW*SW^2)}}, 
 C[-S[30], S[2], -U[2], U[3]] == 
  {{-(EL^2*(CW^2 - SW^2)*GaugeXi[Q])/(4*CW*SW^2)}}, 
 C[S[30], S[2], -U[2], U[4]] == 
  {{(EL^2*(CW^2 - SW^2)*GaugeXi[Q])/(4*CW*SW^2)}}, 
 C[-S[3], S[20], -U[4], U[1]] == {{(EL^2*GaugeXi[Q])/(2*SW)}}, 
 C[S[3], S[20], -U[3], U[1]] == {{-(EL^2*GaugeXi[Q])/(2*SW)}}, 
 C[-S[30], S[2], -U[1], U[3]] == {{(EL^2*GaugeXi[Q])/(2*SW)}}, 
 C[S[30], S[2], -U[1], U[4]] == {{-(EL^2*GaugeXi[Q])/(2*SW)}}, 
 C[-S[30], S[30], -U[4], U[4]] == {{((-I/2)*EL^2*GaugeXi[Q])/SW^2}}, 
 C[-S[30], S[30], -U[3], U[3]] == {{((-I/2)*EL^2*GaugeXi[Q])/SW^2}}, 
 C[S[10], S[10], -U[4], U[4]] == {{((-I/2)*EL^2*GaugeXi[Q])/SW^2}}, 
 C[S[10], S[10], -U[3], U[3]] == {{((-I/2)*EL^2*GaugeXi[Q])/SW^2}}, 
 C[S[20], S[20], -U[4], U[4]] == {{((-I/2)*EL^2*GaugeXi[Q])/SW^2}}, 
 C[S[20], S[20], -U[3], U[3]] == {{((-I/2)*EL^2*GaugeXi[Q])/SW^2}}, 
 C[S[10], S[10], -U[2], U[2]] == {{((-I/2)*EL^2*GaugeXi[Q])/(CW^2*SW^2)}}, 
 C[S[20], S[20], -U[2], U[2]] == {{((-I/2)*EL^2*GaugeXi[Q])/(CW^2*SW^2)}}, 
 C[S[30], -S[30], -U[2], U[2]] == 
  {{((I/2)*EL^2*(CW^2 - SW^2)^2*GaugeXi[Q])/(CW^2*SW^2)}}, 
 C[S[30], -S[30], -U[2], U[1]] == 
  {{(I*EL^2*(CW^2 - SW^2)*GaugeXi[Q])/(CW*SW)}}, 
 C[S[30], -S[30], -U[1], U[2]] == 
  {{(I*EL^2*(CW^2 - SW^2)*GaugeXi[Q])/(CW*SW)}}, 
 C[S[30], -S[30], -U[1], U[1]] == {{(-2*I)*EL^2*GaugeXi[Q]}}, 
 C[-S[30], S[10], -U[4], U[2]] == {{((I/2)*EL^2*GaugeXi[Q])/CW}}, 
 C[S[30], S[10], -U[3], U[2]] == {{((I/2)*EL^2*GaugeXi[Q])/CW}}, 
 C[-S[30], S[10], -U[2], U[3]] == {{((I/2)*EL^2*GaugeXi[Q])/CW}}, 
 C[S[30], S[10], -U[2], U[4]] == {{((I/2)*EL^2*GaugeXi[Q])/CW}}, 
 C[-S[30], S[10], -U[4], U[1]] == {{((I/2)*EL^2*GaugeXi[Q])/SW}}, 
 C[S[30], S[10], -U[3], U[1]] == {{((I/2)*EL^2*GaugeXi[Q])/SW}}, 
 C[-S[30], S[10], -U[1], U[3]] == {{((I/2)*EL^2*GaugeXi[Q])/SW}}, 
 C[S[30], S[10], -U[1], U[4]] == {{((I/2)*EL^2*GaugeXi[Q])/SW}}, 
 C[-S[30], S[20], -U[4], U[2]] == {{(EL^2*GaugeXi[Q])/(2*CW)}}, 
 C[S[30], S[20], -U[3], U[2]] == {{-(EL^2*GaugeXi[Q])/(2*CW)}}, 
 C[-S[30], S[20], -U[2], U[3]] == {{(EL^2*GaugeXi[Q])/(2*CW)}}, 
 C[S[30], S[20], -U[2], U[4]] == {{-(EL^2*GaugeXi[Q])/(2*CW)}}, 
 C[-S[30], S[20], -U[4], U[1]] == {{(EL^2*GaugeXi[Q])/(2*SW)}}, 
 C[S[30], S[20], -U[3], U[1]] == {{-(EL^2*GaugeXi[Q])/(2*SW)}}, 
 C[-S[30], S[20], -U[1], U[3]] == {{(EL^2*GaugeXi[Q])/(2*SW)}}, 
 C[S[30], S[20], -U[1], U[4]] == {{-(EL^2*GaugeXi[Q])/(2*SW)}}}


GaugeXi[ V[1 | 2 | 3] ] = GaugeXi[Q];
GaugeXi[ V[10 | 30 | 30] ] = GaugeXi[bg];
GaugeXi[ S[1 | 10] ] = 1;
GaugeXi[ S[2 | 3] ] = GaugeXi[Q];   
GaugeXi[ S[20 | 30] ] = GaugeXi[bg];
GaugeXi[ U[1 | 2 | 3 | 4] ] = GaugeXi[Q]


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


M$LastModelRules = {}


(* some short-hands for excluding classes of particles *)

QEDOnly = ExcludeParticles -> {F[1], V[2], V[3], S, SV, U[2], U[3], U[4]}

NoGeneration1 = ExcludeParticles -> F[_, {1}]

NoGeneration2 = ExcludeParticles -> F[_, {2}]

NoGeneration3 = ExcludeParticles -> F[_, {3}]

NoElectronHCoupling =
  ExcludeFieldPoints -> {
    FieldPoint[_][-F[2, {1}], F[2, {1}], S],
    FieldPoint[_][-F[2, {1}], F[1, {1}], S] }

NoLightFHCoupling =
  ExcludeFieldPoints -> {
    FieldPoint[_][-F[2], F[2], S],
    FieldPoint[_][-F[2], F[1], S],
    FieldPoint[_][-F[3, {1}], F[3, {1}], S],
    FieldPoint[_][-F[3, {2}], F[3, {2}], S],
    FieldPoint[_][-F[4], F[4], S],
    FieldPoint[_][-F[4], F[3, {1, ___}], S],
    FieldPoint[_][-F[4], F[3, {2, ___}], S] }

NoQuarkMixing =
  ExcludeFieldPoints -> {
    FieldPoint[_][-F[4, {1}], F[3, {2}], S[3 | 30]],
    FieldPoint[_][-F[4, {1}], F[3, {2}], V[3 | 30]],
    FieldPoint[_][-F[4, {1}], F[3, {3}], S[3 | 30]],
    FieldPoint[_][-F[4, {1}], F[3, {3}], V[3 | 30]],
    FieldPoint[_][-F[4, {2}], F[3, {1}], S[3 | 30]],
    FieldPoint[_][-F[4, {2}], F[3, {1}], V[3 | 30]],
    FieldPoint[_][-F[4, {2}], F[3, {3}], S[3 | 30]],
    FieldPoint[_][-F[4, {2}], F[3, {3}], V[3 | 30]],
    FieldPoint[_][-F[4, {3}], F[3, {1}], S[3 | 30]],
    FieldPoint[_][-F[4, {3}], F[3, {1}], V[3 | 30]],
    FieldPoint[_][-F[4, {3}], F[3, {2}], S[3 | 30]],
    FieldPoint[_][-F[4, {3}], F[3, {2}], V[3 | 30]] }


(* The following definitions of renormalization constants
   are for the on-shell renormalization of the Standard Model in
   the background-field formulation in the scheme of A. Denner,
   S. Dittmaier, and G. Weiglein, Nucl Phys B440 (1995) 95.

   The renormalization constants are not directly used by
   FeynArts, and hence do not restrict the generation of diagrams
   and amplitudes in any way. *)

Clear[RenConst]

RenConst[ dMf1[type_, j1_] ] := MassRC[F[type, {j1}]]

RenConst[ dZfL1[type_, j1_, j2_] ] :=
  FieldRC[F[type, {j1}], F[type, {j2}]][[1]]

RenConst[ dZfR1[type_, j1_, j2_] ] :=
  FieldRC[F[type, {j1}], F[type, {j2}]][[2]]

RenConst[ dMZsq1 ] := MassRC[V[20]]

RenConst[ dMWsq1 ] := MassRC[V[30]]

RenConst[ dMHsq1 ] := MassRC[S[10]]

RenConst[ dZAA1 ] := FieldRC[V[10]]

RenConst[ dZAZ1 ] := 2 CW/SW (dMWsq1/MW^2 - dMZsq1/MZ^2)

RenConst[ dZZA1 ] := 0

RenConst[ dZZZ1 ] := dZAA1 - (CW^2 - SW^2)/(2 CW SW) dZAZ1

RenConst[ dZW1 ] := dZAA1 - CW/SW/2 dZAZ1

RenConst[ dZH1 ] := dZW1 + dMWsq1/MW2

RenConst[ dTH1 ] := TadpoleRC[S[10]]

RenConst[ dZe1 ] := -1/2 dZAA1

RenConst[ dWFZ1 ] := FieldRC[V[20]] - dZZZ1

RenConst[ dWFAZ1 ] := FieldRC[V[10], V[20]] - dZAZ1

RenConst[ dWFW1 ] := FieldRC[V[30]] - dZW1

