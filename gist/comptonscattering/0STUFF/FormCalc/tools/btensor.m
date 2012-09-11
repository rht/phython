(*
	btensor.m
		explicit decompositions of the two-point
		tensor-coefficient functions
		this file is part of FormCalc
		last modified 3 Oct 09 th
*)


Clear[A0, A00, B0, B1, B00, B11, B001, B111, DB0, DB1, DB00, DB11]


A0[0] = 0

A00[m_] = 1/4 m (A0[m] + m/2)


B0[0, 0, 0] = BAD[B0]	(* divergent, must cancel *)

B0[0, m_, m_] := A0[m]/m - 1 /; m =!= 0

B0[m_, 0, m_] := A0[m]/m + 1 /; m =!= 0

B0[m_, m_, 0] := A0[m]/m + 1 /; m =!= 0

(* B0[m_, 0, 0] := A0[m]/m + 1 + I Pi /; m =!= 0 *)

B0[0, m1_, m2_] := (A0[m1] - A0[m2])/(m1 - m2)


DeltaB0[p_, m1_, m2_] := (B0[p, m1, m2] - B0[0, m1, m2])/p

DeltaB1[p_, m1_, m2_] := (B1[p, m1, m2] - B1[0, m1, m2])/p

DeltaB00[p_, m1_, m2_] := (B00[p, m1, m2] - B00[0, m1, m2])/p


B1[0, m_, m_] := -B0[0, m, m]/2

B1[0, m1_, m2_] := -B0[0, m1, m2]/2 -
  (m1^2 - m2^2 + 2 (m2 A0[m1] - m1 A0[m2]))/(4 (m1 - m2)^2)

B1a[0, m1_, m2_] := 1/2/(m1 - m2) (
  A0[m2] - m1 B0[0, m1, m2] - (m1 + m2)/2 )

B1b[0, m1_, m2_] := -1/2 (
  B0[0, m1, m2] + (m1 - m2) DB0[0, m1, m2] )

B1[p_, m1_, m2_] := 1/(2 p) (
  A0[m1] - A0[m2] - (p + m1 - m2) B0[p, m1, m2] )

B1a[p_, m1_, m2_] := -1/2 (
  B0[p, m1, m2] + (m1 - m2) DeltaB0[p, m1, m2] )


B00[0, m1_, m2_] := 1/4 (A0[m2] + m1 B0[0, m1, m2] + (m1 + m2)/2)

B00a[0, m1_, m2_] := (A00[m1] - A00[m2])/(m1 - m2)

B00[p_, m1_, m2_] :=
  -(p - 3 (m1 + m2))/18 + m1 B0[p, m1, m2]/3 +
  (A0[m2] + (p + m1 - m2) B1[p, m1, m2])/6


B11[0, m1_, m2_] := 1/18 - 2/3 B1[0, m1, m2] -
  2/3 (m1 - m2) DB1[0, m1, m2] - 1/3 m1 DB0[0, m1, m2]

B11a[p_, m1_, m2_] := 1/18 - 2/3 B1[p, m1, m2] -
  2/3 (m1 - m2) DeltaB1[p, m1, m2] - 1/3 m1 DeltaB0[p, m1, m2]

B11[p_, m1_, m2_] := 1/(3 p) (
  (p - 3 (m1 + m2))/6 + A0[m2] - m1 B0[p, m1, m2] -
  2 (p + m1 - m2) B1[p, m1, m2] )


B001[0, m_, m_] := -1/2 B00[p, m1, m2]

B001[0, m1_, m2_] := -1/(m1 - m2) (
  (m1 + m2)^2/36 + (m1 m2)/6 (B0[p, m1, m2] + 1/3) +
  (m1 - 2 m2)/3 B00[p, m1, m2] )

B001[p_, m1_, m2_] := 1/8 (
  2 m1 B1[p, m1, m2] - A0[m2] + (p + m1 - m2) (B11[p, m1, m2] + 1/6) -
  (m1 + m2)/2 )

B001a[p_, m1_, m2_] := 1/(2 p) (
  A00[m1] - A00[m2] - (p + m1 - m2)*B00[p, m1, m2] )

B001b[p_, m1_, m2_] :=
  -1/2 B00[p, m1, m2] - (m1 - m2)/2 DeltaB00[p, m1, m2]

B001c[p_, m1_, m2_] := 1/48 (
 -5 m1 - 7 m2 + 2 p - 6 m1 B0[p, m1, m2] +
 12 (m2 - p) B1[p, m1, m2] + 6 (m1 - m2 - 2 p) B11[p, m1, m2] )


B111[0, m_, m_] := -1/4 B0[0, m, m]

B111[0, m1_, m2_] :=
  3/16 + 1/4 A0[m2]/(m1 - m2) (m1/(m1 - m2) + 1) +
  1/2 m1/(m1 - m2) (m1/(m1 - m2) (B1[0, m1, m2] - 1/2) - 1/6)

B111[p_, m1_, m2_] := -1/(4 p) (
  A0[m2] + (p + m1 - m2) (3 B11[p, m1, m2] + 1/6) +
  2 m1 B1[p, m1, m2] - (m1 + m2)/2 )

B111a[p_, m1_, m2_] := -1/(2 p) (
  A0[m2] + (p + m1 - m2) B11[p, m1, m2] + 4 B001[p, m1, m2] )


Derivative[1, 0, 0][B0] = DB0;
Derivative[1, 0, 0][B1] = DB1;
Derivative[1, 0, 0][B00] = DB00;
Derivative[1, 0, 0][B11] = DB11

DB0[0, m_, m_] := 1/(6 m)

DB0[0, m1_, m2_] := 1/(m1 - m2)^2 (
  (m1 + m2)/2 - A0[m2] + m2 B0[0, m1, m2] )


DB1[0, m_, m_] = -1/(12 m)

DB1[0, m1_, m2_] =
  (2 m2/(m1 - m2) (B1[0, m1, m2] - B1[0, m2, m2]) - 1/3)/(m1 - m2)

DB1[p_, m1_, m2_] = D[B1[p, m1, m2], p]//Simplify


(*DB00[p_, m1_, m2_] = D[B00[p, m1, m2], p]//Simplify*)
DB00[p_, m1_, m2_] :=
  1/6 (2 m1 DB0[p, m1, m2] + B1[p, m1, m2] +
    (p + m1 - m2) DB1[p, m1, m2] - 1/3)


DB11[p_, m1_, m2_] = D[B11[p, m1, m2], p]//Simplify


B0[p_, m1_, m2_] := B0[p, m2, m1] /; !OrderedQ[{m1, m2}];
DB0[p_, m1_, m2_] := DB0[p, m2, m1] /; !OrderedQ[{m1, m2}]

