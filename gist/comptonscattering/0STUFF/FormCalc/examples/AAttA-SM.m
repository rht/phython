(*
	AAttA-SM.m
		generates the Fortran code for
		gamma gamma -> t-bar t gamma in the electroweak SM
		this file is part of FormCalc
		last modified 9 May 08 th

Reference: W. Walter, Diploma thesis, Würzburg 1997.

*)


Needs["FeynArts`"]

Needs["FormCalc`"]


time1 = SessionTime[]

CKM = IndexDelta

Neglect[ME] = Neglect[ME2] = 0


process = {V[1], V[1]} -> {-F[3, {3}], F[3, {3}], V[1]}

name = "AAttA-SM"

SetOptions[InsertFields, Model -> "SM", Restrictions -> NoLightFHCoupling]

SetOptions[CalcFeynAmp, Dimension -> 4]


SetOptions[Paint, PaintLevel -> {Classes}, ColumnsXRows -> {4, 5}]

(* take the comments out if you want the diagrams painted
$PaintSE = MkDir[name <> ".diagrams"];
DoPaint[diags_, file_, opt___] := Paint[diags, opt,
  DisplayFunction -> (Export[ToFileName[$PaintSE, file <> ".ps"], #]&)]
*)


Print["Born"]

tops = CreateTopologies[0, 2 -> 3];
ins = InsertFields[tops, process];
DoPaint[ins, "born"];
born = CalcFeynAmp[CreateFeynAmp[ins]]


col = ColourME[All, born]


abbr = OptimizeAbbr[Abbr[]]


dir = SetupCodeDir[name <> ".fortran", Drivers -> name <> ".drivers"]

WriteSquaredME[born, {}, col, abbr, dir]

WriteRenConst[{}, dir]


Print["time used: ", SessionTime[] - time1]

