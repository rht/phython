(*
	eettA-SM.m
		generates the Fortran code for
		e^+ e^- -> t-bar t gamma in the electroweak SM
		this file is part of FormCalc
		last modified 9 May 08 th

Reference: W. Beenakker, S.C. van der Marck, and W. Hollik,
           Nucl. Phys. B365 (1991) 24.

*)


Needs["FeynArts`"]

Needs["FormCalc`"]


time1 = SessionTime[]

CKM = IndexDelta

Neglect[ME] = Neglect[ME2] = 0


process = {-F[2, {1}], F[2, {1}]} -> {-F[3, {3}], F[3, {3}], V[1]}

name = "eettA-SM"

SetOptions[InsertFields, Model -> "SM", Restrictions -> NoLightFHCoupling]


SetOptions[Paint, PaintLevel -> {Classes}, ColumnsXRows -> {4, 5}]

(* take the comments out if you want the diagrams painted
$PaintSE = MkDir[name <> ".diagrams"];
DoPaint[diags_, file_, opt___] := Paint[diags, opt,
  DisplayFunction -> (Export[ToFileName[$PaintSE, file <> ".ps"], #]&)]
*)


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

