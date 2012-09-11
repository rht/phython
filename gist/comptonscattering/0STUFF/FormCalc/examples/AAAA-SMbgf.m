(*
	AAAA-SMbgf.m  
		generates the Fortran code for
		gamma gamma -> gamma gamma in the electroweak SM
		using the background-field method
		this file is part of FormCalc
		last modified 9 May 08 th

Reference: M. Böhm, R. Schuster, Z. Phys. C63 (1994) 219.

*)


Needs["FeynArts`"]

Needs["FormCalc`"]


time1 = SessionTime[]

CKM = IndexDelta


process = {V[10], V[10]} -> {V[10], V[10]}

name = "AAAA-SMbgf"

SetOptions[InsertFields, Model -> "SMbgf", GenericModel -> "Lorentzbgf"]


SetOptions[Paint, PaintLevel -> {Classes}, ColumnsXRows -> {4, 5}]

(* take the comments out if you want the diagrams painted
$PaintSE = MkDir[name <> ".diagrams"];
DoPaint[diags_, file_, opt___] := Paint[diags, opt,
  DisplayFunction -> (Export[ToFileName[$PaintSE, file <> ".ps"], #]&)]
*)


Print["Boxes"]

tops = CreateTopologies[1, 2 -> 2, BoxesOnly];
ins = InsertFields[tops, process];
DoPaint[ins, "box"];
box = CalcFeynAmp[CreateFeynAmp[ins]]


box = Abbreviate[box, 5,
  Preprocess -> OnSize[100, Simplify, 500, DenCollect]]

abbr = OptimizeAbbr[Abbr[]]

subexpr = OptimizeAbbr[Subexpr[]]

dir = SetupCodeDir[name <> ".fortran", Drivers -> "AAAA-SM.drivers"]

WriteSquaredME[{}, box, abbr, subexpr, dir]

WriteRenConst[{}, dir]


Print["time used: ", SessionTime[] - time1]

