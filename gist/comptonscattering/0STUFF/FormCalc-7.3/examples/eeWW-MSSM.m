(*
	eeWW-MSSM.m
		generates the Fortran code for
		e^+ e^- -> W^+ W^- in the MSSM
		this file is part of FormCalc
		last modified 9 May 08 th

Reference: T. Hahn, Nucl. Phys. B609 (2001) 344 [hep-ph/0007062].
*)


Needs["FeynArts`"]

Needs["FormCalc`"]


time1 = SessionTime[]

CKM = IndexDelta

Neglect[ME] = Neglect[ME2] = 0


process = {-F[2, {1}], F[2, {1}]} -> {-V[3], V[3]}

name = "eeWW-MSSM"

SetOptions[InsertFields, Model -> "MSSM", Restrictions -> NoLightFHCoupling]


SetOptions[Paint, PaintLevel -> {Classes}, ColumnsXRows -> {4, 5}]

(* take the comments out if you want the diagrams painted
$PaintSE = MkDir[name <> ".diagrams"];
DoPaint[diags_, file_, opt___] := Paint[diags, opt,
  DisplayFunction -> (Export[ToFileName[$PaintSE, file <> ".ps"], #]&)]
*)


Print["Counter terms"]

tops = CreateCTTopologies[1, 2 -> 2,
  ExcludeTopologies -> {TadpoleCTs, WFCorrectionCTs}];
ins = InsertFields[tops, process, Model -> "SM"];
DoPaint[ins, "counter"];
counter = CreateFeynAmp[ins] /. SW -> -SW


Print["Born"]

tops = CreateTopologies[0, 2 -> 2];
ins = InsertFields[tops, process];
DoPaint[ins, "born"];
born = CalcFeynAmp[CreateFeynAmp[ins]]


Print["Self energies"]

tops = CreateTopologies[1, 2 -> 2, SelfEnergiesOnly];
ins = InsertFields[tops, process];
ins = DiagramSelect[ins, FreeQ[#, Field[5|6] -> S]&];
DoPaint[ins, "self"];
self = CalcFeynAmp[
  CreateFeynAmp[ins],
  Select[counter, DiagramType[#] == 2 &]]



Print["Vertices"]

tops = CreateTopologies[1, 2 -> 2, TrianglesOnly];
ins = InsertFields[tops, process];
ins = DiagramSelect[ins, FreeQ[#, Field[5] -> S]&];
DoPaint[ins, "vert"];
vert = CalcFeynAmp[
  CreateFeynAmp[ins],
  Select[counter, DiagramType[#] == 1 &]]


Print["Boxes"]

tops = CreateTopologies[1, 2 -> 2, BoxesOnly];
ins = InsertFields[tops, process];
DoPaint[ins, "box"];
box = CalcFeynAmp[
  CreateFeynAmp[ins],
  Select[counter, DiagramType[#] == 0 &]]


amps = {born, self, vert, box}

{born, self, vert, box} = Abbreviate[amps, 6,
  Preprocess -> OnSize[100, Simplify, 500, DenCollect]]

abbr = OptimizeAbbr[Abbr[]]

subexpr = OptimizeAbbr[Subexpr[]]

dir = SetupCodeDir[name <> ".fortran", Drivers -> name <> ".drivers"]

WriteSquaredME[born, {self, vert, box}, abbr, subexpr, dir]

WriteRenConst[amps, dir]


Print["time used: ", SessionTime[] - time1]

