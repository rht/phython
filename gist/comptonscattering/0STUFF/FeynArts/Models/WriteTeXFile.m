(*
	WriteTeXFile.m
		writes out the couplings of a FeynArts
		model file in TeX form
		last modified 17 Mar 09 th

	Usage:	WriteTeXFile["model.mod"]

NOTE: This program must currently be run in Mathematica versions
      4 or lower.
*)


BeginPackage["WriteTeXFile`", "FeynArts`"]

WriteTeXFile::usage = "WriteTeXFile[\"model.mod\"] writes the Feynman
rules in model.mod in TeX format to model.tex."

TeXFile::usage = "TeXFile is an option of WriteTeXFile.  It specifies
the output filename to use.  If set to Automatic, the name of the model
file with extension \".tex\" is used."

MaxLeaf::usage = "MaxLeaf is an option of WriteTeXFile.  It specifies
the maximum permissible leaf count per line."

Numbers::usage = "Numbers is an option of WriteTeXFile.  It specifies
whether to print the position of the coupling in the M$CouplingMatrices
array next to the coupling."

Sym::usage = "Sym[s, sub, sup] prints as symbol s with subscript sub and
superscript sup, where sub and sup are optional."

SymRules::usage = "SymRules specifies a list of rules to transform model
symbols."

WidthRules::usage = "WidthRules specifies a list of rules which are
supposed to replace symbols by objects with a leaf count that better
corresponds to their width.  This is used for computing line breaks."

ConjSym::usage = "ConjSym[f] specifies a conjugation symbol for field f.  
Choices are:\n
  ConjSym[f] = \"-\"\n
  -- f prints as f^-, anti-f prints as f^+,\n
  ConjSym[f] = \"\\\\nodagger\"\n
  -- f prints as f, anti-f prints as f^\\dagger,\n
  ConjSym[f] = Null\n
  -- f prints as f, anti-f prints as \\bar f."

Class::usage = "Class associates each field a class which is used to
group the couplings."


Begin["`Private`"]

template = (# <> "\n")&/@ ReadList[
  StringReplace[System`Private`FindFile[$Input], ".m" -> ".tex"],
  String ]


MakeList[Null] = {}

MakeList[{}] = {}

MakeList[s_Symbol] = {s}

MakeList[i__] :=
  Rest[Flatten[ {",", #}&/@ DeleteCases[Flatten[{i}], Null] ]]


Sym[sym_, sub_:{}, sup_:{}] := TeXSym[
  Sequence@@ Flatten[{sym}],
  Flatten[MakeList[sub]],
  Flatten[MakeList[sup]] ]

Format[TeXSym[sym__, sub_, sup_], TeXForm] :=
  Wrap[Identity][sym, 
    Wrap[Subscript][sub],
    Wrap[Superscript][sup /. "\\nodagger" -> Sequence[] /.
      {a__, ","} -> {a}]]

Wrap[_][{}] = Sequence[]

Wrap[h_][l_List] := h[SequenceForm@@ l]

Wrap[h_][s__] := h[SequenceForm[s]]

Wrap[h_][s_] := h[s]


TeXSym/: TeXSym[x__, sup_]^n_Integer?Positive :=
  TeXSym[x, Append[sup, n]]

TeXSym/: TeXSym[x__, sup_]^n_Integer :=
  1/TeXSym[x, Append[sup, -n]] /; n < -1

TeXSym/: Conjugate[TeXSym[x__, sup_]] := TeXSym[x, Append[sup, "*"]]

TeXSym/: Bar[TeXSym[x__, {r___, "-"}]] := TeXSym[x, {r, "+"}]

TeXSym/: Bar[TeXSym[x__, {r___, "+"}]] := TeXSym[x, {r, "-"}]

TeXSym/: Bar[TeXSym[x__, {r___, "\\nodagger"}]] :=
  TeXSym[x, {r, "\\dagger"}]

TeXSym/: Bar[TeXSym[x__, {r___, "\\dagger"}]] :=
  TeXSym[x, {r, "\\nodagger"}]


Unprotect[Times];
Format[n_Integer?Negative x_, TeXForm]:= SequenceForm["-", -n x];
Protect[Times]


Attributes[MakeFieldRules] = {Listable}

MakeFieldRules[f_ == desc_] :=
Block[ {Index = Identity, i = Indices /. desc, j, jn, Part},
  jn = j[[#]]&/@ Range[Length[i]];
  RuleDelayed@@ {
    If[ Length[jn] === 0, f, Append[f, j_] ],
    ToSym[ PropagatorLabel /. desc /. Thread[i -> jn],
      AddConj[SelfConjugate /. desc, ConjSym[f]] ]}
]

Attributes[ToSym] = {Listable}

ToSym[ComposedChar[s_, sub_, sup_, bar_], h_] :=
  h[bar <> " " <> s, sub, sup]

ToSym[ComposedChar[s_, sub_:Null, sup_:Null], h_] :=
  h[s, sub, sup]

ToSym[other_, h_] := h[other, Null, Null]

AddConj[True, _] = AddConj[False, Null] = Sym

AddConj[False, c_][s_, sub_, sup_] := Sym[s, sub, {sup, c}]


ToBar[-f_] := Bar[f]

ToBar[f_] = f

Format[Bar[f_], TeXForm] := SequenceForm["\\Bar{", f, "}"]


CoupVec[a_] = a

Format[CoupVec[a_, b__], TeXForm] :=
  TeXEnv["CoupVec", "\\\\\n\\\\[-3ex]\n", {a, b}]

Format[ZPlusB[a__], TeXForm] :=
  TeXEnv["PlusB", "\\\\\n", AddSigns[a]]

Format[ZPlusN[a__], TeXForm] :=
  TeXEnv["PlusN", "\\\\\n", AddSigns[a]]

Format[TeXEnv[name_, vspace_, a_], TeXForm] := SequenceForm[" ",
  "\n\\begin{" <> name <> "}\n",
  Sequence@@ Rest[Flatten[{vspace, #}&/@ a]],
  "\n\\end{" <> name <> "}"]


AddSigns[a_, b___] := Partition[
  { a,
    If[ MinusInFrontQ[#], {"\,-", -#}, {"\,+", #} ]&/@ {b}, 
    "" }//Flatten, 2 ]


(* splitting up long expressions *)

SmallEnough[expr__] := LeafCount[{expr} /. WidthRules] < maxleaf


Attributes[SplitLongPlus] = {HoldAll}

SplitLongPlus[x__] := Plus[x] /; SmallEnough[x]

SplitLongPlus[x__] :=
Block[ {cb, maxleaf = maxleaf - 5},
  ZPlusB@@ Flatten[ Operate[Coalesce, Plus[x]] ]
]


Coalesce[h_][a_, b_, r___] :=
  Coalesce[h][{a, b}, r] /; SmallEnough[a, b]

Coalesce[h_][a_, r___] := cb[ h@@ Flatten[{a}], Coalesce[h][r] ]

Coalesce[_][] = Sequence[]


ZPlusB/: lhs_ == ZPlusB[rhs__] := lhs == ZPlusN[rhs]


(* ordering inside a product *)

MinusInFrontQ[p_Plus] := MinusInFrontQ[ p[[1]] ]

MinusInFrontQ[_?Negative _.] = True

MinusInFrontQ[_] = False


TimesS[r__] :=
Block[ {t = Times[r]},
  If[ MinusInFrontQ[t],
    TimesO[-1, -t],
    TimesO[1, t] ]
]

TimesO[x_, r_. p_Plus] := If[ MinusInFrontQ[p],
  TimesO[-x ZPlusA@@ -p, r],
  TimesO[x ZPlusA@@ p, r] ]

TimesO[x_, r_. p_Plus^n_?Positive] := If[ MinusInFrontQ[p],
  TimesO[x (-1)^n (-p)^n, r],
  TimesO[x p^n, r] ]

TimesO[x_, r_. p_Plus^n_?Negative] :=
  TimesO[x (-1)^n, r (-p)^n] /; MinusInFrontQ[p]

TimesO[x_, r_. ZPlusB[a_, b__]] := If[ MinusInFrontQ[a],
  TimesO[-x ZPlusB@@ -{a, b}, r],
  TimesO[x ZPlusB[a, b], r] ]

TimesO[x_, r_. t_TeXSum] := TimesO[x t, r]

TimesO[x_, r_] := x HoldForm[r] /; Denominator[r] =!= 1

TimesO[x_, r_] := x r


Format[ZPlusA[p__], TeXForm] := Plus[p]


Format[TeXSum[a___, p_Plus], TeXForm] :=
  SequenceForm[a, "\\left(", p, "\\right)"]

TeXSum[s__, expr_] := -TeXSum[s, -expr] /; MinusInFrontQ[expr]

TeXSum[s__, ZPlusB[a_, b__]] := -TeXSum[s, ZPlusB@@ -{a, b}] /; 
  MinusInFrontQ[a]

Format[TeXSum[a__], TeXForm] := SequenceForm[a]


MakeSum[] = TeXSum

MakeSum[{v1_, r__}, a___, {v2_, r__}, b___] :=
  MakeSum[{{v1, v2}, r}, a, b]

MakeSum[{var_, from_:1, to_}, a___][b___, expr_] :=
  MakeSum[a][b,
    TeXSym["\\sum\\limits", {MakeList[var], "=", from}//Flatten, to],
    expr ]


(* sorting into classes of couplings *)

ToTeX[c_, lhs_, rhs_, n_] := class[c] = {class[c],
  If[ numbers, {"\\nbox{", ToString[n], "}$\n"}, "$\n" ],
  ToString[ToBar/@ lhs == CoupVec@@ rhs /.
    Plus -> SplitLongPlus //.
    SymRules /.
    FieldRules /.
    Conjugate[(t:Times)[a__]] :> Conjugate[Sym[{"(", a, ")"}]] /.
    Conjugate[x_] :> Conjugate[Sym[x]] /.
    IndexSum :> (MakeSum[##2][#1]&) /.
    Times -> TimesS, TeXForm],
  "\n$\\\\\n\\bigskip\n"}


Plural[p_Plus] := Plural/@ List@@ p

Plural[_[s_String]] := {" -- ", s}

Plural[n_ _[s_]] := {" -- ", ToString[n], " ", s} /;
  StringTake[s, -1] === "s"

Plural[n_ _[s_]] := {" -- ", ToString[n], " ", s, "s"}


AddCoup[lhs_ == rhs:{__List}, {n_}] :=
Block[ {cv = Transpose[rhs], name, h1, h2},
  If[ cto >= 0 && cto < Length[cv] && !MatchQ[cv = cv[[cto + 1]], {(0)...}],
    name = Class/@ lhs;
    h1 = ToString[Head[#]]&/@ List@@ name;
    h2 = Rest[Flatten[Plural[Plus@@ name]]];
    ToTeX["\\Class{" <> h1 <> "}{" <> h2 <> "}\n", lhs, cv, n] ]
]

AddCoup[c_, _] := Message[WriteTeXFile::badcoup, c]


WriteTeXFile::badcoup = "Warning: `` is not recognized as a coupling."

Options[WriteTeXFile] = {
  TeXFile -> Automatic,
  MaxLeaf -> 60,
  Numbers -> True,
  CTOrder -> 0 }

WriteTeXFile[model_, opt___Rule] :=
Block[ {texfile, maxleaf, numbers, cto,
mod, FieldRules, class, couplings, hh},

  {texfile, maxleaf, numbers, cto} = {TeXFile, MaxLeaf, Numbers, CTOrder} /.
    {opt} /. Options[WriteTeXFile];

  Check[ mod = LoadModel[model], Abort[] ];

  FieldRules = MakeFieldRules[M$ClassesDescription];

  _class = {};
  MapIndexed[AddCoup, M$CouplingMatrices];
  _class =.;

  couplings = {#[[1, 1, 1]], #[[2]]}&/@ DownValues[class];

  If[ texfile === Automatic, texfile = mod <> ".tex" ];
  If[ cto > 0,
    mod = mod <> " (" <> ToString[cto] <> "-loop counter terms)"];

  hh = OpenWrite[texfile];
  WriteString[hh,
    StringReplace[
      StringJoin[template /. "COUPLINGS\n" -> couplings],
      {"MODEL" -> StringJoin[mod], "\n \n" -> "\n"} ]
  ];
  Close[hh]
]

End[]

EndPackage[]


(* here come the model-dependent things *)

Class[_. S[11|12, _]] = S["Slepton"];
Class[_. S[13|14, _]] = S["Squark"];
Class[_. S[__]] = S["Higgs"];
Class[_. SV[_]] = SV["Scalar--Vector"];
Class[_. F[1|2, _]] = F["Lepton"];
Class[_. F[3|4, _]] = F["Quark"];
Class[_. F[11, _]] = F["Neutralino"];
Class[_. F[12, _]] = F["Chargino"];
Class[_. F[15, _]] = F["Gluino"];
Class[_. V[5, _]] = V["Gluon"];
Class[_. V[_]] = V["Gauge Boson"];
Class[_. U[__]] = U["Ghost"]

(* we want i and coupling constants to print in front,
   thus the silly name "AAA" *)

Format[AAA[x_], TeXForm] := x

AAA/: AAA[x_]^n_ := AAA[x^n]

ferm[1] = "\\nu";
ferm[2] = "e";
ferm[3] = "u";
ferm[4] = "d"

sferm[i_] := "\\tilde " <> ferm[i]

ConjSym[F[1|2|3|4]] = ConjSym[_U] = Null

ConjSym[S[11|12|13|14]] = "\\nodagger"

ConjSym[_] = "-"

SymRules = {
  $HKSign -> 1 (* "(-)" *),
  Complex[0, im_] -> im AAA["\\i"],
  EL :> AAA[Sym["e"]],
  GS :> AAA[Sym["g", "s"]],
  SW :> Sym["s", "W"],
  CW :> Sym["c", "W"],
  SA :> Sym["s", "\\alpha"],
  CA :> Sym["c", "\\alpha"],
  SB :> Sym["s", "\\beta"],
  CB :> Sym["c", "\\beta"],
  TB :> Sym["t", "\\beta"],
  SAB :> Sym["s", "\\alpha+\\beta"],
  SBA :> Sym["s", "\\beta-\\alpha"],
  CAB :> Sym["c", "\\alpha+\\beta"],
  CBA :> Sym["c", "\\beta-\\alpha"],
  S2A :> Sym["s", "2\\alpha"],
  S2B :> Sym["s", "2\\beta"],
  C2A :> Sym["c", "2\\alpha"],
  C2B :> Sym["c", "2\\beta"],
  MUE :> Sym["\\mu"],
  Lambda5 :> Sym["\\lambda", "5"],
  Yuk1 :> Sym["Y", "1"],
  Yuk2 :> Sym["Y", "2"],
  Yuk3 :> Sym["Y", "3"],
  MW :> Sym["M", "W"],
  MZ :> Sym["M", "Z"],
  MH :> Sym["M", "H"],
  Mh0 :> Sym["M", S[1]],
  MHH :> Sym["M", S[2]],
  MA0 :> Sym["M", S[3]],
  MHp :> Sym["M", S[5]],
  MLE[j_] :> Sym["m", F[2, {j}]],
  MQU[j_] :> Sym["m", F[3, {j}]],
  MQD[j_] :> Sym["m", F[4, {j}]],
  Mass[f_] :> Sym["m", f],
  UHiggs[h__] :> Sym["U", {h}, "H"],
  ZHiggs[h__] :> Sym["Z", {h}, "H"],
  ZNeu[n__] :> Sym["Z", {n}],
  VCha[c__] :> Sym["V", {c}],
  UCha[c__] :> Sym["U", {c}],
  USf[t_, g_][s__] :> Sym["U", {s}, {sferm[t], g}],
  UASf[t_][as__] :> Sym["R", {as}, sferm[t]],
  Af[t_, i__] :> Sym["A", {i}, ferm[t]],
  CKM[g__] :> Sym["\\Mvariable{CKM}", {g}],
  SUNTSum[c1_, c2_, c3_, c4_] :>
    Sym[{SUNT[x, c1, c2], SUNT[x, c3, c4]}],
  SUNT[g_, c1_, c2_] :> Sym["T", {c1, c2}, g],
  SUNT[g__, c1_, c2_] :>
    Sym[{"(", Sym["T", Null, #]&/@ {g}, ")"}, {c1, c2}],
  SUNF[g1_, g2_, g3_, g4_] :>
    Sym[{SUNF[g1, g2, "x"], SUNF["x", g3, g4]}],
  SUNF[g1_, g2_, g3_] :> Sym["f", Null, {g1, g2, g3}],
  d_IndexDelta :> 1 /; !FreeQ[d, o1|o2|o3|o4],
  IndexDelta[n_Integer, i_] :> Sym["\\delta", {i, n}],
  IndexDelta[i__] :> Sym["\\delta", {i}],
  GaugeXi[v_] :> Sym["\\xi", v],
  dZe1 :> Sym["\\delta Z", "e"],
  dMHsq1 :> Sym["\\delta M", "H", "2"],
  dMWsq1 :> Sym["\\delta M", "W", "2"],
  dMZsq1 :> Sym["\\delta M", "Z", "2"],
  dMf1[t_, j_] :> Sym["\\delta m", j, ferm[t]],
  dMf1[t_, j_] :> Sym["\\delta m", j, ferm[t]],
  dSW1 :> Sym["\\delta s", "W"],
  dCW1 :> Sym["\\delta c", "W"],
  dZH1 :> Sym["\\delta Z", "H"],
  dZW1 :> Sym["\\delta Z", "W"],
  dZAA1 :> Sym["\\delta Z", "\\gamma\\gamma"],
  dZZA1 :> Sym["\\delta Z", "Z\\gamma"],
  dZAZ1 :> Sym["\\delta Z", "\\gamma Z"],
  dZZZ1 :> Sym["\\delta Z", "ZZ"],
  dUW1 :> Sym["\\delta U", "W"],
  dUAA1 :> Sym["\\delta U", "\\gamma\\gamma"],
  dUZA1 :> Sym["\\delta U", "Z\\gamma"],
  dUAZ1 :> Sym["\\delta U", "\\gamma Z"],
  dUZZ1 :> Sym["\\delta U", "ZZ"],
  dZG01 :> Sym["\\delta Z", Sym["G", Null, "0"]],
  dZGp1 :> Sym["\\delta Z", "G"],
  dZfL1[t_, j1_, j2_] :> Sym["\\delta Z", {j1, j2}, {ferm[t], "L"}],
  dZfR1[t_, j1_, j2_] :> Sym["\\delta Z", {j1, j2}, {ferm[t], "R"}],
  dCKM1[j1_, j2_] :> Sym["\\delta\\Mvariable{CKM}", {j1, j2}],
  dTad1 :> Sym["\\delta T"]
}

WidthRules = {
  s:CBA|SBA|CAB|SAB|C2A|S2A|C2B|S2B -> s[1],  (* count as 2 chars *)
  m:Mh0|MHH|MA0|MHp -> m[1],  (* ditto *)
  _Mass -> m[1],
  h_[__][i__] -> h[i, 1],  (* upper indices of e.g. USf don't count *)
  (d:dZfL1|dZfR1)[t_, j1_, j2_] -> d[j1, j2],
  -1 -> +1,  (* a - b counts same as a + b *)
  Conjugate -> Identity,  (* USf^* counts as USf *)
  Power -> (If[IntegerQ[#2], #1, #1^#2]&)  (* M_W^2 counts as M_W *)
}

Null

