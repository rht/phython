(*
	HMix.mod
		Add-on model file which adds two new particles:
		S[0,  {h}] = Sum[UHiggs[h,i] S[i], {i, 3}]
		S[10, {h}] = Sum[ZHiggs[h,i] S[i], {i, 3}]
		this file is part of FeynArts
		last modified 8 Sep 08 th
*)


IndexRange[ Index[Higgs] ] = Range[3]

M$ClassesDescription = Flatten[{
  M$ClassesDescription,
  S[0] == {
	SelfConjugate -> True,
	Indices -> {Index[Higgs]},
	Mass -> MHiggs,
	PropagatorLabel -> ComposedChar["H", Index[Higgs]],
	PropagatorType -> ScalarDash,
	PropagatorArrow -> None },
  S[10] == {
	SelfConjugate -> True,
	Indices -> {Index[Higgs]},
	InsertOnly -> {External},
	Mass -> MHiggs,
	PropagatorLabel -> ComposedChar["H", Index[Higgs], Null, "\\hat"],
	PropagatorType -> ScalarDash,
	PropagatorArrow -> None }
}]


Block[ {NewCoup, UZPerm, coup, oldcoup, newcoup},

NewCoup[c_ == rhs_] :=
Block[ {p = Position[c, S[1|2|3]], sign},
  sign = If[ ToGeneric[c] === C[S, S, V], Signature, 1 & ];
  Block[ {newi = Take[{h1, h2, h3, h4}, Length[p]]},
    (coup[#1] += #2 rhs)&[
      ReplacePart[c, S[0, {#}]&/@ newi, p, Array[List, Length[p]]],
      Plus@@ (sign[#] Times@@ MapThread[UHiggs, {newi, #}]&)/@
        Permutations[Apply[c[[#, 1]]&, p, 1]]
    ];
    {}
  ] /; Length[p] > 0
];

NewCoup[other_] = other;


UZPerm[_[_[c_]], rhs_] := FoldList[
  ReplaceAll,
  c == rhs (* (Simplify[rhs] /. CB TB -> SB) *),
  Cases[c, S[0, {h_}] -> {UHiggs[h, i_] -> ZHiggs[h, i],
                          S[0, {h}] -> S[10, {h}]}] ];

_coup = 0;
oldcoup = NewCoup/@ M$CouplingMatrices;
_coup =.;
newcoup = Apply[UZPerm, DownValues[coup], 1];


If[ TrueQ[$JustNewCouplings], oldcoup = {} ];

M$CouplingMatrices = Flatten[{oldcoup, newcoup}];

]


DownValues[RenConst] = DownValues[RenConst] /. S[h:1|2|3] -> S[0, {h}]

