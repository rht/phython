(*
	ModelMaker.m
		derives the Feynman rules from the Lagrangian
		and splices them into a template model file
		last modified 28 Sep 07 th
*)


(* symbols from the model files live in Global` *)

{ MetricTensor, FourVector, ScalarProduct,
  ChiralityProjector, DiracMatrix, DiracSlash }

Attributes[MetricTensor] = {Orderless}


BeginPackage["ModelMaker`", {"FeynArts`", "Global`"}]

FunctionalD::usage = "FunctionalD[expr, fields] takes the functional
derivative of expr with respect to fields."

FeynmanRules::usage = "FeynmanRules[lagr] returns the list of Feynman
rules derived from the Lagrangian lagr.  These rules still contain the
kinematic parts, i.e. the couplings have not yet been decomposed into a
kinematic and a coupling vector."

CouplingVector::usage = "CouplingVector[rul] extracts the coupling vector
of each Feynman rule in rul according to the corresponding kinematic
vector in the currently initialized generic model file.  The output is in
a form suitable for the classes model file."

WriteModelFile::usage = "WriteModelFile[rul, \"MODEL.mmod\"] writes the
Feynman rules rul to the model file MODEL.mod using the template model
file MODEL.mmod.  The location in MODEL.mmod where the Feynman rules are
to be inserted must be marked by \"<* M$CouplingMatrices *>\"."

NFields::usage = "NFields[term] returns the number of fields in term."


Begin["`Private`"]

FunctionalD[expr_, f___] :=
  FuncD[
    Expand[expr /. d_Dot :> Distribute[ExpandAll[d]], Field] /.
                   fi_Field^n_ :> Dot@@ Table[fi, {n}],
    f ] /. IndexDelta[_] -> 1


(* version 1: no second arg, expr is FuncD'd wrt all fields *)

FuncD[p_Plus] := Flatten[FuncD/@ List@@ p]

FuncD[c_] := (Message[FeynmanRules::nocontrib, c]; {}) /;
  FreeQ[c, Field]

FuncD[expr_] :=
Block[ {vert, perm},
  If[ $GenericModel === "", InitializeModel[] ];
  vert = Sort[Cases[expr, _Field, Infinity],
    OrderedQ[ToClasses[First/@ {##}]]&];
  perm = FindVertex[ToGeneric[First/@ vert], Generic];
  If[ perm === $Failed,
    Message[FeynmanRules::nogeneric, t];
    Return[0] ];
  vert = MapIndexed[IndexRename, vert[[perm]], {2}];
  {C@@ vert == FuncD[expr, vert]}
]

IndexRename[s_. (f:P$Generic)[t_, i_List], {n_, 1}] :=
Block[ {char = FromCharacterCode[96 + n]},
  s f[t, Array[ToExpression[char <> ToString[#]]&, Length[i]]]
]

IndexRename[f_, {_, 1}] = f

IndexRename[ki_List, {n_, _}] := Through[Take[KIs, Length[ki]][n]]

IndexRename[_, {n_, _}] := Mom[n]


(* version 2: "ordinary" FuncD wrt second arg *)

FuncD[p_Plus, f_] := FuncD[#, f]&/@ p

FuncD[expr_, f_List] := 0 /; Count[expr, _Field, Infinity] =!= Length[f]

FuncD[expr_, f_List] := ReduceDeltas[Fold[FuncD, expr, f]]

FuncD[expr_, f_Field] :=
Block[ {NCSign, s, times, texpr},
  If[ !FreeQ[f, P$NonCommuting],
    s = -1;
    NCSign[P$NonCommuting] := s = -s ];
  _NCSign = 1;
  texpr = times[expr];
  Plus@@ (ReduceDeltas[IndexD[Extract[texpr, #], f] Delete[texpr, #]]&)/@
    Position[texpr, _Field] /. times -> Times
]


IndexD[ _[s_. (f:P$Generic)[t_, i1___], mom1___, ki1___List],
        _[s_. (f:P$Generic)[t_, i2___], mom2___, ki2___List] ] :=
  NCSign[f] ClassDelta[i1][i2] MomRename[mom1, mom2] KinDelta[ki1][ki2]

IndexD[__] = 0


ClassDelta[][] = KinDelta[][] = 1

ClassDelta[i1_][i2_] := Inner[ClassRename, i1, i2, Times] 

KinDelta[ki1_][ki2_] := Inner[MetricTensor, ki1, ki2, Times]

MomRename[] = MomRename[_] =
ClassRename[] = ClassRename[_] = 1


ReduceDeltas[MomRename[i_, j_] r_.] :=
  ReduceDeltas[r /. i -> j]

ReduceDeltas[ClassRename[i_, j_] IndexDelta[i_] r_.] :=
  IndexDelta[i, j] ReduceDeltas[r]

ReduceDeltas[ClassRename[i_, j_] r_.] :=
  IndexDelta[j] ReduceDeltas[r /. i -> j]

ReduceDeltas[other_] = other


MetricTensor/: MetricTensor[mu_, nu_] p_Plus :=
  (MetricTensor[mu, nu] #)&/@ p /; !FreeQ[p, nu]

MetricTensor/: MetricTensor[mu_, nu_] h_[a___, nu_, b___] :=
  h[a, mu, b]


NFields[p_Plus] := Block[{n = NFields/@ List@@ p}, n[[1]] /; SameQ@@ n]

NFields[_Plus] = Indeterminate

NFields[term_] := Exponent[term /. _Field -> Field /. Dot -> Times, Field]


FeynmanRules::nocontrib = "Warning: The term `` does not contain any
fields and thus does not contribute to the Feynman rules."

FeynmanRules::nogeneric = "No field point corresponding to the vertex
`` was found in the generic model file."

FeynmanRules[lagr_] :=
Block[ {cpl},
  cpl = Apply[AddTo, Hold@@ FunctionalD[lagr], 1];
  cpl = Block[ {C},
    _C = 0;
    ReleaseHold[cpl];
    DownValues[C] ];
  Flatten[ Apply[ToRule, cpl, 1] ]
]


ToRule[_, 0] = {}

ToRule[_[c_], rhs_] := First/@ c == I rhs


CouplingVector::nomatch = "`` in `` does not match any component of the
kinematic vector."

Attributes[CouplingVector] = {Listable}

CouplingVector[rul:_ == _List] = rul

CouplingVector[fi_ == cpl_] :=
Block[ {kc, kv, cv},
  kc = Kinalyze[Expand[cpl]] /. Kin -> KinExpand;
  kv = KinematicVector@@ ToGeneric[fi];
  cv = UnDot[kc, kv] //Flatten;
  If[ cv[[-1, 1]] =!= 0,
	(* try momentum conservation *)
    kc = kc /. Kin[Mom[n_][mu_]] :> -Plus@@ (Kin[Mom[#][mu]]&)/@ 
      Drop[Range[Length[fi]], {n}];
    cv = UnDot[kc, kv] //Flatten;
    If[ cv[[-1, 1]] =!= 0,
      Message[CouplingVector::nomatch, cv[[-1, 1]], fi] ]
  ];
  fi == List/@ Drop[cv, -1]
]

CouplingVector[other_] := CouplingVector[FeynmanRules[other]]


fermobjs = NonCommutative | ChiralityProjector | DiracMatrix | DiracSlash

bosobjs = _Mom | MetricTensor | FourVector | ScalarProduct

kinobjs = Join[bosobjs, fermobjs]

KinCat[f:fermobjs[__]] := {f, {}, {}}

KinCat[b:bosobjs[__]] := {{}, b, {}}

KinCat[other_] := {{}, {}, other}


Kinalyze[expr_] := expr Kin[1] /; FreeQ[expr, kinobjs]

Kinalyze[p_Plus] := Kinalyze/@ p

Kinalyze[t_Times] := Kinalyze@@ t

Kinalyze[t__] :=
Block[ {nc, kin, coeff},
  {nc, kin, coeff} = Flatten/@
    Transpose[KinCat/@ ({t} /. Dot -> NonCommutative)];
  Kin[If[Length[nc] === 0, 1, Flatten[NonCommutative@@ nc]] Times@@ kin] *
    Times@@ coeff
]


Kin[_[ga_DiracMatrix]] :=
  Kin[NonCommutative[ga, ChiralityProjector[+1]]] +
  Kin[NonCommutative[ga, ChiralityProjector[-1]]]


KinExpand[expr_] :=
  Distribute[Kin[
    Expand[expr /. FourVector[k_, mu_] :> (k /. m_Mom -> m[mu])]
  ]] /. Kin[n_?NumberQ x_] -> n Kin[x]

KinCoeff[c_, n_. k_Kin + r_.] :=
  {ReleaseHold[Coefficient[c, k]]/n, KinCoeff[c /. k -> 0, r]}

UnDot[c_, {k_, r___}] :=
Block[ {v = Flatten[ KinCoeff[c, KinExpand[k]] ]},
  If[ SameQ@@ Expand[Drop[v, -1]],
    {v[[1]], UnDot[v[[-1, 1]], {r}]},
    {0, UnDot[c, {r}]} ]
]


WriteModelFile[rulz_, template_] :=
Block[ {M$CouplingMatrices = rulz},
  Splice[template, FormatType -> InputForm]
]


End[]

EndPackage[]

