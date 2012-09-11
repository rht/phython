(*
	Analytic.m
		Translation of InsertFields output into
		analytic expressions
		last modified 27 Mar 12 th
*)

Begin["`Analytic`"]

Options[ CreateFeynAmp ] = {
  AmplitudeLevel -> InsertionLevel,	(* i.e. taken from InsertFields *)
  GaugeRules -> {_GaugeXi -> 1},
  PreFactor -> -I (2 Pi)^(-4 LoopNumber),
  Truncated -> False,
  MomentumConservation -> True,
  GraphInfoFunction -> (1 &)
}

(* for D dimensions use
   PreFactor -> -I (Global`Mu^(4 - $D)/(2 Pi)^$D)^LoopNumber *)

CreateFeynAmp::nolevel =
"Warning: Level `1` is not contained in this insertion."

CreateFeynAmp::mtrace =
"Different MatrixTraceFactors inside one loop.  \
Involved fields are `1`.  \
Please check the classes model and try again."

CreateFeynAmp::noprop =
"Cannot resolve propagator of field `1`."

CreateFeynAmp::ambig =
"Warning: `1` contains more than two noncommuting fields, hence the \
application of the flipping rules is not unambiguous and may give \
wrong results."

CreateFeynAmp::nocoupl =
"Cannot resolve coupling `1` for kinematical object `2`."

CreateFeynAmp::counter =
"Counter-term order `2` is not defined in coupling of `1`."


(* CreateFeynAmp invokes a hierarchy of functions:
     CreateFeynAmp[TopologyList]	select levels, init model
       CreateAmpTop[Topology]		add momenta, make fermion chains
         CreateAmpGraph[FeynmanGraph]	add indices, make generic expr
           CreateAmpIns[FeynmanGraph -> Ins] add GM replacement rules *)

CreateFeynAmp[ top:(P$Topology -> _), opt___Rule ] :=
  CreateFeynAmp[ TopologyList[][top], opt ]

CreateFeynAmp[ TopologyList[tops__], opt___Rule ] :=
  CreateFeynAmp[ TopologyList[][tops], opt ]

CreateFeynAmp[ tops:TopologyList[info___][___], options___Rule ] :=
Block[ {alevel, pref, next, gaugeru, truncru, momcons, graphinfo, toplist,
amps, head, topnr = 1, opt = ActualOptions[CreateFeynAmp, options]},

  If[ (alevel = ResolveLevel[AmplitudeLevel /. opt /. {info} /.
        Options[InsertFields]]) === $Failed,
    Return[$Failed] ];

  If[ InitializeModel[ Model /. {info} /. Options[InsertFields],
    GenericModel -> (GenericModel /. {info} /. Options[InsertFields]),
    Reinitialize -> False ] === $Failed, Return[$Failed] ];

  FAPrint[2, ""];
  FAPrint[2, "creating amplitudes at level(s) ", alevel];

  next = Plus@@ Length/@ (Process /. {info});
  pref = PreFactor /. opt;
  gaugeru = GaugeRules /. opt;
  truncru = If[ TrueQ[Truncated /. opt], M$TruncationRules, {} ];
  momcons = TrueQ[MomentumConservation /. opt];
  graphinfo = GraphInfoFunction /. opt;
  toplist = TopologyList[info];

  amps = PickLevel[alevel][tops];
  Scan[ If[FreeQ[amps, #], Message[CreateFeynAmp::nolevel, #]]&, alevel ];

  amps = CreateAmpTop/@ ( amps //.
    (_ -> Insertions[_][]) :> Seq[] /.
    (Field[i_] -> fi_?AtomQ) -> (Field[i] -> fi[Index[Generic, i]]) );
  FAPrint[1, "in total: ",
    Statistics[{Insertions[Generic]@@ amps}, alevel, " amplitude"]];

  head = FeynAmpList[info] /.
    (Process -> iorule_) :> (Process ->
      MapIndexed[{#1, iomom@@ #2, TheMass[#1, External]}&, iorule, {2}]) /.
    (InsertionLevel -> _) :> (AmplitudeLevel -> alevel);

  amps = head@@ amps /. _MTF -> 1 /. gaugeru //. M$LastModelRules;

  If[ Length[alevel] === 1, PickLevel[ alevel[[1]] ][amps], amps ]
]


iomom[ 1, n_ ] = FourMomentum[Incoming, n]

iomom[ 2, n_ ] = FourMomentum[Outgoing, n]


CreateAmpTop[ P$Topology ] = Sequence[]

CreateAmpTop[ top:P$Topology -> ins_ ] :=
Block[ {momtop, imom, oldmom, amp, c, toppref, mtf, mc = 0, gennr = 0},

	(* append momenta and enforce momentum conservation for
	   every vertex.  For economical reasons, external momentum
	   conservation (i.e. elimination of one momentum) is not
	   carried out. *)
  c[_] = 0;
  momtop = AppendMomentum/@ top;
  If[ momcons,
	(* since we don't touch the external momenta, it's
	   important to go through the vertices last to first
	   because the first are always those that connect the
	   external particles *)
    momtop = Catch[
      Fold[
        MomConservation,
        momtop,
        Reverse[RemoveDups[
          Cases[top, Vertex[n__][_] /; {n} =!= {1}, {2}]] ]] ]
  ];
	(* renumber the internal momenta *)
  oldmom = Union[ Cases[momtop, FourMomentum[_ZZZ, _], Infinity] ];
  imom = RenumberMom@@@ oldmom;
  momtop = momtop /. Thread[oldmom -> imom];

  toppref = pref /. LoopNumber :> Genus[top];

  amp = Sequence@@ (CreateAmpGraph[momtop, #]&)/@ ins;
  FAPrint[2, "> Top. ", topnr++, ": ",
    Statistics[{Insertions[Generic][amp]}, alevel, " amplitude"]];
  amp
]


RemoveDups[ li_ ] :=
Block[ {f},
  f[x_] := (f[x] = Sequence[]; x);
  f/@ li
]


(* loop number using Euler's relation: *)

Genus[ top_ ] := 
Block[ {c, vn = {}},
  c[n_] := (AppendTo[vn, n]; 0);
  ++c[ #[[0, 1]] ]&/@ Union[Cases[top, Vertex[__][_], {2}]];
  (Plus@@ ((# - 2) c[#] &)/@ vn)/2 + 1
]


FourMomentum[ type_, n_Integer mom_ ] :=
  -FourMomentum[type, -n mom] /; n < 0


	(* ZZZ[priority] gives a ranking for eliminating momenta:
	   of the sorted list of momenta at a vertex the last is
	   eliminated, hence ZZZ[4] is eliminated before ZZZ[3] etc.
	   Since momenta on tree propagators (Propagator[Internal])
	   can always be expressed by the external momenta, they
	   have the highest priority. *)

	(* in case a momentum is given from outside: *)
AppendMomentum[ Propagator[type_][from_, to_, fi_, mom_] ] :=
  Propagator[type][ from, to, fi,
    FourMomentum[type /. {Loop[_] -> ZZZ[1], Internal -> ZZZ[3]},
      mom /. FourMomentum[_, t_] -> t] ]

AppendMomentum[ pr:Propagator[Outgoing][from_, __] ] :=
  Append[ pr,
    If[from[[0, 1]] === 1, -1, 1] FourMomentum[Outgoing, ++c[Outgoing]] ]

AppendMomentum[ pr:Propagator[Incoming | External][from_, __] ] :=
  Append[ pr,
    If[from[[0, 1]] === 1, 1, -1] FourMomentum[Incoming, ++c[Incoming]] ]

AppendMomentum[ pr:Propagator[Loop[_]][__] ] :=
  Append[ pr, FourMomentum[ZZZ[2], ++c[Internal]] ]

AppendMomentum[ pr:Propagator[Internal][__] ] :=
  Append[ pr, FourMomentum[ZZZ[4], ++c[Internal]] ]


RenumberMom[ _, _Integer ] := FourMomentum[Internal, ++mc]

RenumberMom[ _, id_ ] = FourMomentum[Internal, id]


MomConservation[ top_, vert_ ] := Throw[top] /; FreeQ[top, ZZZ]

MomConservation[ top_, vert_ ] :=
Block[ {eq},
  eq = Plus@@ (IncomingMomentum[vert, #]&)/@ top;
  If[ eq === 0 || (Head[eq] === Plus && FreeQ[eq, ZZZ]), top,
    top /. If[ Head[eq] =!= Plus, eq -> 0,
      Solve[ eq == 0,
        Sort[Cases[{eq}, _FourMomentum, Infinity]][[-1]] ][[1]] ]
  ]
]

IncomingMomentum[ v_, _[v_, v_, ___] ] = 0

IncomingMomentum[ v_, _[v_, _, ___, m_] ] = -m

IncomingMomentum[ v_, _[_, v_, ___, m_] ] = m

IncomingMomentum[ __ ] = 0


app[ fi_ ] = fi


CreateAmpGraph[ top_, gr:FeynmanGraph[s_, ___][__] -> ins_ ] :=
Block[ {amp, gm, rawgm, orig, anti},
	(* must save Field[n] information to be able to subsequently
	   apply the insertion rules of deeper levels *)
  amp = CreateAmpGraph[ top,
    gr /. (n_ -> x_. fi_[i__]) :> (n -> x fi[i, orig[x, n]]) ];
  gm = Append[
    Union[ Cases[amp[[3]], P$InsertionObjects, Infinity, Heads -> True] ],
    RelativeCF ];
  rawgm = gm /.
    s1_. _[__, orig[s2_, fi_], k___] :>
      app[ If[s1 === s2, fi, anti[fi]], k ];
  Append[amp, gm -> (CreateAmpIns[top, rawgm, s mtf, #]&)/@ ins] /.
    orig[__] :> Seq[]
]

FieldNumber[ fi_ ] := Sequence@@ Cases[fi, Field[n_] -> n, Infinity, 1] /;
  !FreeQ[fi, orig]


(* Create the basic amplitude *)

CreateAmpGraph[ top_, FeynmanGraph[s_, ___][ru__] ] :=
Block[ {c, res, props, vert, faden, prden = {},
scalars = {RelativeCF, toppref, 1/s}},

  c[_] = 0;
  res = AddKinematicIndices/@ (List@@ top /. {ru});
  mtf = 1;
  If[ $FermionLines, res = MakeFermionChains[res] ];

	(* props contains the propagators not involved in gmcs *)
  props = Cases[res, Propagator[_][__]];
  vert = Vertices[props];

	(* insert the vertices in fermion chains first.  Note that
	   MidVertex and ToChain modify vert *)
  res = res /. ch:_dot | _tr :> ResolveChain[ch];

	(* now the remaining vertices *)
  vert = ResolveGeneric/@ vert;

	(* TakeNC does the multiplication business.  It also updates
	   scalars and prden *)
  res = Join[vert, res /. Propagator -> ResolvePropagator /. gaugeru] /.
    PV -> TakeNC;

  FeynAmp[
    GraphID[Topology == topnr, Generic == ++gennr],
    Integral@@ imom,
    Times@@ DeleteCases[Flatten[scalars], 1] *
      LoopPD[Expand[Times@@ Flatten[prden], PropagatorDenominator]] *
      Times@@ res /.
      truncru /. Mass -> TheMass ] //. M$LastGenericRules
]


AddKinematicIndices[
  Propagator[type_][vert__, s_. fi_[ind___], mom_] ] :=
Block[ {ki = KinematicIndices[fi], kin},
  If[ Length[ki] === 0, kin = Seq[],
    kin = If[ FreeQ[{vert}, Vertex[1]],
      Rule@@ Transpose[{Index[#, ++c[#]], Index[#, ++c[#]]}&/@ ki],
    (* else *)
      Index[#, ++c[#]]&/@ ki
    ] ];
  Propagator[type][vert, s fi[ind, mom, kin]]
]


(* Building fermion chains.
   Fermionic (noncommuting) objects need to be organized into chains
   if they don't carry explicit spinor indices (which is the default).
   This works unambiguously (i.e. correctly) only if fermion chains
   never "touch", i.e. if there are at most two fermions at each vertex.
   This is always the case in renormalizable theories.  Effective theories
   may however contain 4- or more-fermion vertices, typically as a result
   of integrating out heavy bosons, such as in the Fermi model.  In such a
   case one must set $FermionLines = False and give the fermion fields an
   explicit spinor index with which it is possible (outside of FeynArts)
   to find the correct ordering of the fermionic objects. *)

ReverseProp[ pr_[from_, to_, part_] ] :=
  pr[to, from, AntiParticle[part]]

Attributes[ BuildChain ] = {Flat, Orderless}

BuildChain[ c1:_[___, _[_, v_, _]], c2:_[_[v_, __], ___] ] :=
  Join[c1, c2]

BuildChain[ c1:_[___, _[_, v_, _]], c2:_[___, _[_, v_, _]] ] :=
  Join[c1, Reverse[ReverseProp/@ c2]]

BuildChain[ c1:_[_[v_, __], ___], c2:_[_[v_, __], ___] ] :=
  Join[Reverse[ReverseProp/@ c1], c2]


Fixgmc[ c__ ] := tr[c] /; FreeQ[{c}, Vertex[1]]

	(* a Dirac fermion, and it's in the right place *)
Fixgmc[ c:_[__, -_], r___ ] := dot[c, r]

	(* assuming that the front end must be a Majorana fermion, then *)
Fixgmc[ r___, c:_[__, -_] ] := dot[r, c]

	(* in principle there is no convention how to order Majorana 
	   lines; whenever Truncated -> True is used, however, the 
	   ordering, i.e. the type of spinor, does become important,
	   e.g. when computing counter terms from the self energies *)
Fixgmc[ c1:_[__, _?SelfConjugate], r___, c2:_[__, _?SelfConjugate] ] :=
  dot[c1, r, c2] /;
  OrderedQ[{c2[[0, 1]][ c2[[2, 1]] ], c1[[0, 1]][ c1[[1, 1]] ]}]

Fixgmc[ c__ ] := Reverse[ ReverseProp/@ dot[c] ]


MakeFermionChains[ top_ ] := top /; FreeQ[top, P$NonCommuting]

MakeFermionChains[ top_ ] :=
Block[ {res, ext},
  res = Append[
    Select[top, FreeQ[#[[3]], P$NonCommuting]&],
    BuildChain@@ gmc/@ Select[top, !FreeQ[#[[3]], P$NonCommuting]&] /.
      gmc -> Fixgmc
  ] /. BuildChain -> Sequence;

	(* Since fermion chains are always traversed opposite to the
	   fermion flow, we need the sign of the permutation that gets
	   the list of external fermions into _descending_ order.
	   Since Signature gives the sign for _ascending_ order, we
	   need another (-1)^(Length[ext]/2).
	   (Actually, the factor is (-1)^(len (len - 1)/2), but since
	   the number of external fermions is always even, (-1)^(len/2)
	   gives the same result.) *)
  ext = Flatten[ Cases[res, d_dot :> Extract[d, {{1, 1}, {-1, 2}}, Leg]] ];
  AppendTo[ scalars,
    Signature[ext] (-1)^(Count[res, _tr] + Length[ext]/2) ];

  mtf = Times@@ Cases[res, t_tr :> MTF[ Union[Cases[t, Field[_], {-2}]] ]];

  res
]


Leg[ Vertex[1][n_] ] = n

Leg[ _ ] = {}


MTF[ {} ] = 1

MTF[ fi_List ] :=
Block[ {res},
  res = Union[MatrixTraceFactor/@ fi];
  If[ Length[res] === 1,
    res[[1]] /. i_Index :> MTFIndex[i, fi],
  (* else *)
    Message[CreateFeynAmp::mtrace, fi]; 1 ]
] /; FreeQ[fi, Field]


Attributes[ MTFIndex ] = {Listable}

MTFIndex[ i_, _. fi_[n_, ind_] ] :=
  i /. Thread[Indices[fi[n]] -> ind]


	(* useful for Truncated -> True *)
NonCommutative[ 1 ] = Sequence[]


SignedMixers[ fi_ ] := MixingPartners[fi][[-1]] /; FreeQ[fi, Generic]

SignedMixers[ -fi_[x__] ] := -MixingPartners[AntiParticle[fi]][[-1]][x]

SignedMixers[ fi_[x__] ] := MixingPartners[fi][[-1]][x]


ResolveGeneric[ vert:Vertex[_, cto_:0][_], chainprops___ ] :=
Block[ {v, perm},
  v = SignedMixers/@ TakeInc[vert]/@ Flatten[{chainprops, props}];
  perm = FindVertex[ToGeneric[v], Generic];
  If[ perm === $Failed, Return[{}] ];
  v = v[[perm]];
  If[ cto < 0,
    I PV[ If[FreeQ[v, P$NonCommuting], Identity, NonCommutative][
            VertexFunction[-cto]@@ v ] ],
  (* else *)
    AnalyticalCoupling[cto]@@ v ]
]


ResolveChain[ c_[props__] ] := ToChain[c]@@ MidVertex[props]


Attributes[ MidVertex ] = {Flat}

MidVertex[ p1:Propagator[_][__], p2:Propagator[_][v_, __] ] := (
  vert = DeleteCases[vert, v];
  MidVertex[p1, ResolveGeneric[v, p1, p2], p2]
)


LeftVertex[ p1:_[Vertex[1][_], __] ] = p1

LeftVertex[ p1:_[v_, __] ] := (
  vert = DeleteCases[vert, v];
  Seq[ResolveGeneric[v, p1], p1]
)


RightVertex[ p2:_[_, Vertex[1][_], ___] ] = p2

RightVertex[ p2:_[_, v_, ___] ] := (
  vert = DeleteCases[vert, v];
  Seq[p2, ResolveGeneric[v, p2]]
)


ToChain[dot][ p1_, pr___, p2_ ] :=
  FermionChain[LeftVertex[p1], pr, RightVertex[p2]]

ToChain[tr][ p1:_[v_, __], pr___, p2_ ] := (
  vert = DeleteCases[vert, v];
  MatrixTrace[p1, pr, p2, ResolveGeneric[v, p2, p1]]
)

	(* the single-propagator version applies to tadpoles: *)
ToChain[tr][ p1:_[v_, __] ] := (
  vert = DeleteCases[vert, v];
  MatrixTrace[p1, ResolveGeneric[v, p1]]
)


FermionChain[] = MatrixTrace[] = 1	(* e.g. for ghosts *)


	(* if AnalyticalPropagator with the exact type is not defined
	   use a more generic type.  The replacement must be limited to
	   the head or else the kinematical information is altered. *)
ResolvePropagator[type_][ _, _, part_ ] :=
Block[ {res},
  res = MapAt[
    # /. {_Loop -> Internal, Incoming | Outgoing -> External} &,
    AnalyticalPropagator[type][part],
    0 ];
  If[ Head[res] === PV,
    res /. Mass[fi_] -> Mass[fi, ResolveType[type]],
  (* else *)
    Message[CreateFeynAmp::noprop, part]; Propagator[part] ]
]


TakeNC[ f_List ] :=
Block[ {s},
  s = Select[f, FreeQ[#, NonCommutative]&];
  AppendTo[scalars, Select[s, FreeQ[#, PropagatorDenominator]&]];
  AppendTo[prden, Select[s, !FreeQ[#, PropagatorDenominator]&]];
  Sequence@@ Select[f, !FreeQ[#, NonCommutative]&]
]

TakeNC[ f_Times ] := TakeNC[List@@ f]

TakeNC[ f_ ] := TakeNC[{f}]


LoopPD[ p_ ] := p /; FreeQ[p, Internal]

LoopPD[ p_Plus ] := LoopPD/@ p

LoopPD[ p_Times ] :=
  Select[p, FreeQ[#, PropagatorDenominator]&] *
    LoopPD@@ Cases[p, _PropagatorDenominator]

LoopPD[ p__PropagatorDenominator ] :=
  Times@@ Select[{p}, FreeQ[#, FourMomentum[Internal, _]]&] *
    FeynAmpDenominator@@ SortPD/@
      Select[{p}, !FreeQ[#, FourMomentum[Internal, _]]&]


SortPD[ PropagatorDenominator[mom_, mass__] ] :=
  PropagatorDenominator[ Expand[-mom], mass ] /;
  !FreeQ[mom, -FourMomentum[Internal, _]]

SortPD[ p_ ] = p


Attributes[ FeynAmpDenominator ] = {Orderless}

FeynAmpDenominator[ ] = 1


SumOver[ i_, {}, ext___ ] := SumOver[i, 0, ext]

SumOver[ i_, NoUnfold[l_], ext___ ] := SumOver[i, l, ext]

SumOver[ i_, r:{___, l_Integer}, ext___ ] :=
  SumOver[i, l, ext] /; r === Range[l]


CreateAmpIns[ top_, gm_, sgen_, gr_ -> ins_ ] :=
  CreateAmpIns[top, gm, sgen, gr] ->
    (CreateAmpIns[top, gm, sgen, #]&)/@ ins

CreateAmpIns[ top_, gm_, sgen_, gr:FeynmanGraph[s_, ___][ru__] ] :=
Block[ {ext, int, ins, deltas},
  ins = ReplacePart[gm, sgen/s, -1] /. {ru} /.
    anti -> AntiParticle /.
    app[ x_. (fi:P$Generic)[n__], k__ ] :> x fi[n, k];
  deltas = DeleteCases[ Union@@ CouplingDeltas/@
    Union[ Cases[ins, G[_][cto_][fi__][_] :> FieldPoint[cto][fi]] ],
    _Integer ];
  ins = ins /. G -> GtoC /. Mass -> TheMass /. gaugeru /.
    _MTF -> 1 /. Thread[deltas -> 1];
  ext = Union[Cases[Take[{ru}, next], _Index, Infinity]];
  int = Complement[Cases[Drop[{ru}, next], _Index, Infinity], ext];
  ins[[-1]] *=
    graphinfo[FeynmanGraph[ru], top, toplist] *
    Times@@ deltas *
    Times@@ (SumOver[#, IndexRange[Take[#, 1]], External]&)/@ ext *
    Times@@ (SumOver[#, IndexRange[Take[#, 1]]]&)/@ int;
  ins
]


(* about G -> C replacement: 
   GtoC tries to replace the head "G" by "TheC" (the classes permutation
   is first resolved by applying the appropriate mapping of kinematical
   indices to all G-expressions).  Failing that, it will try the negative
   kinematical expression (for a G[-]).  If neither method resolves TheC,
   it will issue a warning and return C[cto][fields][kinpart]. *) 

GtoC[ sym_ ][ cto_ ][ fi__ ][ k_ ] :=
Block[ {vert, perm, ferm, kin, cv, cvr},
  vert = Last/@ MixingPartners/@ {fi};

  perm = FindVertex[ToClasses[vert], Classes];
  If[ perm === $Failed, Return[C[cto][fi][k]] ];

  ferm = Cases[vert, _. _F];
  vert = vert[[perm]];
  kin = k /. MapIndexed[KinRule, perm];
  If[ Length[M$FlippingRules] > 0 && ferm =!= Cases[vert, _. _F],
    If[ Length[ferm] > 2, Message[CreateFeynAmp::ambig, vert] ];
    kin = kin /. M$FlippingRules ];

  cv = SignResolve[sym, cvr = TheC[kin]@@ vert];

  If[ !FreeQ[cv, TheC],
    Message[CreateFeynAmp::nocoupl, vert, kin];
    VertexDebug[{"vert" -> vert, "kin" -> kin, "cv" -> cvr,
      "fi" -> {fi}, "cto" -> cto}];
    Return[C[cto][fi][k]] ];

	(* check requested counter-term order *)
  If[ Length[cv] <= cto,
    Message[CreateFeynAmp::counter, vert, cto];
    Return[C[cto][fi][k]] ];

  cv[[cto + 1]]
]


SignResolve[ -1, TheC[kin_][cv__] ] := -TheC[-kin][cv]

SignResolve[ _, c_ ] = c


KinRule[i_, {i_}] = Sequence[]

KinRule[i_, {j_}] = (obj:Alternatives@@ Prepend[KIs, Mom])[i] -> obj[j]


FindVertex::novert =
"Cannot find vertex `1`."

	(* at classes level there may be several definitions for
	   fermionic vertices, e.g. C[F, -F, ...] and C[-F, F, ...].
	   Thus, an exact match is attempted first and only then the
	   sorted form is used. *)
FindVertex[ v_, Classes ] := Range[Length[v]] /;
  MemberQ[ReferenceOrder[Classes], v]

FindVertex[ v_, lev_ ] :=
Block[ {pos, fp},
  pos = Position[FieldPointList[lev], FieldPoint@@ v, 1, 1];
  If[ Length[pos] === 0,
    Message[FindVertex::novert, v];
    Return[$Failed] ];
  fp = ReferenceOrder[lev][[ pos[[1, 1]] ]];
  fp[[ Ordering[fp] ]] = Ordering[v];
  fp
]


Format[ G[sym_][cto_][fi__] ] :=
  DisplayForm[
    SubsuperscriptBox[ "G",
      StringJoin@@ ToString/@ ToGeneric[{fi}],
      "(" <> ToString[cto] <> ")" ] ]

Format[ NonCommutative[nc__] ] := Dot[nc]

Format[ MatrixTrace ] = "tr"

Format[ FermionChain[a__] ] := Dot[a]

Format[ PropagatorDenominator[p_, m_, d___] ] :=
  Block[ {x = p^2 - m^2}, 1/x^d /; x =!= 0 ]

Format[ FourMomentum[Incoming, i_Integer] ] := SequenceForm["p", i]

Format[ FourMomentum[Outgoing, i_Integer] ] := SequenceForm["k", i]

Format[ FourMomentum[Internal, i_Integer] ] := SequenceForm["q", i]

Format[ FourMomentum[_, s_Symbol] ] = s

Format[ Index[type_, i_] ] :=
  SequenceForm[StringTake[ToString[type], 3], i]


PickLevel[ _ ][ tops_TopologyList ] = tops

PickLevel[ lev_ ][ tops:TopologyList[___][___] ] :=
Block[ {Rule, levels, res},
  _ -> Insertions[_][] = Sequence[];
	(* Generic is always kept *)
  res = Switch[ {FreeQ[lev, Classes], FreeQ[lev, Particles]},
    {True, True},
      tops /. (x_ -> Insertions[Classes | Particles][__]) -> x,
    {False, True},
      tops /. (x_ -> Insertions[Particles][__]) -> x,
    {True, False},
      tops /. gr:Insertions[Classes][__] :>
        Insertions[Particles]@@ Join@@ TakeIns/@ gr,
    _,
      tops ];
  levels = DeleteCases[Flatten[{lev}], Generic];
  If[ Length[levels] =!= 0,
    res = MapAt[ Select[#, ContainsQ[#, levels]&]&, res,
      Array[{#, 2}&, Length[res]] ] ];
  res
]

PickLevel::nolevel =
"Warning: FeynAmps have already been picked at a different level, \
`1` level cannot be extracted."

PickLevel[ lev_ ][ amps:FeynAmpList[___][___] ] :=
Block[ {n = 0, c = 0, warn = True},
  LevelPick[lev]/@ (amps /. Number == _ :> Seq[]) /.
    (AmplitudeLevel -> _) -> (AmplitudeLevel -> {lev})
]

PickLevel[ lev_ ][ amp_FeynAmp ] :=
Block[ {n = 0, c = 0, warn = True},
  FeynAmpList[AmplitudeLevel -> {lev}][
    LevelPick[lev][amp /. Number == _ :> Seq[]] ]
]


LevelPick[ Generic ][ amp_ ] :=
Block[ {RelativeCF = 1},
  If[ Length[amp] =!= 3 || MatchQ[amp[[1]], GraphID[__, Generic == _]],
    Insert[ Take[amp, 3], Number == ++n, {1, -1} ],
  (* else *)
    If[warn, Message[PickLevel::nolevel, Generic]; warn = False];
      Seq[] ]
]

LevelPick[ lev:Classes | Particles ][ amp_ ] := (
  Sequence@@ (Insert[#, Number == ++n, {1, -1}]&)/@
    If[ Length[amp] === 3,
      If[ MatchQ[amp[[1]], GraphID[__, lev == _]], amp,
        If[warn, Message[PickLevel::nolevel, lev]; warn = False]; {} ],
    (* else *)
      ApplyGMRules[Take[amp, 3], amp[[-1]], lev] ]
)


ApplyGMRules[ amp_, gm_ -> Insertions[lev_][ru__], lev_ ] :=
Block[ {n = 0},
  Insert[#, lev == ++n, {1, -1}]&/@
    (amp /. (Thread[gm -> TakeGraph[#]]&)/@ {ru})
]

ApplyGMRules[ amp_, gm_ -> Insertions[Classes][ru__], Particles ] :=
Block[ {partru},
  partru = Flatten[TakeIns/@ Insertions[Particles][ru]];
  If[ Length[partru] === 0, {},
    ApplyGMRules[ Insert[amp, Classes == ++c, {1, -1}],
      gm -> partru, Particles ] ]
]


ExpandRanges[n__] := List/@ Union[Flatten[ {n} /.
  a_Integer (Repeated | RepeatedNull)[b_] :>
    Range@@ Sort[Floor[{b, a}]] ]]


AndLower[this:{r__, _, _}, o___] := AndLower[{r, 1}, this, o]

DiagramExtract[ tops:TopologyList[info__][__], n__ ] :=
Block[ {lev, p, Rule},
  lev = Alternatives@@ ResolveLevel[InsertionLevel /. {info}];
  p = Position[tops, FeynmanGraph[_, lev == _][__]];
  Check[ 
    p = Complement[p,
      Level[AndLower/@ Extract[p, ExpandRanges[n]], {2}]],
    Return[$Failed] ];
  Rule[_] := Sequence[];
  Delete[tops, p] /.
    (FeynmanGraph[__][__] -> _[]) :> Seq[] /.
    (Topology[__][__] -> _[]) :> Seq[]
]

DiagramExtract[ other_, n__ ] := Extract[other, ExpandRanges[n]]


DiagramDelete[ tops:TopologyList[info__][__], n__ ] :=
Block[ {lev, p, Rule},
  lev = Alternatives@@ ResolveLevel[InsertionLevel /. {info}];
  p = Position[tops, FeynmanGraph[_, lev == _][__]];
  Check[ p = Extract[p, ExpandRanges[n]], Return[$Failed] ];
  Rule[_] := Sequence[];
  Delete[tops, p] /.
    (FeynmanGraph[__][__] -> _[]) :> Seq[] /.
    (Topology[__][__] -> _[]) :> Seq[]
]

DiagramDelete[ other_, n__ ] := Delete[other, ExpandRanges[n]]


DiagramMap[ foo_, tops:(h_TopologyList)[__] ] :=
Block[ {lev, Rule},
  lev = ResolveLevel[InsertionLevel /. List@@ h][[-1]];
  Rule[_] := Sequence[];
  Apply[ #1 -> (#2 /. gr:FeynmanGraph[__, lev == _][__] :> foo[gr, #1, h])&,
      tops, 1 ] /.
    (FeynmanGraph[__][__] -> _[]) :> Seq[] /.
    (Topology[__][__] -> _[]) :> Seq[]
]

DiagramMap[ foo_, other_ ] := foo/@ other


DiagramSelect[ tops:TopologyList[___][__], crit_ ] :=
Block[ {sel},
  sel[g_, r__] := g /; crit[g, r];
  _sel = Sequence[];
  DiagramMap[sel, tops]
]

DiagramSelect[ other_, crit_ ] := Select[other, crit]


DiagramGrouping[ tops:TopologyList[___][___], foo_ ] :=
Block[ {Rule, tag, group, res, c = 0},
  Rule[_] := Sequence[];
  tag[new_] := tag[new] = group[++c];
  res = DiagramMap[tag[foo[##]][#1]&, tops];
  Cases[ DownValues[tag], _[_[_[obj_]], g_group] :>
    obj -> (res /. {g[gr_] :> gr, group[_][_] :> Seq[]}) ] /.
    (FeynmanGraph[__][__] -> _[]) :> Seq[] /.
    (Topology[__][__] -> _[]) :> Seq[]
]


Attributes[ prop ] = {Orderless}

Attributes[ merge ] = {Flat, Orderless}

merge[ prop[i_, j_], prop[j_, k_] ] := prop[i, k]

FermionRouting[ gr_:{}, top:P$Topology, ___ ] := Level[
  merge@@ Apply[ prop[ #1[[1]], #2[[1]] ]&,
    Select[AddFieldNo[top] /. List@@ gr, !FreeQ[#, P$NonCommuting]&],
    1 ],
  {-1} ]


toins[ Generic, gr_ ] := Insertions[Generic][gr]

toins[ lev_, gr_ ] := Insertions[Generic][gr -> Insertions[lev][gr]]

_FeynAmpCases[ _[], ___ ] = {}

(fac_FeynAmpCases)[ gr:FeynmanGraph[__, lev_ == _][__],
    top:P$Topology, h_ ] :=
  fac @ FeynAmpExpr[gr, top, h]

FeynAmpCases[ patt_, lev_:Infinity ][ amp_ ] := Cases[amp, patt, lev]


FeynAmpExpr[ gr:FeynmanGraph[__, lev_ == _][__], top:P$Topology, h_ ] :=
Block[ {$Verbose = 0},
  CreateFeynAmp[ h[top -> toins[lev, gr]],
    AmplitudeLevel -> {lev} ]
]


Unprotect[Exponent]

Exponent[ FeynAmp[_, _, amp_], sym_ ] :=
  Exponent[amp /. FermionChain | MatrixTrace -> Times, sym]

Exponent[ FeynAmp[_, _, amp_, gm_ -> ins_], sym_ ] :=
Block[ {tamp = amp /. FermionChain | MatrixTrace -> Times},
  Exponent[tamp /. gm -> #, sym]&/@ ins
]

Protect[Exponent]


DiagramComplement[ tops:TopologyList[info__][___],
  more:TopologyList[__][___].. ] :=
Block[ {lev = ResolveLevel[InsertionLevel /. {info}][[-1]]},
  Fold[DiagramRemove[#2]/@ #1 &, tops, Level[{more}, {2}]]
]

DiagramRemove[ top_ -> rem_ ][ top_ -> ins_ ] := top -> (
Block[ {Rule, FeynmanGraph},
  Rule[_] := Sequence[];
  Cases[rem,
    FeynmanGraph[_, lev == _][gr__] :>
      (FeynmanGraph[_, lev == _][gr] = Sequence[]),
    Infinity];
  ins
] /. (FeynmanGraph[__][__] -> _[]) :> Seq[] ) /.
  (Topology[__][__] -> _[]) :> Seq[]

DiagramRemove[ _ ][ t_ ] = t


ToJoin[ h:Topology == _, r__ ] := {h, ToJoin[r]}

ToJoin[ r__, h:Number == _ ] := {ToJoin[r], h}


ToFA1Conventions[ expr_ ] :=
Block[ {GraphID, FourMomentum, Conjugate, Global`PolarizationVector,
Global`DiracSpinor, Index, Integral = Sequence, FermionChain = Dot,
NonCommutative = Dot, MatrixTrace = Global`DiracTrace},

  GraphID[ id__ ] := Global`GraphName@@
    Apply[StringJoin, Flatten[{"", ToJoin[id]}] /.
      lev_ == n_ :> {StringTake[ToString[lev], 1], ToString[n]}, 1];

  FourMomentum[ Incoming | External, n_Integer ] :=
    FourMomentum[ Incoming | External, n ] =
      ToExpression["p" <> ToString[n]];
  FourMomentum[ Outgoing, n_Integer ] :=
    FourMomentum[ Outgoing, n ] =
      ToExpression["k" <> ToString[n]];
  FourMomentum[ Internal, n_Integer ] :=
    FourMomentum[ Internal, n ] =
      ToExpression["q" <> ToString[n]];

  Conjugate[Global`PolarizationVector][ args__ ] :=
    Conjugate[ Global`PolarizationVector[args] ];
  Global`PolarizationVector[ _, mom_, li_ ] =
    Global`PolarizationVector[mom, li];

  (* Global`DiracSpinor[ mom_, mass_, ___ ] := FeynArts`Spinor[mom, mass]; *)

  Index[ Global`Lorentz, n_ ] := Index[Global`Lorentz, n] =
    ToExpression["li" <> ToString[n]];

  expr
]

End[]

