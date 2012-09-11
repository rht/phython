(*

This is FeynArts, Version 3.7
Copyright by Sepp Kueblbeck, Hagen Eck, and Thomas Hahn 1991-2012
last modified 11 Jul 12 by Thomas Hahn

Release notes:

FeynArts is free software, but is not in the public domain.
Instead it is covered by the GNU Lesser General Public License.
In plain English this means:

1. We don't promise that this software works.
   (But if you find any bugs, please let us know!)

2. You can use this software for whatever you want.
   You don't have to pay us.

3. You may not pretend that you wrote this software.
   If you use it in a program, you must acknowledge
   somewhere in your publication that you've used
   our code.

If you're a lawyer, you can find the legal stuff at
http://www.fsf.org/copyleft/lgpl.html.

The user guide for this program can be found at
http://www.feynarts.de.

If you find any bugs, or want to make suggestions, or
just write fan mail, address it to:
	Thomas Hahn
	Max-Planck-Institute for Physics
	Foehringer Ring 6
	D-80805 Munich, Germany
	e-mail: hahn@feynarts.de

There exists a low-traffic mailing list where updates will be
announced.  Contact hahn@feynarts.de to be added to this list. 

Have fun!

*)


Print[""];
Print["FeynArts 3.7"];
Print["by Hagen Eck, Sepp Kueblbeck, and Thomas Hahn"];
Print["last revised 11 Jul 12"]


BeginPackage["FeynArts`"]

(* definitions for Utilities.m *)

FAPrint::usage =
"FAPrint[l, s] prints s if l <= $Verbose."

ActualOptions::usage =
"ActualOptions[sym, options] returns a list of options of sym with the
valid options of sym replaced by their actual values."

ResolveLevel::usage =
"ResolveLevel[lev] returns a full set of levels selected by lev.  For
example, ResolveLevel[Particles] gives {Generic, Classes, Particles}."

ResolveType::usage =
"ResolveType[t] returns an abridged class of propagator with values
External, Internal, or Loop."

ContainsQ::usage =
"ContainsQ[expr, items] gives True if expr contains every element in
items."

ToGeneric::usage =
"ToGeneric[expr] returns expr with all classes and particle fields
replaced by their generic fields.  Mind that this procedure removes the
signs of the fields."

ToClasses::usage =
"ToClasses[expr] returns expr with all particle fields replaced by their
classes fields."

Seq::usage =
"Seq is almost identical to Sequence except that it is not expanded
automatically."

TakeGraph::usage =
"TakeGraph[ins -> graph] returns graph."

TakeIns::usage =
"TakeIns[ins -> graph] returns ins."

Subst::usage =
"Subst[expr, i, j] substitutes the elements of i by the corresponding
elements of j, where j may be shorter than i."

PSort::usage =
"PSort[p] sorts the first two elements of a propagator p."

VSort::usage =
"VSort[v] sorts vertex v into canonical order."

Vertices::usage =
"Vertices[top] returns a list of all vertices of topology top."

AddFieldNo::usage =
"AddFieldNo[top] adds numbers of the form Field[n] to the propagators in
topology top."

Compare::usage =
"Compare[t] is the `pure' compare function to eliminate equivalent
topologies of a TopologyList t."

ProcessName::usage =
"ProcessName[amp] constructs a string suitable as filename for the
inserted topology or amplitude list amp which is unique to the model
and particle selection."

Pluralize::usage =
"Pluralize is an internal function."

Statistics::usage =
"Statistics is an internal function."

Alph::usage =
"Alph[n] gives the nth lowercase letter."

UCAlph::usage =
"UCAlph[n] gives the nth uppercase letter."

Greek::usage =
"Greek[n] gives the nth lowercase greek letter."

UCGreek::usage =
"UCGreek[n] gives the nth uppercase greek letter."


(* definitions for Topology.m *)

Topology::usage =
"Topology is the head of a topology data structure.  Topology[s] is
the head of a topology with combinatorial factor 1/s."

TopologyList::usage =
"TopologyList is the head of a list of topologies."

Propagator::usage =
"Propagator[v1, v2] is an (undirected) propagator joining the two 
vertices v1 and v2.\n
Propagator[v1, v2, f] is a (directed) propagator transporting field f
from vertex v1 to v2.\n
Propagator[t][...] is the representation of a propagator of type t.
Possible types are: Incoming, Outgoing, External, Internal and Loop[n]."

Incoming::usage =
"Propagator[Incoming][...] denotes an incoming propagator."

Outgoing::usage =
"Propagator[Outgoing][...] denotes an outgoing propagator."

External::usage =
"Propagator[External][...] denotes an external propagator.\n
External is also used in SumOver[index, range, External] to indicate
that the index to be summed belongs to an external particle."

Internal::usage =
"Propagator[Internal][...] denotes an internal propagator.\n
ExcludeTopologies -> Internal excludes topologies that contain
Propagator[Internal], i.e. the one-particle reducible ones."

Loop::usage =
"Propagator[Loop[n]][...] denotes a propagator on loop n."

Vertex::usage =
"Vertex[e][n] is the representation of a vertex with e propagators in
a topology.  Vertex[e, cto][n] is the representation of a vertex of
counter-term order cto in a topology."

CreateTopologies::usage =
"CreateTopologies[l, i -> o] returns a TopologyList of topologies with
i incoming and o outgoing legs and l loops (at the moment l = 0...3)."

ExcludeTopologies::usage =
"ExcludeTopologies is an option of CreateTopologies which specifies
filters for excluding topologies.  You may use the built-in filters:
Tadpoles, TadpoleCTs, SelfEnergies, SelfEnergyCTs, WFCorrections,
WFCorrectionCTs, Triangles, TriangleCTs, Boxes[n], BoxCTs[n], AllBoxes,
AllBoxCTs, or define new ones using $ExcludeTopologies."

$ExcludeTopologies::usage =
"$ExcludeTopologies[filt] is the function corresponding to filt in
ExcludeTopologies -> {filt, ...}.  It must be defined as a pure function
(i.e. func[#]&).  Given a topology, this function must return True if
the topology shall not be discarded."

Loops::usage =
"Loops[n] is a filter used with ExcludeTopologies.  It excludes
topologies containing loops that are connected to the rest of the graph
with n propagators."

CTs::usage =
"CTs[n] is a filter used with ExcludeTopologies.  It excludes the
counter-term topologies corresponding to Loops[n]."

Tadpoles::usage =
"Tadpoles is a filter used with ExcludeTopologies.  It excludes
topologies containing loops that are connected to the rest of the graph
with one propagator."

TadpoleCTs::usage =
"TadpoleCTs is a filter used with ExcludeTopologies.  It excludes
counter-term topologies corresponding to Tadpoles."

SelfEnergies::usage =
"SelfEnergies is a filter used with ExcludeTopologies.  It excludes
topologies containing loops that are connected to the rest of the graph
with two propagators."

SelfEnergyCTs::usage =
"SelfEnergyCTs is a filter used with ExcludeTopologies.  It excludes
counter-term topologies corresponding to SelfEnergies."

WFCorrections::usage =
"WFCorrections is a filter used with ExcludeTopologies.  It excludes
wave-function-correction topologies, i.e. self-energy insertions and
tadpoles on external legs.  It can also be used in the form
WFCorrections[patt], which limits the selection to only those external 
lines matching patt."

WFCorrectionCTs::usage =
"WFCorrectionCTs is a filter used with ExcludeTopologies.  It excludes
wave-function-correction counter-term topologies, i.e. the counter-terms
corresponding to WFCorrections.  It can also be used in the form
WFCorrectionCTs[patt], which limits the selection to only those external
lines matching patt."

Triangles::usage =
"Triangles is a filter used with ExcludeTopologies.  It excludes
topologies containing loops that are connected to the rest of the graph
with three propagators."

TriangleCTs::usage =
"TriangleCTs is a filter used with ExcludeTopologies.  It excludes the
counter-term topologies corresponding to Triangles."

Boxes::usage =
"Boxes is a filter used with ExcludeTopologies.  It excludes
topologies containing loops that are connected to the rest of the graph
with four propagators."

BoxCTs::usage =
"BoxCTs is a filter used with ExcludeTopologies.  It excludes the
counter-term topologies corresponding to Boxes."

Pentagons::usage =
"Pentagons is a filter used with ExcludeTopologies.  It excludes
topologies containing loops that are connected to the rest of the graph
with five propagators."

PentagonCTs::usage =
"PentagonCTs is a filter used with ExcludeTopologies.  It excludes the
counter-term topologies corresponding to Pentagons."

Hexagons::usage =
"Hexagons is a filter used with ExcludeTopologies.  It excludes
topologies containing loops that are connected to the rest of the graph
with six propagators."

HexagonCTs::usage =
"HexagonCTs is a filter used with ExcludeTopologies.  It excludes the
counter-term topologies corresponding to Hexagons."

AllBoxes::usage =
"AllBoxes is a filter used with ExcludeTopologies.  It excludes
topologies containing loops that are connected to the rest of the graph
with four or more propagators."

AllBoxCTs::usage =
"AllBoxCTs is a filter used with ExcludeTopologies.  It excludes
counter-term topologies corresponding to AllBoxes."

TadpolesOnly::usage =
"TadpolesOnly is a short-cut used with CreateTopologies to keep only
topologies containing tadpoles."

TadpoleCTsOnly::usage =
"TadpoleCTsOnly is a short-cut used with CreateCTTopologies to keep only
topologies containing one-point counter-terms."

SelfEnergiesOnly::usage =
"SelfEnergiesOnly is a short-cut used with CreateTopologies to keep only
topologies containing self-energies on internal lines."

SelfEnergyCTsOnly::usage =
"SelfEnergyCTsOnly is a short-cut used with CreateCTTopologies to keep
only topologies containing two-point counter-terms on internal lines."

TrianglesOnly::usage =
"TrianglesOnly is a short-cut used with CreateTopologies to keep only
topologies containing triangles."

TriangleCTsOnly::usage =
"TriangleCTsOnly is a short-cut used with CreateCTTopologies to keep
only topologies containing three-point counter-terms."

BoxesOnly::usage =
"BoxesOnly is a short-cut used with CreateTopologies to keep only
topologies containing boxes."

BoxCTsOnly::usage =
"BoxCTsOnly is a shortcut used with CreateCTTopologies to keep only
topologies containing four-point counter-terms."

PentagonsOnly::usage =
"PentagonsOnly is a short-cut used with CreateTopologies to keep only
topologies containing pentagons."

PentagonCTsOnly::usage =
"PentagonCTsOnly is a shortcut used with CreateCTTopologies to keep
only topologies containing five-point counter-terms."

HexagonsOnly::usage =
"HexagonsOnly is a short-cut used with CreateTopologies to keep only
topologies containing hexagons."

HexagonCTsOnly::usage =
"HexagonCTsOnly is a shortcut used with CreateCTTopologies to keep
only topologies containing six-point counter-terms."

ToTree::usage =
"ToTree[top] returns top with the loops shrunk to points named
Centre[adj][n] where adj is the adjacency of loop n."

Centre::usage =
"Centre[adj][n] represents the remains of loop n with adjacency adj
after being shrunk to a point by ToTree."

FreeWFQ::usage =
"FreeWFQ[top, patt1, patt2] determines if the topology top is free of
one-point vertices specified by patt1 and two-point vertices specified
by patt2 on external legs.  For example, the WFCorrections filter uses
FreeWFQ[ToTree[top], Centre[1], Centre[2]]&."

LoopFields::usage =
"LoopFields[top] returns a list of the fields that are part of any
loop in the topology top.  LoopFields[rul, top] first substitutes
the insertion rules rul into the bare topology top before proceeding. 
This function is typically used as a filter for DiagramSelect or
DiagramGrouping, as in
DiagramSelect[diags, FreeQ[LoopFields[##], V[1]]&]."

WFCorrectionFields::usage =
"WFCorrectionFields[top] extracts the fields external to any
wave-function correction from topology top.
WFCorrectionFields[rul, top] first substitutes the insertion rules rul
into the bare topology top before proceeding.  This function is 
typically used as a filter for DiagramSelect, as in  
DiagramSelect[diags, UnsameQ@@ WFCorrectionFields[##] &]."

WFCorrectionCTFields::usage =
"WFCorrectionCTFields[top] extracts the fields external to any
wave-function-correction counter-term from topology top.
WFCorrectionFields[rul, top] first substitutes the insertion rules rul
into the bare topology top before proceeding.  This function is 
typically used as a filter for DiagramSelect, as in  
DiagramSelect[diags, UnsameQ@@ WFCorrectionCTFields[##] &]."

StartingTopologies::usage =
"StartingTopologies is an option of CreateTopologies.  It specifies a
pattern for selecting the starting topologies.  The latter are defined
in Topology.m, e.g. at two-loop level there are three starting
topologies: Theta, Eight, and Bicycle."

StartTop::usage =
"StartTop[l, cto] is the list of starting topologies for topologies with
l loops and counter-term order cto.  The starting topologies are defined
in Topology.m."

Theta::usage =
"Theta is the name of the two-loop starting topology that looks like the
greek letter theta."

Eight::usage =
"Eight is the name of the two-loop starting topology that looks like the
number 8."

Bicycle::usage =
"Bicycle is the name of the two-loop starting topology that looks (sort
of) like a bicycle."

Three::usage =
"Three[n] is the name for the irreducible three-loop starting
topologies, where n = 1...8."

ThreeRed::usage =
"ThreeRed[n] is the name for the reducible three-loop starting
topologies, where n = 1...7."

Adjacencies::usage =
"Adjacencies is an option of CreateTopologies.  Its setting is a list
{e1, e2, ...} of integers (ei > 2) of allowed adjacencies of vertices. 
The adjacency of a vertex is the number of propagators ending at that
vertex.  The two special cases ei = 1 and 2 are for external particles
and counter terms, respectively, and are taken care of by
CreateTopologies."

CTOrder::usage =
"CTOrder is an option of CreateTopologies that specifies for which order
counter terms shall be generated."

CreateCTTopologies::usage =
"CreateCTTopologies[l, i -> o] generates all counter-term topologies
needed for the l-loop diagrams generated by CreateTopologies[l, i -> o]."

CreateVFTopologies::usage =
"CreateVFTopologies[l, i -> o] generates all topologies with 1PI vertex
functions whose total loop order is l."

TopologySort::usage = "TopologySort[top] sorts the topology top into a
(more or less) canonical order."

TopologyOrdering::usage =
"TopologyOrdering[top] returns the topology top sorted into a (more or
less) canonical order, together with the permutation that brings it into
this order, i.e. a list of the positions at which the propagators in the
sorted version appeared in the unsorted version.  A negative integer in
the permutation indicates that also the vertices in the respective
propagator were exchanged with respect to the original."

SymmetryFactor::usage =
"SymmetryFactor[top] returns the symmetry factor for the topology top. 
This value is needed if you want to enter new starting topologies."


(* definitions for Insert.m *)

Insertions::usage =
"Insertions is the head of an insertion list.  Insertions[lev] specifies
insertions at level lev.  Insertion lists are returned by InsertFields
as a rule of the form \"topology -> insertionlist\"."

FeynmanGraph::usage =
"FeynmanGraph is the head of a list of field replacement rules.  The
elements of an insertion list are FeynmanGraphs."

Field::usage =
"Field[n] denotes the nth field in a topology."

InsertFields::usage =
"InsertFields[top, {inc1, inc2, ...} -> {out1, out2, ...}] constructs
all Feynman diagrams for the Topology or TopologyList top with incoming
fields inc1, inc2, ... and outgoing fields out1, out2, ..."

Restrictions::usage =
"Restrictions is an option of InsertFields.  It contains shorthands to
exclude vertices or particles defined in the corresponding model file."

LastSelections::usage =
"LastSelections is an option of InsertFields.  It is given as a list of
symbols which must (or must not, if preceded by \"!\") appear in the
Insertions."

InsertionLevel::usage =
"InsertionLevel is an option of InsertFields and CreateFeynAmp. 
Possible values are Generic, Classes, or Particles.  Just as with the
usual Mathematica level specification, e.g. {Particles} means \"only
Particles level\" whereas Particles means \"down to Particles level\". 
By default, CreateFeynAmp uses the same level as InsertFields."

Generic::usage =
"Generic denotes the generic (general field types) level of insertion."

Classes::usage =
"Classes denotes the classes (multiplets) level of insertion."

Particles::usage =
"Particles denotes the particles (members of classes) level of
insertion."

Model::usage =
"Model -> \"MOD\" is an option of InsertFields to select the classes
model MOD.  The model information is taken from the file MOD.mod."

GenericModel::usage =
"GenericModel -> \"GEN\" is an option of InsertFields and
InitializeModel to select the generic model GEN.  The model information
is taken from the file GEN.gen."

Process::usage =
"Process is returned by InsertFields as an option of TopologyList.  
It specifies the process as a rule \"inparticles -> outparticles\"."

FieldPoints::usage =
"FieldPoints[top] returns a list of the field points contained in
the topology top.  FieldPoints[rul, top] first substitutes the
insertion rules rul into the bare topology top before proceeding. 
This function is typically used as a filter for DiagramSelect or
DiagramGrouping, as in
DiagramSelect[diags, MemberQ[FieldPoints[##], (some field point)]&]."

TakeInc::usage =
"TakeInc[v][p] returns the incoming particle from vertex v in
propagator p."

IndexDelta::usage =
"IndexDelta[i1, i2] is a symbol in the definition of a classes coupling
indicating that the coupling is diagonal in the indices i1 and i2."

IndexEps::usage =
"IndexEps[i1, i2, i3] is the totally antisymmetric symbol in the indices
i1, i2, i3."


(* definitions for Initialize.m *)

LoadGenericModel::usage =
"LoadGenericModel[genname] loads the generic model file(s) genname.gen. 
LoadGenericModel[genname, ext] specifies an explicit extension, i.e.
loads genname.ext."

LoadModel::usage =
"LoadModel[modname] loads the classes model file(s) modname.mod. 
LoadModel[modname, ext] specifies an explicit extension, i.e. loads
modname.ext."

DumpGenericModel::usage =
"DumpGenericModel[genfile] saves the generic model file presently in
memory in genfile.  DumpGenericModel[genfile, syms] includes the symbols
syms in the variables to be saved in genfile."

DumpModel::usage =
"DumpModel[modfile] saves the classes model file presently in memory in
modfile.  DumpModel[modfile, syms] includes the symbols syms in the
variables to be saved in modfile."

InitializeModel::usage =
"InitializeModel[modname] initializes the classes model for the model
modname and the generic model given by the GenericModel option.  The
model information is taken from the file modname.mod.  InitializeModel[]
initializes only the generic model."

Reinitialize::usage =
"Reinitialize is an option of InitializeModel.  InitializeModel will
reinitialize the current model only if Reinitialize is set to True."

GenericModelEdit::usage =
"GenericModelEdit is an option of InitializeModel.  It specifies code
that will be executed directly after loading the generic model, i.e.
before the actual initialization."

ModelEdit::usage =
"ModelEdit is an option of InitializeModel.  It specifies code that will
be executed directly after loading the classes model, i.e. before the
actual initialization."

RestrictCurrentModel::usage =
"RestrictCurrentModel[args] applies a number of ExcludeFieldPoints and
ExcludeParticles restrictions to the current model. 
RestrictCurrentModel[] removes all currently active restrictions."

ExcludeParticles::usage =
"ExcludeParticles is an option of RestrictCurrentModel.  It specifies a
list of fields to exclude from the current model.  Excluding a field at
a particular level automatically excludes all derived fields at lower
levels.  For example, excluding F[3] also excludes all F[3, {...}]. 
Further, the exclusion of a particle always implies exclusion of its
antiparticle."

ExcludeFieldPoints::usage =
"ExcludeFieldPoints is an option of RestrictCurrentModel.  It specified
a list of field points to exclude from the current model.  Excluding a
field point at a particular level automatically excludes derived field
points on lower levels.  For example, excluding FieldPoint[F[3], -F[3],
V[1]] also excludes FieldPoint[F[3, {1}], -F[3, {1}], V[1]] etc. 
Further, the exclusion of a field point always implies exclusion of its
charge-conjugate field point."

FieldMatchQ::usage =
"FieldMatchQ[f, fpatt] returns True if the field f matches the pattern
fpatt and False otherwise.  The matching takes into account field levels,
e.g. F[1] matches F."

FieldMemberQ::usage =
"FieldMemberQ[flist, fpatt] returns True if an element of flist matches
fpatt in the sense that FieldMatchQ returns True."

FieldPointMatchQ::usage =
"FieldPointMatchQ[fp, fppatt] returns True if the field point fp matches
the pattern fppatt and False otherwise.  The matching takes into account
field levels, e.g. F[1] matches F."

FieldPointMemberQ::usage =
"FieldPointMemberQ[fplist, fppatt] returns True if an element of fplist
matches fppatt in the sense that FieldPointMatchQ returns True."

ExcludedQ::usage =
"ExcludedQ[vertlist] gives True if vertlist contains any vertices
currently excluded at Particles level."

PossibleFields::usage =
"PossibleFields[cto][t, fp] returns all possible fields of type t that
fit fp at counter-term order cto.  t may be 0 (allowing all fields) or a
generic field."

CheckFieldPoint::usage =
"CheckFieldPoint[fp] yields True if fp is a valid field point in the
current model and False otherwise."

AntiParticle::usage =
"AntiParticle[f] returns the antiparticle of f."

TheMass::usage =
"TheMass[p] gives the value of the mass of particle p (if specified in
the model file).  TheMass[p, t] gives the mass for particle p running
on a propagator of type t (External, Internal, Loop)."

Indices::usage =
"Indices[c] gives a list of index names of class c."

IndexRange::usage =
"IndexRange[i] gives a list of possible values of index i."

NoUnfold::usage =
"NoUnfold, when wrapped around the right hand side of an IndexRange
assignment, prevents InsertFields from \"unfolding\" that index at
particles level, i.e. InsertFields then does not generate an extra
diagram for every value the index can take on."

IndexSum::usage =
"IndexSum[expr, {i, range}] represents the unevaluated sum of expr in
the index i over range.  To execute the sum, replace IndexSum by Sum."

AddHC::usage =
"AddHC[mat] extends mat by its Hermitian conjugate part, e.g.
AddHC[A[1, i, j]] returns (A[1, i, j] + Conjugate[A[1, j, i]])/2.  
AddHC[mat, w] forms the weighted sum with weight function w, e.g.
AddHC[A[1, i, j], w] returns (w[i, j] A[1, i, j] + Conjugate[w[j, i]]
Conjugate[A[1, j, i]])/2."

RenConst::usage =
"RenConst[rc] := ... defines the renormalization constant rc."

ReferenceOrder::usage =
"ReferenceOrder[x] gives a list of all field points of the current model
in (unsorted) list form.  x can be Generic or Classes."

FieldPoint::usage =
"FieldPoint[cto][f1, f2, ...] is the representation of a field point of
counter-term order cto with incoming fields f1, f2, ..."

FieldPointList::usage =
"FieldPointList[cto] returns a list of field points of counter-term order
cto in the current model, FieldPointList[Classes] returns a list of field
points of all orders in the current classes model, and
FieldPointList[Generic] returns a list of all generic field points."

KinematicVector::usage =
"KinematicVector[fi] returns the kinematic vector of the coupling of the
fields fi."

CouplingDeltas::usage =
"CouplingDeltas[v] returns a list of the IndexDeltas in which the entire
vertex v is diagonal."

F$Generic::usage =
"F$Generic gives the list of generic fields of the current model.  
Its contents may change with every call to RestrictCurrentModel."

F$Classes::usage =
"F$Classes gives the list of classes fields of the current model.  
Its contents may change with every call to RestrictCurrentModel."

F$Particles::usage =
"F$Particles gives the list of particles fields of the current model. 
Its contents may change with every call to RestrictCurrentModel."

F$AllGeneric::usage =
"F$AllGeneric gives the list of generic fields of the current model."

F$AllClasses::usage =
"F$AllClasses gives the list of classes fields of the current model."

F$AllParticles::usage =
"F$AllParticles gives the list of particles fields of the current
model."

F$AllowedFields::usage =
"F$AllowedFields is the list of all fields at all three levels that are
in the current model.  Its contents may change with every call to
RestrictCurrentModel."

L$CTOrders::usage =
"L$CTOrders is the list of counter-term orders of the current model."

Mom::usage =
"Mom[n] is the momentum of the nth field in the kinematic vector."

KI1::usage = KI2::usage = KI3::usage = KI4::usage =
"KIi[n] is the ith kinematic index of the nth field in the kinematic
vector."

KIs = {KI1, KI2, KI3, KI4}

SI::usage = "SI[n] is the nth summation index in a component of the
kinematic vector."

SIs = {SI1_, SI2_, SI3_, SI4}

CI::usage = "CI[n] is the classes index of the nth field in the
kinematic vector."

PV::usage =
"PV is the head of a general analytical expression in FeynArts.  
The letters stand for Propagator/Vertex."

F::usage =
"F is a fermion field."

S::usage =
"S is a scalar field."

SV::usage =
"SV is a scalar-vector mixing field."

VS::usage =
"VS is a vector-scalar mixing field.  This field is for internal
purposes only.  Do not use it in a model file."

U::usage =
"U is a ghost field (Grassmann-valued scalar field)."

V::usage =
"V is a vector boson field."

T::usage =
"T is a tensor field."

$SVMixing::usage =
"$SVMixing determines whether mixing of scalar and vector fields is
allowed or not."

$CounterTerms::usage =
"$CounterTerms determines whether or not counter-term couplings are
initialized during InitializeModel.  Note that once a model is
initialized, if you change the value of $CounterTerms you have to
re-initialize the model for this to take effect."

$ModelPath::usage =
"$ModelPath is a list of directory names which is searched for model
files."

$Model::usage =
"$Model is the currently initialized classes model."

$GenericModel::usage =
"$GenericModel is the currently initialized generic model."

ModelDebug::usage =
"ModelDebug can be wrapped around a model name to enable debugging for
that model."

$ModelDebug::usage =
"$ModelDebug determines whether changes introduced by add-on model files
will be reported as a model is initialized.  It can be set to True, in
which case debugging output will be generated for all add-on model
files, or to the name (or list of names) of the add-on model file(s) to
be debugged."

$ModelDebugForm::usage =
"$ModelDebugForm specifies the output form for debugging output when
$ModelDebug = True is set."

$ModelAdded::usage =
"$ModelAdded contains the couplings added by the most recently
initialized add-on model file if $ModelDebug = True."

$ModelChanged::usage =
"$ModelAdded contains the couplings changed by the most recently
initialized add-on model file if $ModelDebug = True."

$ModelRemoved::usage =
"$ModelRemoved contains the couplings removed by the most recently
initialized add-on model file if $ModelDebug = True."

$ExcludedFPs::usage =
"$ExcludedFPs is the list of currently excluded Generic- and
Classes-level field points."

$ExcludedParticleFPs::usage =
"$ExcludedParticleFPs is the list of currently excluded Particles-level
field points."


(* definitions for the model files *)

M$GenericPropagators::usage =
"M$GenericPropagators is the list of propagators and their analytical
expressions in the generic model file."

M$GenericCouplings::usage =
"M$GenericCouplings is the list of couplings and their analytical
expressions in the generic model file."

M$ClassesDescription::usage =
"M$ClassesDescription is the list of classes properties in the classes
model file."

M$CouplingMatrices::usage =
"M$CouplingMatrices is the list of explicit coupling matrices of the
current model."

M$TruncationRules::usage =
"M$TruncationRules is a set of rules that is applied by CreateFeynAmp if
Truncated -> True.  It is defined in the generic model file.  Typically,
it removes external wavefunctions."

M$FlippingRules::usage =
"M$FlippingRules is a set of rules defined in the generic model file
that specify how noncommutative structures in a coupling change when the
order of the corresponding fields changes.  For example, it specifies
how to derive the coupling C[F, -F, ...] from C[-F, F, ...]."

M$LastModelRules::usage =
"M$LastModelRules is a set of rules that is applied by CreateFeynAmp at
the very end.  It is defined in the model file and can e.g. contain
mappings of symbols to certain Contexts or to special symbols of other
packages."

M$LastGenericRules::usage =
"M$LastGenericRules is a set of rules that is applied by CreateFeynAmp
after creation of the generic amplitudes.  It is defined in the generic
model file and can contain e.g. mappings of symbols to certain Contexts
or to special symbols of other packages."

AnalyticalCoupling::usage =
"AnalyticalCoupling[vertex] == G[_][__] . {__} is the form a vertex is
specified in a generic model file."

AnalyticalPropagator::usage =
"AnalyticalPropagator[type][field] == expr is the form a propagator is
specified in a generic model file."

KinematicIndices::usage =
"KinematicIndices[fi] is a list of the kinematic indices the generic
field fi carries.  For example, KinematicIndices[V] = {Lorentz} in
Lorentz.gen."

FieldNumber::usage =
"FieldNumber[f] can be used in the AnalyticalPropagator and
AnalyticalCoupling definitions of the Generic model file to find out
the ordinal number of a field in the diagram being inserted."

MatrixTraceFactor::usage =
"MatrixTraceFactor -> n is an optional entry for fermions in the 
M$ClassesDescription list.  A MatrixTrace (a closed loop of fermions) is 
multiplied by the MatrixTraceFactor of its consituents.  Typical usage 
is MatrixTraceFactor -> 3 for quarks.  The factor may contain objects of 
the form Index[type], which are substituted by the actual indices 
carried by the fields in the loop."

SelfConjugate::usage =
"SelfConjugate -> True | False is an entry in the M$ClassesDescription
list.\n
SelfConjugate[p] is True if field p is self-conjugate and False
otherwise."

Mixture::usage =
"Mixture -> lc is an optional entry in the M$ClassesDescription list
which specifies that the field is the linear combination lc of other
fields."

InsertOnly::usage =
"InsertOnly is an entry in the M$ClassesDescription list.  It specifies
the types of progators the particle may be inserted into.  If not
explicitly specified, the particle may be inserted into all types of
propagators.\n
InsertOnly[p] returns the types of propagators in which field p may be
inserted into."

Mass::usage =
"Mass[p] denotes the mass of particle p.  Mass[p, t] denotes the mass of
particle p running on a propagator of type t (External, Internal, Loop). 
It is just a symbol carrying no further information.  FeynArts defines
the function TheMass that returns the explicit symbol of the mass (if
specified by the model file)."

MixingPartners::usage =
"MixingPartners -> {...} is an entry in the M$ClassesDescription list
for mixing fields and specifies their mixing partners.\n
MixingPartners[p] returns the partners of mixing field p."

QuantumNumbers::usage =
"QuantumNumbers -> {...} is an entry in the M$ClassesDescription list. 
It lists the quantum numbers a particle possesses.  Identifiers for the
quantum numbers can be chosen freely, e.g. Charge.  The quantum numbers
are needed only when 1PI vertex-function insertions are generated (i.e.
for topologies created with CreateVFTopologies) to weed out such vertex
function that violate the conservation of a quantum number, as
determined by ViolatesQ."

ViolatesQ::usage =
"ViolatesQ is a function defined in the model file that determines
whether a vertex function violates quantum-number conservation.  It is
called by InsertFields for every 1PI vertex-function insertion (i.e.
when inserting topologies created with CreateVFTopologies).  It receives
as arguments the quantum numbers of the involved fields (times -1 for
antiparticles) and must return True if the vertex violates the
conservation of those quantum numbers."

Compatibles::usage =
"Compatibles[p] is a list of particles that are compatible with p."

SVCompatibles::usage =
"SVCompatibles[p] is a list of SV particles that are compatible with p."

Index::usage =
"Index is the head of an index name (i.e. Index[Generation])."

PropagatorType::usage =
"PropagatorType is an option used in the M$ClassesDescription list of a
model file and specifies the type of propagator for a particle. 
Possible values are Straight, Sine, Cycles, ScalarDash, and GhostDash. 
For a mixing propagtor, a list with two propagator types may be
specified."

Straight::usage =
"Straight selects a straight line in the PropagtorType option."

Sine::usage =
"Sine selects a wavy line in the PropagatorType option."

Cycles::usage =
"Cycles selects a cycloid in the PropagatorType option."

ScalarDash::usage =
"ScalarDash selects a dashed line in the PropagatorType option."

GhostDash::usage =
"GhostDash selects a dotted line in the PropagatorType option."

PropagatorArrow::usage =
"PropagatorArrow -> Forward, Backward, None is an option used in the
M$ClassesDescription list in the model file."

Forward::usage =
"Forward selects a forward arrow in the PropagatorArrow option."

Backward::usage =
"Backward selects a backward arrow in the PropagatorArrow option."

PropagatorLabel::usage =
"PropagatorLabel is an option used in the M$ClassesDescription list in
the model file.  It is translated into TheLabel statements during model
initialization."

TheLabel::usage =
"TheLabel[p] returns the PropagatorLabel of particle p."

IndexStyle::usage =
"IndexStyle[i] gives the rendering information for the index i.  For
example, IndexStyle[Index[Lorentz, i_]] := Greek[i + 11] makes Lorentz
indices appear as \"\\mu\", \"\\nu\", etc."

TheCoeff::usage =
"TheCoeff[f] gives the list of component fields of which f is a
linear combination together with the respective coefficients.  
For a non-composite field, TheCoeff returns the field itself."

TheC::usage =
"TheC is an internal symbol for storing the coupling matrices."

C::usage =
"C[cto][fields][kinpart] is the symbolic form of a coupling that is
returned when CreateFeynAmp fails to resolve a classes coupling.\n
C[fi] == coupl used as an entry in the M$CouplingMatrices list in the
classes model file defines the coupling of the fields fi."

CC::usage =
"CC[fields] == coupl may be used in the M$CouplingMatrices list in the
classes model file to define the coupling of the fields fi, together
with the conjugate coupling in one line.  The conjugation is performed
with the function ConjugateCoupling which must be defined accordingly. 
If no definition is made, the final amplitudes will simply contain the
symbol ConjugateCoupling."

ConjugateCoupling::usage =
"ConjugateCoupling[coupl] defines how the charge-conjugated coupling is
derived from coupl.  Typically, one I multiplying the coupling constant
must not be conjugated because it derives from the exponent of the path
integral.  If no definition is made for ConjugateCoupling, the final
amplitudes will contain this symbol."

GetCouplings::usage =
"GetCouplings[C[fi], ...] returns all couplings of the form C[fi] == _
in M$CouplingMatrices."

ReplaceCouplings::usage =
"ReplaceCouplings[C[fi] == coup] replaces all couplings matching C[fi]
in M$CouplingMatrices.  All couplings in M$CouplingMatrices, including
the former versions of the ones being replaced, may be used on the
r.h.s. in the form C[fields]."

G::usage =
"G[sym][cto][fields][kinpart] is a generic coupling matrix of
counter-term order cto for fields corresponding to the kinematical
object kinpart.  G is symmetric for sym = 1 and antisymmetric for sym =
-1."


(* definitions for Analytic.m *)

CreateFeynAmp::usage =
"CreateFeynAmp[tops] creates a list of Feynman amplitudes (Head:
FeynAmpList) of all insertions of Topology or TopologyList tops."

AmplitudeLevel::usage =
"AmplitudeLevel is an option of CreateFeynAmp and specifies for which
levels amplitudes are to be created."

GaugeRules::usage =
"GaugeRules is an option of CreateFeynAmp.  It is given as a list of
rules to select a particular gauge.  For example, Feynman gauge is
obtained by GaugeRules -> {_GaugeXi -> 1}."

PreFactor::usage =
"PreFactor is an option of CreateFeynAmp and specifies the prefactor for
each diagram.  In the expression for the prefactor the symbol LoopNumber
can be used which will subsequently be replaced by the actual loop
number."

LoopNumber::usage =
"LoopNumber is used with the option PreFactor of CreateFeynAmp and is
replaced by the loop number of the topology."

Truncated::usage =
"Truncated is an option of CreateFeynAmp that determines whether to
apply the M$TruncationRules of the current generic model."

MomentumConservation::usage =
"MomentumConservation is an option of CreateFeynAmp.  It specifies
whether momentum conservation at each vertex is enforced.  If set to
False, every propagator will carry its own momentum."

GraphInfoFunction::usage =
"GraphInfoFunction is an option of CreateFeynAmp.  It specifies a
function with which every diagram is multiplied.  The function receives
two arguments, f[rul, top], where rul are the insertion rules and top
the corresponding topology of the diagram.  This function can be used
to add graph information to the amplitude."

VertexDebug::usage =
"VertexDebug[debuginfo] is a function invoked whenever a vertex cannot
be resolved.  It is used for debugging FeynArts."

VertexFunction::usage =
"VertexFunction[o][f1, f2, ...] represents the 1PI vertex function of
loop-order o with external fields f1, f2, ..."

FindVertex::usage =
"FindVertex[vert, lev] looks up vertex vert in the currently initialized
model where lev = Generic or Classes.  If a match is found, the
permutation is returned which brings the vertex vert into the reference
order."

PickLevel::usage =
"PickLevel[lev][amp] constructs the concrete amplitudes at level lev
from the CreateFeynAmp result amp.  PickLevel[lev][tops] picks out the
diagrams at level lev from TopologyList tops.  Note that in topology
lists you can never delete the Generic level."

DiagramExtract::usage =
"DiagramExtract[expr, n] extracts diagrams by number from the topology
or amplitude list expr.  n may be of the form 3, 42, 17...28 which
selects diagrams 3, 42, and 17 through 28."

DiagramDelete::usage =
"DiagramDelete[expr, n] discards diagrams by number from the topology
or amplitude list expr.  n may be of the form 3, 42, 17...28 which
selects diagrams 3, 42, and 17 through 28."

Discard[args__] := (
  Message[Discard::obsalt, Discard, DiagramDelete];
  DiagramDelete[args] )

DiagramMap::usage =
"DiagramMap[foo, diags] maps foo over all Feynman diagrams in diags."

DiagramSelect::usage =
"DiagramSelect[diags, crit] selects the diagrams from diags for which
crit is true.  For example, for selecting diagrams where field #5 is not
a S[1], one could use DiagramSelect[diags, FreeQ[#, Field[5] -> S[1]]&]."

DiagramGrouping::usage =
"DiagramGrouping[tops, foo] returns a list of parts of the inserted
topologies tops, grouped according to the output of foo.  For example,
DiagramGrouping[tops, FermionRouting] returns a list of
fermion-flow-ordered diagrams."

FermionRouting::usage =
"FermionRouting[top] finds out which external lines are connected
through fermion lines for the inserted topology top. 
FermionRouting[rul, top] first substitutes the insertion rules rul
into the bare topology top before proceeding.  The output is a list
of integers of which every successive two denote the end-points of
a fermion line in the diagram.  This function is typically used as
a filter for DiagramSelect or DiagramGrouping, as in
DiagramSelect[diags, FermionRouting[##] == {1, 4, 2, 3} &]."

FeynAmpCases::usage =
"FeynAmpCases[patt, lev][amp] finds the instances of pattern patt in the
amplitude amp.  FeynAmpCases[patt, lev][rul, top, h] does the same for
the amplitude resulting from topology top inserted with insertion rules
rul using topology-list header h.  The level specification lev is
optional and defaults to Infinity.  This function is typically used as
a filter for DiagramSelect or DiagramGrouping, as in
DiagramGrouping[diags, FeynAmpCases[_[Index[Colour | Gluon, _], ___]]]."

FeynAmpExpr::usage =
"FeynAmpExpr[rul, top, h] returns the amplitude resulting from topology
top inserted with insertion rules rul using topology-list header h. 
This function is typically used as a filter for DiagramSelect or
DiagramGrouping."

DiagramComplement::usage =
"DiagramComplement[diagall, diag1, diag2, ...] gives all diagrams in
diagall which are not in any of the diagi."

ToFA1Conventions::usage =
"ToFA1Conventions[expr] converts expr back to FeynArts 1 conventions. 
Note that this conversion only renames some symbols.  The output may
thus not be 100% FeynArts 1 compatible since certain kinds of
expressions (e.g. Generic insertions) could not be generated with
FeynArts 1 at all."

$FermionLines::usage =
"$FermionLines is a FeynArts system constant that can be True or False
indicating whether CreateFeynAmp should collect Grassmann-valued fields
(typically F and U) in lines or not.
Fermion line construction has to be turned off for generic models that
contain couplings involving more than two fermions (it is for example
impossible to reliably determine how two fermion lines run through a
four-fermion vertex).  In such a case the fermions should carry a
kinematical (e.g. Dirac) index.
Note that if $FermionLines = False, the M$ClassesDescription option
MatrixTraceFactor has no effect on fermionic classes and also no
additional minus signs are inserted for odd permutations of external
fermions."

$SparseCouplings::usage =
"$SparseCouplings is a FeynArts system constant.  If set to True, a
model is initialized such that all zero-components of the coupling
vector are handled by a single default rule, which speeds up lookup of
Feynman rules in the presence of large but sparse coupling vectors.
Use with care, as individual checks of the kinematical components are
not possible anymore then."

FeynAmp::usage =
"FeynAmp is the head of a data structure that represents the analytical
expression of a single Feynman graph.  Its members are: graph name, list
of integration momenta, analytical expression of the amplitude."

FeynAmpList::usage =
"FeynAmpList is the head of a list of FeynAmps."

GraphID::usage =
"GraphID is the head of a data structure identifying a single FeynAmp in
a FeynAmpList."

Integral::usage =
"Integral[q] is a member of a FeynAmp data structure and contains the
integration momenta of the amplitude (empty for tree graphs)."

RelativeCF::usage =
"RelativeCF is the relative combinatorial factor of a graph with respect
to the generic graph it was created from.  This symbol appears in
generic amplitudes and its value is defined for every Classes or
Particles insertion."

SumOver::usage =
"SumOver[i, r] indicates that the amplitude it is multiplied with is to
be summed in the index i over the range r.  If r is an integer, it
represents the range {1, 2, ..., r}.  For an index belonging to an
external particle there is a third argument, External."

PropagatorDenominator::usage =
"PropagatorDenominator[p, m] stands for the expression 1/(p^2 - m^2)
that is contained in internal propagators of a Feynman graph. 
PropagatorDenominator[p, m, d] is the denominator raised to the
power d."

FeynAmpDenominator::usage =
"FeynAmpDenominator[d1, d2, ...] contains the PropagatorDenominators 
d1, d2, ... that belong to a loop."

GaugeXi::usage =
"GaugeXi[s] is a gauge parameter with index s."

FourMomentum::usage =
"FourMomentum[s, n] is the nth momentum of type s.  Allowed types are
Incoming, Outgoing, External, and Internal."

NonCommutative::usage =
"NonCommutative is the head of noncommuting objects in a Feynman rule."

MatrixTrace::usage =
"MatrixTrace is the head of a trace of noncommuting objects (i.e. of
symbols with head NonCommutative in the Feynman rules) in closed fermion
loops."

FermionChain::usage =
"FermionChain is the head of a trace of noncommuting objects (i.e. of
symbols with head NonCommutative in the Feynman rules) in open fermion
chains."


(* definitions for Graphics.m *)

Paint::usage =
"Paint[tops] draws a list of inserted or bare topologies."

AutoEdit::usage =
"AutoEdit is an option of Paint.  If True, the topology editor is called
when an unshaped topology is found."

ColumnsXRows::usage =
"ColumnsXRows is an option of Paint.  It specifies how many diagrams are
arranged on a page.  For example, ColumnsXRows -> 3 draws 3 rows of 3
diagrams each, or ColumnsXRows -> {5, 8} draws 8 rows of 5 diagrams
each."

PaintLevel::usage =
"PaintLevel is an option of Paint and specifies for which levels
diagrams are drawn."

SheetHeader::usage =
"SheetHeader is an option of Paint.  Automatic or True selects the
default header (the process for inserted or #in -> #out for bare
topologies), False disables headers, everything else is taken literally
as a header (e.g. SheetHeader -> \"Self-energy diagrams\")."

Numbering::usage =
"Numbering is an option of Paint.  It can take the values Full, Simple,
or None.  Full generates numbering of the form T1 C8 N15 (= topology 1,
classes insertion 8, running number 15), Simple generates just a running
number, and None omits the numbering altogether."

Full::usage =
"Full is a possible choice for the Numbering option of Paint.  It
specifies that the numbering of the diagrams be of the form T1 C8 N15
(= topology 1, classes insertion 8, running number 15)."

Simple::usage =
"Simple is a possible choice for the Numbering option of Paint.  It
specifies that the numbering of the diagrams be just a running number."

FieldNumbers::usage =
"FieldNumbers is an option of Paint.  It it meaningful only for bare
topologies where it specifies whether the field numbers (the n in
Field[n]) are used for labelling the propagators.  This can be helpful
for selecting diagrams."

PropagatorGraphics::usage =
"PropagatorGraphics[type, arrow, label][from, to, height, labelpos] is 
a graphics object representing a propagator."

VertexGraphics::usage =
"VertexGraphics[cto][xy] is a graphics object representing a vertex."

DiagramGraphics::usage =
"DiagramGraphics[title][graphics] is a graphics object representing a
Feynman diagram."

FeynArtsGraphics::usage =
"FeynArtsGraphics[title][sheet1, sheet2, ...] is a graphics object
representing Feynman diagrams.  The individual sheets are matrices of
DiagramGraphics objects with the rows and columns of the matrix
representing the graphical arrangement of the diagrams into rows and
columns on the output page."

Render::usage =
"Render[g, format] renders the FeynArtsGraphics object g.  The output
is a string for the \"PS\", \"EPS\", and \"TeX\" formats and a regular
Graphics object for all other formats."

ComposedChar::usage =
"ComposedChar[t, sub, sup, bar] represents a label t with subscript sub,
superscript sup, and accent bar.  Typically, bar is something like
\"\\\\bar\", \"\\\\tilde\", etc.  The arguments sub, sup, and bar are
optional but their position is significant.  For example,
ComposedChar[t, Null, sup] is a label with a superscript only."

Shape::usage =
"Shape[tops] edits the shapes of the topologies tops."

ShapeData::usage =
"ShapeData[topcode] is the database of shapes currently in memory.  
It is indexed by the three strings given by TopologyCode."

TopologyCode::usage =
"TopologyCode[top] returns a list of three strings identifying the
topology.  This code is unique as far painting the topology is
concerned."

{TopBottom, LeftRight}	(* visible for printing *)


(* FeynArts system constants *)

$FeynArts::usage =
"$FeynArts contains the version of FeynArts."

$Verbose::usage =
"$Verbose is an integer that determines the extent of run-time messages
in FeynArts.  It ranges from 0 (no messages) to 2 (normal)."

$FADebug::usage =
"$FADebug = True causes various functions to print internal information
for debugging purposes."

$FeynArtsDir::usage =
"$FeynArtsDir points to the directory from which FeynArts was loaded."

$FeynArtsProgramDir::usage =
"$FeynArtsProgramDir points to the directory which contains the FeynArts
program files."

$ShapeDataDir::usage =
"$ShapeDataDir points to the directory which contains the data for
drawing Feynman diagrams."

P$Options::usage =
"P$Options is a pattern for options."

P$Topology::usage =
"P$Topology is the pattern for a topology."

P$Generic::usage =
"P$Generic is the pattern for generic fields."

P$NonCommuting::usage =
"P$NonCommuting is the pattern for the non-commuting generic fields."

P$InsertionObjects::usage =
"P$InsertionObjects matches the objects in the generic amplitude that
will be taken for insertions."


P$Topology = Topology[__] | Topology[_][__]

P$Generic = F | S | V | T | U | SV

P$NonCommuting = F | U

P$InsertionObjects = G[_][_][__][__] | _Mass | _GaugeXi |
  VertexFunction[_][__]

P$Options = (_Rule | _RuleDelayed)...


$FeynArts = 3.7

$FeynArtsDir = DirectoryName[ File /.
  FileInformation[System`Private`FindFile[$Input]] ]

$FeynArtsProgramDir = ToFileName[{$FeynArtsDir, "FeynArts"}]


Get[ ToFileName[$FeynArtsDir, "Setup.m"] ]

If[ FileType["Setup.m"] === File, Get["Setup.m"] ]


Block[ {$Path = {$FeynArtsProgramDir}},
<< Utilities`;
<< Topology`;
<< Initialize`;
<< Insert`;
<< Analytic`;
<< Graphics`
]

EndPackage[]


Null

