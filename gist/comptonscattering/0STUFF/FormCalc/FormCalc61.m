(*

This is FormCalc, Version 6.1
Copyright by Thomas Hahn 1996-2010
last modified 20 Jul 10 by Thomas Hahn

Release notes:

FormCalc is free software, but is not in the public domain.
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
http://www.feynarts.de/formcalc.

If you find any bugs, or want to make suggestions, or
just write fan mail, address it to:
	Thomas Hahn
	Max Planck Institute for Physics
	Foehringer Ring 6
	D-80805 Munich, Germany
	e-mail: hahn@feynarts.de

There exists a low-traffic mailing list where updates will be
announced.  Contact hahn@feynarts.de to be added to this list.

Have fun!

*)

Print[""];
Print["FormCalc 6.1"];
Print["by Thomas Hahn"];
Print["last revised 20 Jul 10"]


(* symbols from FeynArts *)

BeginPackage["FeynArts`"]

{ FeynAmp, FeynAmpList, Process, GraphID,
  Generic, Classes, Particles, S, F, U, V, T,
  Insertions, G, Mass, GaugeXi, VertexFunction,
  PropagatorDenominator, FeynAmpDenominator,
  FourMomentum, Internal, External, TheMass,
  Index, IndexDelta, IndexEps, IndexSum, SumOver,
  MatrixTrace, FermionChain, NonCommutative,
  CreateTopologies, ExcludeTopologies, Tadpoles,
  InitializeModel, $Model, Model, GenericModel, Reinitialize,
  InsertFields, InsertionLevel, ExcludeParticles,
  ExcludeFieldPoints, LastSelections, Restrictions,
  CreateFeynAmp, Truncated, RenConst, ProcessName,
  Paint, DiagramGrouping }

EndPackage[]


(* symbols from LoopTools *)

BeginPackage["LoopTools`"]

A0i (* need this internally *)

{ bb0, bb1, bb00, bb11, bb001, bb111, dbb0, dbb1, dbb00, dbb11,
  cc0, cc00, cc001, cc002, cc0000, cc0011, cc0022, cc0012,
  dd0, dd0000, dd00001, dd00002, dd00003, ee0, ff0 }

A0::usage =
"A0[m] is the one-point one-loop scalar integral.  m is the mass
squared."

A00::usage =
"A00[m] is the one-point tensor coefficient of g_{mu nu}.  m is the mass
squared."

B0i::usage =
"B0i[id, p, m1, m2] is the generic two-point one-loop integral which
includes scalar and tensor coefficients as well as their derivatives
with respect to p, specified by id.  For example, B0i[bb0, ...] is the
scalar function B_0, B0i[bb11, ...] the tensor coefficient function
B_{11} etc.  p is the external momentum squared and m1 and m2 are the
masses squared."

Bget::usage =
"Bget[p, m1, m2] computes all two-point coefficients in LoopTools."

Bval::usage =
"Bval[id, n] is the array containing the cached two-point integrals in
LoopTools, where id is the coefficient label and n is the index returned
by Bget."

C0i::usage =
"C0i[id, p1, p2, p1p2, m1, m2, m3] is the generic three-point one-loop
integral which includes both scalar and tensor coefficients, specified
by id.  For example, C0i[cc0, ...] is the scalar function C_0,
C0i[cc112, ...] the tensor coefficient function C_{112} etc.  p1, p2,
and p1p2 are the external momenta squared and m1, m2, m3 are the masses
squared."

Cget::usage =
"Cget[p1, p2, p1p2, m1, m2, m3] computes all three-point coefficients in
LoopTools."

Cval::usage =
"Cval[id, n] is the array containing the cached three-point integrals in
LoopTools, where id is the coefficient label and n is the index returned
by Cget."

D0i::usage =
"D0i[id, p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4] is the generic
four-point one-loop integral which includes both scalar and tensor
coefficients, specified by id.  For example, D0i[dd0, ...] is the scalar
function D_0, D0i[dd1233, ...] the tensor function D_{1233} etc. 
p1...p4 are the external momenta squared, p1p2 and p2p3 are the squares
of external momenta (1+2) and (2+3), respectively, and m1...m4 are the
masses squared."

Dget::usage =
"Dget[p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4] computes all
four-point coefficients in LoopTools."

Dval::usage =
"Dval[id, p] is the array containing the cached four-point integrals in
LoopTools, where id is the coefficient label and n is the index returned
by Dget."

E0i::usage =
"E0i[id, p1, p2, p3, p4, p5, p1p2, p2p3, p3p4, p4p5, p5p1, m1, m2, m3,
m4, m5] is the generic five-point one-loop integral which includes both
scalar and tensor coefficients, specified by id.  For example,
E0i[ee0, ...] is the scalar function E_0, E0i[ee1244, ...] the tensor
function E_{1244} etc.  p1...p5 are the external momenta squared,
p1p2...p5p1 are the squares of external momenta (1+2)...(5+1),
respectively, and m1...m5 are the masses squared."

Eget::usage =
"Eget[p1, p2, p3, p4, p5, p1p2, p2p3, p3p4, p4p5, p5p1, m1, m2, m3, m4,
m5] computes all five-point coefficients in LoopTools."

Eval::usage =
"Eval[id, n] is the array containing the cached five-point integrals in
LoopTools, where id is the coefficient label and n is the index returned
by Eget."

F0i::usage =
"F0i[id, p1, p2, p3, p4, p5, p6, p1p2, p2p3, p3p4, p4p5, p5p6, p6p1,
p1p2p3, p2p3p4, p3p4p5, m1, m2, m3, m4, m5, m6] is the generic six-point
one-loop integral which includes both scalar and tensor coefficients,
specified by id.  For example, F0i[ff0, ...] is the scalar function F_0,
F0i[ff1244, ...] the tensor function F_{1244} etc.  p1...p6 are the
external momenta squared, p1p2...p6p1 are the squares of external
momenta (1+2)...(6+1), respectively, p1p2p3...p3p4p5 are the external
momenta (1+2+3)...(3+4+5) squared, and m1...m6 are the masses squared."

Fget::usage =
"Fget[p1, p2, p3, p4, p5, p6, p1p2, p2p3, p3p4, p4p5, p5p6, p6p1,
p1p2p3, p2p3p4, p3p4p5, m1, m2, m3, m4, m5, m6] computes all six-point
coefficients in a future version of LoopTools."

Fval::usage =
"Fval[id, n] is the array containing the cached six-point integrals in a
future version of LoopTools, where id is the coefficient label and n is
the index returned by Fget."


(* compatibility functions *)

B0::usage =
"B0[p, m1, m2] is the two-point one-loop scalar integral."

B1::usage =
"B1[p, m1, m2] is the coefficient of k_mu in the two-point one-loop
tensor integral B_mu."

B00::usage =
"B00[p, m1, m2] is the coefficient of g_{mu nu} in the two-point
one-loop tensor integral B_{mu nu}."

B11::usage =
"B11[p, m1, m2] is the coefficient of k_mu k_nu in the two-point
one-loop tensor integral B_{mu nu}."

B001::usage =
"B001[p, m1, m2] is the coefficient of g_{mu nu} k_rho in the two-point
one-loop tensor integral B_{mu nu rho}."

B111::usage =
"B111[p, m1, m2] is the coefficient of k_mu k_nu k_rho in the two-point
one-loop tensor integral B_{mu nu rho}."

DB0::usage =
"DB0[p, m1, m2] is the derivative of B0[p, m1, m2] with respect to p."

DB1::usage =
"DB1[p, m1, m2] is the derivative of B1[p, m1, m2] with respect to p."

DB00::usage =
"DB00[p, m1, m2] is the derivative of B00[p, m1, m2] with respect to p."

DB11::usage =
"DB11[p, m1, m2] is the derivative of B11[p, m1, m2] with respect to p."

C0::usage =
"C0[p1, p2, p1p2, m1, m2, m3] is the three-point scalar one-loop
integral."

D0::usage =
"D0[p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4] is the four-point scalar
one-loop integral."

E0::usage =
"E0[p1, p2, p3, p4, p5, p1p2, p2p3, p3p4, p4p5, p5p1, m1, m2, m3, m4,
m5] is the five-point scalar one-loop integral."

F0::usage =
"F0[p1, p2, p3, p4, p5, p6, p1p2, p2p3, p3p4, p4p5, p5p6, p6p1, p1p2p3,
p2p3p4, p3p4p5, m1, m2, m3, m4, m5, m6] is the six-point scalar one-loop
integral."

PaVe::usage =
"PaVe[ind, {pi}, {mi}] is the generalized Passarino-Veltman function
used by FeynCalc.  It is converted to B0i, C0i, D0i, E0i, or F0i in
FormCalc."

ToOldBRules::usage =
"ToOldBRules is a list of rules for converting two-point functions to
the old (LoopTools 2.1) conventions."

ToNewBRules::usage =
"ToNewBRules is a list of rules for converting two-point functions to
the new (LoopTools 2.2) conventions."

EndPackage[]


(* symbols from CutTools *)

BeginPackage["CutTools`"]

Acut::usage =
"Acut[rank, num, numtilde, m] is the one-point CutTools integral with
numerator functions num (4-dimensional part) and numtilde
((D-4)-dimensional part), where num contains rank powers of the
integration momentum.  m is the mass squared."

Bcut::usage =
"Bcut[rank, num, numtilde, p, m1, m2] is the two-point CutTools
integral with numerator functions num (4-dimensional part) and
numtilde ((D-4)-dimensional part), where num contains rank powers of
the integration momentum.  p is the external momentum and m1 and m2
are the masses squared."

Ccut::usage =
"Ccut[rank, num, numtilde, p1, p2, m1, m2, m3] is the three-point
CutTools integral with numerator functions  num (4-dimensional part)
and numtilde ((D-4)-dimensional part), where num contains rank
powers of the integration momentum.  p1 and p2 are the external
momenta and m1...m3 are the masses squared."

Dcut::usage =
"Dcut[rank, num, numtilde, p1, p2, p3, m1, m2, m3, m4] is the
four-point CutTools integral with numerator functions num
(4-dimensional part) and numtilde ((D-4)-dimensional part), where
num contains rank powers of the integration momentum.  p1...p3
are the external momenta and m1...m4 are the masses squared."

Ecut::usage =
"Ecut[num, p1, p2, p3, p4, m1, m2, m3, m4, m5] is the five-point
CutTools integral with numerator functions num (4-dimensional
part) and numtilde ((D-4)-dimensional part), where num contains
rank powers of the integration momentum.  p1...p4 are the
external momenta and m1...m5 are the masses squared."

Fcut::usage =
"Fcut[num, p1, p2, p3, p4, p5, m1, m2, m3, m4, m5, m6] is the
six-point CutTools integral with numerator functions num
(4-dimensional part) and numtilde ((D-4)-dimensional part), where
num contains rank powers of the integration momentum.  p1...p5
are the external momenta and m1...m6 are the masses squared."

EndPackage[]


BeginPackage["Form`"]

Attributes[cutM] = Attributes[qfM] = Attributes[GA] = {HoldAll}

Attributes[addM] = Attributes[mulM] = {Flat, Orderless}

{ q1, GA, intM, paveM, abbM, fmeM, sunM, d$$, e$$, i$$, dummy$$,
  dirM, sM, powM, addM, mulM, cutM, numM, qfM, qcM, helM }

EndPackage[]


(* symbols from the model files live in Global` *)

{ DiracMatrix, DiracSlash, ChiralityProjector,
  DiracSpinor, MajoranaSpinor, DiracObject,
  PolarizationVector, PolarizationTensor,
  MetricTensor, LeviCivita, FourVector, ScalarProduct,
  Lorentz, Lorentz4, EpsilonScalar,
  SUNT, SUNF, SUNTSum, SUNEps, Colour, Gluon }

SelfEnergy::usage =
"SelfEnergy[from -> to, mass] calculates the self-energy with incoming
particle from and outgoing particle to, taken at k^2 = mass^2. 
SelfEnergy[f] calculates the self-energy of particle f on its mass
shell.  SelfEnergy[..., opt] specifies InsertFields options to be used
in the computation."

DSelfEnergy::usage =
"DSelfEnergy[from -> to, mass] calculates the derivative with respect to
k^2 of the self-energy with incoming particle from and outgoing particle
to, taken at k^2 = mass^2.  DSelfEnergy[f] calculates the self-energy of
particle f on its mass shell.  DSelfEnergy[..., opt] specifies
InsertFields options to be used in the computation."

ReTilde::usage =
"ReTilde[expr] takes the real part of loop integrals occurring in expr."

ImTilde::usage =
"ImTilde[expr] takes the imaginary part of loop integrals occurring in
expr."

LVectorCoeff::usage =
"LVectorCoeff[expr] returns the coefficient of DiracChain[6, k[1]]
(= k1slash omega_-) in expr."

RVectorCoeff::usage =
"RVectorCoeff[expr] returns the coefficient of DiracChain[7, k[1]]
(= k1slash omega_+) in expr."

LScalarCoeff::usage =
"LScalarCoeff[expr] returns the coefficient of DiracChain[7] (= omega_-)
in expr."

RScalarCoeff::usage =
"RScalarCoeff[expr] returns the coefficient of DiracChain[6] (= omega_+)
in expr."

SEPart::usage =
"SEPart[p, se] returns part p of self-energy se.  It is applied during
the calculation of a renormalization constant, where p is one of
LVectorCoeff, RVectorCoeff, LScalarCoeff, RScalarCoeff for fermion
self-energies, and Identity otherwise."

MassRC::usage =
"MassRC[f] computes the one-loop mass renormalization constant of field
f.  MassRC[f1, f2] computes the one-loop mass renormalization constant
for the f1-f2 transition.  For fermions the output is a list
{left-handed RC, right-handed RC}.  MassRC[..., opt] specifies
InsertFields options to be used in the computation."

FieldRC::usage =
"FieldRC[f] computes the one-loop field renormalization constant of
field f.  FieldRC[f1, f2] computes the one-loop field renormalization
constant for the f1-f2 transition.  FieldRC[f1, f2, c] subtracts c from
the self-energy entering into the calculation.  For fermions the output
is a list {left-handed RC, right-handed RC}.  FieldRC[..., opt]
specifies InsertFields options to be used in the computation."

TadpoleRC::usage =
"TadpoleRC[f] computes the one-loop tadpole renormalization constant of
field f.  TadpoleRC[..., opt] specifies InsertFields options to be used
in the computation."

WidthRC::usage =
"WidthRC[f] computes the one-loop width of field f.  WidthRC[..., opt]
specifies InsertFields options to be used in the computation."


BeginPackage["FormCalc`",
  {"FeynArts`", "LoopTools`", "CutTools`", "Form`", "Global`",
   "Utilities`FilterOptions`"}]

(* some internal symbols must be visible for FORM/ReadForm *)

{ SUNSum, ReadForm, ReadFormClear, ReadFormDebug,
  FormEval, FormEvalDecl, FormExpr }

(* some internal symbols made visible for debugging *)

{ FormKins, KinFunc, InvSum, PairRules, LastAmps, FormSetup }


(* symbols appearing in the output *)

Amp::usage =
"Amp[proc][expr1, expr2, ...] is the result of the calculation of
diagrams of the process proc.  The result is divided into parts expr1,
expr2, ..., such that index sums (marked by SumOver) apply to the whole
of each part."

Den::usage =
"Den[p, m] stands for 1/(p - m).  Note that in contrast to
PropagatorDenominator, p and m are the momentum and mass *squared*.
Den[p, m, d] is the denominator raised to the power d."

Num::usage =
"Num[expr] contains the numerator of a loop integral as a function of
the loop momentum.  This representation is used when calculating loop
integrals by the CutTools package."

DiracChain::usage =
"DiracChain[objs] is a chain of Dirac matrices contracted with the given
objects.  The integers 1, 5, 6, and 7 appearing as first argument denote
1, gamma_5, (1 + gamma_5)/2, and (1 - gamma_5)/2, respectively."

WeylChain::usage =
"WeylChain[objs] is a chain of sigma matrices contracted with the given
objects.  The integers 6, 7 respectively denote upper and lower indices
at the given position, and -1 stands for epsilon, the spinor metric."

Spinor::usage =
"Spinor[p, m, s] is a spinor with momentum p and mass m, i.e. a solution
of the Dirac equation (pslash + s m) Spinor[p, m, s] = 0.  On screen,
particle spinors (s = 1) are printed as u[p, m], antiparticle spinors
(s = -1) as v[p, m].  Inside a WeylChain, Spinor denotes a 2-dimensional
Weyl spinor.  Whether it corresponds to the upper or lower half of the
4-dimensional Dirac spinor is determined by the index convention of the
WeylChain (fixed by arguments 6 or 7), propagated to the position of the
spinor."

DottedSpinor::usage =
"DottedSpinor[p, m, s] denotes the 2-dimensional conjugated spinor
corresponding to Spinor[p, m, s]."

e::usage =
"e[i] is the ith polarization vector."

ec::usage =
"ec[i] is the conjugate of the ith polarization vector."

z::usage =
"z[i] is the ith polarization vector in D - 4 dimensions."

zc::usage =
"zc[i] is the conjugate of the ith polarization vector in D - 4
dimensions."

eT::usage =
"eT[i] is the ith polarization tensor."

eTc::usage =
"eTc[i] is the conjugate of the ith polarization tensor."

k::usage =
"k[i] is the ith momentum."

SUNN::usage =
"SUNN specifies the N in SU(N), i.e. the number of colours."

S::usage =
"S is the Mandelstam variable s.  If k1 and k2 are the incoming momenta,
S = (k1 + k2)^2."

T::usage =
"T is the Mandelstam variable t.  If k1 denotes the first incoming and
k3 the first outgoing momentum, T = (k1 - k3)^2."

U::usage =
"U is the Mandelstam variable u.  If k2 denotes the first incoming and
k3 the second outgoing momentum, U = (k2 - k3)^2."

Dminus4::usage =
"Dminus4 represents the difference D - 4."


(* DeclareProcess, CalcFeynAmp and their options *)

DeclareProcess::usage =
"DeclareProcess[amps] sets up internal definitions for the computation
of the amplitudes amps."

OnShell::usage =
"OnShell is an option of DeclareProcess.  It specifies whether FORM
should put the external particles on their mass shell, i.e. apply
ki^2 = mi^2."

Invariants::usage =
"Invariants is an option of DeclareProcess.  It specifies whether FORM
should introduce kinematical invariants, like the Mandelstam variables
for a 2 -> 2 process."

Transverse::usage =
"Transverse is an option of DeclareProcess.  It specifies whether FORM
should apply the transversality relations for polarization vectors
(ei.ki = 0)."

Normalized::usage =
"Normalized is an option of DeclareProcess.  It specifies whether FORM
should apply the normalization of polarization vectors (ei.ei^* = -1)."

InvSimplify::usage =
"InvSimplify is an option of DeclareProcess.  It specifies whether FORM
should try to simplify combinations of invariants as much as possible."

Antisymmetrize::usage =
"Antisymmetrize is an option of DeclareProcess.  It specifies whether
Dirac chains are antisymmetrized."

MomElim::usage =
"MomElim is an option of DeclareProcess.  It controls in which way
momentum conservation is used to eliminate momenta.  False performs no
elimination, an integer between 1 and the number of legs substitutes
the specified momentum in favour of the others, and Automatic tries all
substitutions and chooses the one resulting in the fewest terms."

FormAmp::usage =
"FormAmp[proc][amps] contains a preprocessed form of the FeynArts
amplitudes amps for process proc, to be calculated by CalcFeynAmp."

CalcFeynAmp::usage =
"CalcFeynAmp[amps] calculates the Feynman amplitudes given in amps.  The
resulting expression is broken up into categories which are returned in
an Amp object."

CalcLevel::usage =
"CalcLevel is an option of CalcFeynAmp.  It specifies the level (Classes
or Particles) at which to calculate the amplitudes. Automatic takes
Classes level, if available, otherwise Particles."

Dimension::usage =
"Dimension is an option of CalcFeynAmp.  It specifies the space-time
dimension in which to perform the calculation and can take the values D,
where dimensional regularization is used, and 4, where constrained
differential renormalization is used.  The latter method is equivalent
to dimensional reduction at the one-loop level. 
Experimental: Dimension -> 0 retains the Dminus4 terms, i.e. does not
emit local (rational) terms."

FermionChains::usage =
"FermionChains is an option of CalcFeynAmp.  It can take the three
values Chiral, VA, and Weyl, which specify how fermion chains are
returned by CalcFeynAmp.  Chiral and VA both select (4-dimensional)
Dirac chains, where the chiral decomposition is taken for Chiral and the
vector/axial-vector decomposition for VA.  Weyl selects (2-dimensional)
Weyl chains."

Chiral::usage =
"Chiral is a possible value for the FermionChains option of CalcFeynAmp. 
It instructs CalcFeynAmp to return fermion chains as left- and
right-handed (4-dimensional) Dirac chains."

VA::usage =
"VA is a possible value for the FermionChains option of CalcFeynAmp. 
It instructs CalcFeynAmp to return fermion chains as vector and
axial-vector parts of (4-dimensional) Dirac chains."

Weyl::usage =
"Weyl is a possible value for the FermionChains option of CalcFeynAmp. 
It instructs CalcFeynAmp to return fermion chains as (2-dimensional)
Weyl chains."

FermionOrder::usage =
"FermionOrder is an option of CalcFeynAmp.  It determines the ordering
of Dirac spinor chains in the output, i.e. requires FermionChains ->
Chiral or VA.  Possible values are None, Fierz, Automatic, Colour, or
an explicit ordering, e.g. {2, 1, 4, 3}.  None applies no reordering. 
Fierz applies the Fierz identities twice, thus simplifying the chains
but keeping the original order.  Colour applies the ordering of the
external colour indices (after simplification) to the spinors. 
Automatic chooses a lexicographical ordering."

Fierz::usage =
"Fierz is a possible value for the FermionOrder option of CalcFeynAmp. 
It instructs CalcFeynAmp to apply the Fierz identities twice, thus
simplifying the chains but keeping the original spinor order."

Colour::usage =
"Colour is a possible value for the FermionOrder option of CalcFeynAmp. 
It instructs CalcFeynAmp to bring the spinors into the same order as the
external colour indices, i.e. \"Fermion flow follows colour flow\"."

InsertionPolicy::usage =
"InsertionPolicy is an option of CalcFeynAmp.  Begin specifies that the
insertions shall be applied at the beginning of the FORM code (this
ensures that all loop integrals are fully symmetrized).  Default applies
them after simplifying the generic amplitudes (this is fastest).  A
positive integer does the same, except that insertions with a LeafCount
larger than that integer are inserted only after the amplitude comes
back from FORM (this is a workaround for the rare cases where the FORM
code aborts due to very long insertions)."

CutTools::usage =
"CutTools is an option of CalcFeynAmp.  It can take the three values
False, True, and Rational.  False chooses the regular Passarino-Veltman
reduction with LoopTools' tensor-coefficient functions.  True and
Rational select the CutTools functions instead, which receive the
integral's numerator as a function of the loop momentum.  Rational
computes the rational terms analytically while True leaves their
computation to CutTools."

NoExpand::usage =
"NoExpand is an option of CalcFeynAmp.  NoExpand -> {sym1, sym2, ...}
specifies that sums containing any of sym1, sym2, ... are not expanded
during the FORM calculation."

NoBracket::usage =
"NoBracket is an option of CalcFeynAmp.  NoExpand -> {sym1, sym2, ...} 
specifies that sym1, sym2, ... are not collected inside a multiplication 
bracket during the FORM calculation."

AbbrScale::usage =
"AbbrScale is an option of CalcFeynAmp.  The automatically introduced
abbreviations are scaled by the square root of the provided factor for
every momentum they contain.  Thus AbbrScale -> S makes the
abbreviations dimensionless, which can be of advantage in some
applications, e.g. the treatment of resonances."

EditCode::usage =
"EditCode is a debugging option of CalcFeynAmp, HelicityME, and
PolarizationSum.  It edits the temporary file passed to FORM using the
editor command in $Editor."

RetainFile::usage =
"RetainFile is a debugging option of CalcFeynAmp, HelicityME, and
PolarizationSum.  When set to True, the temporary file containing the
FORM code is not removed after running FORM."


(* abbreviationing-related functions *)

Abbr::usage =
"Abbr[] returns a list of all abbreviations introduced so far.
Abbr[patt] returns a list of all abbreviations including the pattern
patt.  Patterns prefixed by ! (Not) are excluded."

GenericList::usage =
"GenericList[] returns a list of the substitutions made for the
computation of generic amplitudes."

ClearProcess::usage =
"ClearProcess[] is necessary to clear internal definitions before
calculating a process with a different kinematical set-up."

ZapFunction::usage =
"ZapFunction is an option of ClearProcess and ClearSubexpr and
determines the function used to clear the definitions of symbols
introduced for abbreviations."

RegisterAbbr::usage =
"RegisterAbbr[abbr] registers a list of abbreviations such that
future invocations of CalcFeynAmp will make use of them.  Note that
abbreviations introduced for different processes are in general not
compatible."

Abbreviate::usage =
"Abbreviate[expr, lev] introduces abbreviations for all subexpressions
in expr (currently only sums are considered), starting at level lev. 
The optional parameter lev defaults to 2. 
Abbreviate[expr, patt] introduces abbreviations for all subexpressions
of expr free of patt."

Deny::usage =
"Deny is an option of Abbreviate.  It specifies items which must not be
included in abbreviations."

Preprocess::usage =
"Preprocess is an option of Abbreviate.  It specifies a function to be
applied to all subexpressions before introducing abbreviations for
them."

$SubPrefix::usage =
"$SubPrefix specifies the prefix for subexpressions introduced by
Abbreviate, i.e. the Sub in Sub123."

Subexpr::usage =
"Subexpr[] returns a list of all subexpressions introduced by
Abbreviate."

ClearSubexpr::usage =
"ClearSubexpr[] clears the internal definitions of the subexpressions
introduced by Abbreviate."

RegisterSubexpr::usage =
"RegisterSubexpr[subexpr] registers a list of subexpressions such that
future invocations of Abbreviate will make use of them."

OptimizeAbbr::usage =
"OptimizeAbbr[abbr] optimizes the abbreviations in abbr by eliminating
common subexpressions."

$OptPrefix::usage =
"$OptPrefix specifies the prefix for additional abbreviations introduced
by OptimizeAbbr, i.e. the Opt in Opt123."

SubstSimpleAbbr::usage =
"SubstSimpleAbbr[{expr, abbr}] removes `simple' abbreviations from abbr
and substitutes them back into expr. 
SubstSimpleAbbr[FortranExpr[var, tmpvar, expr]] removes `simple'
abbreviations from expr and deletes them from the variable lists var
and tmpvar."

ExtractInt::usage =
"ExtractInt[expr] extracts the loop integrals from expr for evaluation
with LoopTools/CutTools.  It output is a list {cint, iint, abbrexpr},
where cint is the list of complex-valued loop functions, iint is the
list of integer-valued loop functions, and abbrexpr is expr with the
loop integrals substituted by the identifiers in cint and iint."

Pair::usage =
"Pair[a, b] represents the contraction of the two four-vectors or
Lorentz indices a and b."

Eps::usage =
"Eps[a, b, c, d] represents -I times the antisymmetric Levi-Civita
tensor epsilon_{abcd}.  The sign convention is epsilon^{0123} = +1."

ToSymbol::usage =
"ToSymbol[s...] concatenates its arguments into a new symbol."

ToArray::usage =
"ToArray[s] turns the symbol s into an array reference by taking it
apart into letters and digits, e.g. Var1234 becomes Var[1234]. 
ToArray[expr, s1, s2, ...] turns all occurrences of the symbols s1NNN,
s2NNN, etc. in expr into s1[NNN], s2[NNN], etc."

Renumber::usage =
"Renumber[expr, var1, var2, ...] renumbers all var1[n], var2[n], ... in
expr."

MaxDims::usage =
"MaxDims[args] returns a list of all distinct functions in args with the
highest indices that appear, e.g. MaxDims[foo[1, 2], foo[2, 1]] returns
{foo[2, 2]}."

Keep::usage =
"Keep[expr, name, path] loads path/name.m if it exists, otherwise it
evaluates expr and stores it (together with the Abbr[] and Subexpr[])
in that file.  path is optional and defaults to $KeepDir. 
Keep[lhs = rhs] is short for lhs = Keep[rhs, \"lhs\"]."

$KeepDir::usage =
"$KeepDir specifies the default directory for storing intermediate
expressions with Keep."


(* miscellaneous functions *)

FormSimplify::usage =
"FormSimplify is a function applied to parts of the FORM output for
simplification."

DenCollect::usage =
"DenCollect[expr] collects terms in expr whose denominators are
identical up to a numerical constant.  DenCollect[expr, wrap] applies
wrap to the collected numerators."

Pool::usage =
"Pool[expr] combines terms with common factors.  Unlike Factor, it looks
at the terms pairwise and can thus do a b + a c + d -> a (b + c) + d
fast.  Unlike Simplify, it does not modify b and c.  Pool[expr, wrap]
applies wrap to the (b + c) part."

ApplyUnitarity::usage =
"ApplyUnitarity[expr, mat, d, simp] simplifies expr by exploiting the
unitarity of the d-dimensional matrix mat.  The optional argument simp
specifies the simplification function to use internally and defaults to
FullSimplify."

OnSize::usage =
"OnSize[n1, f1, n2, f2, ..., fdef][expr] returns f1[expr] if
LeafCount[expr] < n1, f2[expr] if LeafCount[expr] < n2, etc., and
fdef[expr] if the expression is still larger.  fdef can take the
special value Map which means that fdef[expr] recursively applies
the entire OnSize function to the parts of expr.  If omitted, fdef
defaults to Identity."

OffShell::usage =
"OffShell[amps, i -> mi, ...] returns the FeynAmpList amps with the mass
of the ith external particle set to mi.  This will in general take
particle i off its mass shell since now ki^2 = mi^2 is fulfilled with
the new value of mi."

Combine::usage =
"Combine[amp1, amp2, ...] combines the amplitudes amp1, amp2, ... which
can be either FeynAmpList or Amp objects, i.e. Combine works before and
after CalcFeynAmp."

ExpandSums::usage =
"ExpandSums[expr] turns all pieces of expr multiplied with SumOver
into an actual Sum.  ExpandSums[expr, h] uses h instead of Sum."

MultiplyDiagrams::usage =
"MultiplyDiagrams[func][amp] multiplies the diagrams in amp with the
factor returned by the function func.  The latter is invoked for each
diagram either as func[amplitude] (for a fully inserted diagram), or as
func[generic amplitude, insertion]."

TagDiagrams::usage =
"TagDiagrams[amp] tags each diagram in amp with an identifier of the
form Diagram[number], where number runs sequentially through the
diagrams at all levels.  This makes it possible to locate the
contribution of individual diagrams in the final CalcFeynAmp output. 
TagDiagrams[amp, tag] uses tag rather than Diagram."

Diagram::usage =
"Diagram[number] is the identifier used to tag a single diagram by
TagDiagrams."

DiagramType::usage =
"DiagramType[diag] returns the number of denominators not containing
the integration momentum."

FermionicQ::usage =
"FermionicQ[diag] gives True for a diagram containing fermions and
False otherwise."

IndexIf::usage =
"IndexIf[cond, a, b] is identical to If[cond, a, b], except that a
and b are not held unevaluated.  IndexIf[cond, a] is equivalent to
IndexIf[cond, a, 0] and IndexIf[cond1, a1, cond2, a2, ...] is
equivalent to IndexIf[cond1, a1, IndexIf[cond2, a2, ...]]."

IndexDiff::usage =
"IndexDiff[i, j] is the same as 1 - IndexDelta[i, j]."

ToIndexIf::usage =
"ToIndexIf[expr] converts all IndexDeltas and IndexDiffs in expr to
IndexIf, which will be written out as if-statements in the generated
Fortran code.  ToIndexIf[expr, patt] operates only on indices matching
patt.  If patt is a string, e.g. \"Glu*\", it is first expanded to all
matching symbols."

Neglect::usage =
"Neglect[sym] = 0 makes FORM replace sym = 0 except when it appears in
negative powers or in loop integrals."

Square::usage =
"Square[m] = m2 makes FORM replace all m^2 by m2."

NClear::usage =
"NClear[patt] clears the NValues of all symbols matching patt. 
NClear[] is equivalent to NClear[\"Global`*\"]."

ColourSimplify::usage =
"ColourSimplify[expr] simplifies the colour objects in expr.
ColourSimplify[plain, conj] simplifies the colour objects in
(plain conj^*)."

ColourGrouping::usage =
"ColourGrouping[tops] returns a list of parts of the inserted topologies
tops, grouped according to their colour structures."


(* FeynCalc compatibility functions *)

FeynCalcGet::usage =
"FeynCalcGet[mask] reads files produced with FeynCalc.  mask is taken
as input to the Mathematica function FileNames, so it might be
FeynCalcGet[\"file.m\"] or FeynCalcGet[\"*.m\"] or
FeynCalcGet[\"*.m\", \"~/feyncalcfiles\"]."

FeynCalcPut::usage =
"FeynCalcPut[expr, file] writes expr to file in FeynCalc format."


(* finiteness checks *)

UVDivergentPart::usage =
"UVDivergentPart[expr] returns expr with all loop integrals replaced
by their UV-divergent part.  The divergence itself is denoted by
Divergence."

UVSeries::usage =
"UVSeries[expr] expands expr into a Laurent series in Dminus4.
UVSeries[expr, pow] returns the coefficient of Dminus4^pow."

Divergence::usage =
"Divergence represents the dimensionally regularized divergence
2/(4 - D) of loop integrals.  It is used by the function
UVDivergentPart."


(* matrix elements *)

HelicityME::usage =
"HelicityME[plain, conj] calculates the helicity matrix elements for all
combinations of spinor chains that appear in the expression
(plain conj^*).  Terms of this kind arise in the calculation of the
squared matrix element, where typically plain is the one-loop result
and conj the tree-level expression.  The arguments do not necessarily
have to be amplitudes since they are only used to determine which spinor
chains to select from the abbreviations.  The symbol All can be used to
select all spinor chains currently defined in the abbreviations."

ColourME::usage =
"ColourME[plain, conj] calculates the colour matrix elements.  ColourME
is very similar to HelicityME, except that it computes the matrix
elements for SU(N) objects, not for spinor chains."

All::usage =
"All as an argument of HelicityME and ColourME indicates that all spinor
chains or SUNT objects currently defined in the abbreviations should be
used instead of just those appearing in the argument."

Source::usage =
"Source is an option of HelicityME, ColourME, and PolarizationSum. 
It specifies from where the abbreviations used to calculate the matrix
elements are taken."

Hel::usage =
"Hel[i] is the helicity of the ith external particle.  It can take the
values +1, 0, -1, where 0 stands for an unpolarized particle."

s::usage =
"s[i] is the ith helicity reference vector."

Mat::usage =
"Mat[Fi SUNi] is a matrix element in an amplitude, i.e. an amplitude is
a linear combination of Mat objects.\n Mat[Fi, Fj] appears in the
squared matrix element and stands for the product of the two arguments,
Fi Fj^*.  Such expressions are calculated by HelicityME and ColourME."

Lor::usage =
"Lor[i] is a contracted Lorentz index in a product of Dirac chains."

SquaredME::usage =
"SquaredME[plain, conj] returns the matrix element
(plain Conjugate[conj]).  This performs a nontrivial task only for
fermionic amplitudes: the product of two fermionic amplitudes\n
    M1 = a1 F1 + a2 F2 + ... and\n
    M2 = b1 F1 + b2 F2 + ... is returned as\n
    M1 M2^* = a1 b1^* Mat[F1, F1] + a2 b1^* Mat[F2, F1] + ...\n
The special case of plain === conj can be written as SquaredME[plain]
which is of course equivalent to SquaredME[plain, plain]."

RealQ::usage =
"RealQ[sym] is True if sym represents a real quantity which means in
particular that Conjugate[sym] = sym."

PolarizationSum::usage =
"PolarizationSum[expr] sums expr over the polarizations of external
gauge bosons.  It is assumed that expr is the squared amplitude into
which the helicity matrix elements have already been inserted. 
Alternatively, expr may also be given as an amplitude directly, in which
case PolarizationSum will first invoke SquaredME and HelicityME (with
Hel[_] = 0) to obtain the squared amplitude."

GaugeTerms::usage =
"GaugeTerms is an option of PolarizationSum.  With GaugeTerms -> False,
the gauge-dependent terms in the polarization sum, which should
eventually cancel in gauge-invariant subsets of diagrams, are omitted
from the beginning."

eta::usage =
"eta[i] is a vector that defines a particular gauge via Pair[eta[i],
e[i]] = 0.  It is introduced by PolarizationSum for massless particles,
where the sum over e[i][mu] ec[i][nu] is gauge dependent.  Only for
gauge-invariant subsets of diagrams should the dependence on eta[i]
cancel.  eta obeys Pair[eta[i], k[i]] != 0 and Pair[eta[i], e[i]] = 0."


(* writing out Fortran code *)

SetupCodeDir::usage =
"SetupCodeDir[dir] installs the driver programs necessary to compile the
Fortran code generated by WriteSquaredME and WriteRenConst in the
directory dir.  Customized versions of the drivers are taken from the
directory pointed to by the Drivers option and take precedence over the
default versions from $DriversDir.  Drivers already in dir are not
overwritten."

Drivers::usage =
"Drivers is an option of SetupCodeDir.  Drivers points to a directory
containing customized versions of the driver programs necessary for
compiling the generated Fortran code.  This directory need not contain
all driver programs: files not contained therein are taken from the
default directory $DriversDir."

WriteSquaredME::usage =
"WriteSquaredME[tree, loop, me, abbr, ..., dir] writes out Fortran code
to compute the squared matrix element for a process whose tree-level and
one-loop contributions are given in the first and second argument,
respectively.  All further arguments except the last specify the
necessary matrix elements and abbreviations.  The last argument dir
finally gives the path to write the generated code to."

ExtraRules::usage =
"ExtraRules is an option of WriteSquaredME.  Rules given here will be
applied before the loop integrals are abbreviated."

TreeSquare::usage =
"TreeSquare is an option of WriteSquaredME.  It specifies whether to
add the square of the tree-level amplitude, |M_0|^2, to the result 
in the SquaredME subroutine."

LoopSquare::usage =
"LoopSquare is an option of WriteSquaredME.  It specifies whether to
add the square of the 1-loop amplitude, |M_1|^2, to the result in the
SquaredME subroutine.  This term is of order alpha^2 with respect to the
tree-level contribution, |M_0|^2.  Usually one takes into account only
the interference term, 2 Re M_0^* M_1, which is of order alpha."

Folder::usage =
"Folder is an option of WriteSquaredME and WriteRenConst.  It specifies
the folder into which the generated files are written."

FilePrefix::usage =
"FilePrefix is an option of WriteSquaredME and WriteRenConst.  It
specifies a string to be prepended to the filenames of the generated
code."

SymbolPrefix::usage =
"SymbolPrefix is an option of WriteSquaredME and WriteRenConst.  It
specifies a string which is prepended to externally visible symbols in
the generated Fortran code to prevent collision of names when several
processes are linked together."

FileIncludes::usage =
"FileIncludes is an option of WriteSquaredME and WriteRenConst.  It
specifies per-file #include statements (or other declarations)."

SubroutineIncludes::usage =
"SubroutineIncludes is an option of WriteSquaredME and WriteRenConst. 
It specifies per-subroutine #include statements (or other
declarations)."

FileHeader::usage =
"FileHeader is an option of WriteSquaredME and WriteRenConst and
specifies a file header.  This string may contain %f, %d, and %t, which
are substituted at file creation by file name, description, and time
stamp, respectively."


(* renormalization constants *)

FindRenConst::usage =
"FindRenConst[expr] returns a list of all renormalization constants
found in expr including those needed to compute the former."

CalcRenConst::usage =
"CalcRenConst[expr] calculates the renormalization constants appearing
in expr."

WriteRenConst::usage =
"WriteRenConst[expr, dir] calculates the renormalization constants
appearing in expr and generates a Fortran program from the results. 
The resulting files (the Fortran program itself and the corresponding
declarations) are written to the directory dir.  The names of the files
are determined by the RenConstFile option."

InsertFieldsHook::usage =
"InsertFieldsHook[tops, proc] is the function called by SelfEnergy and
DSelfEnergy to insert fields into the topologies tops for the process
proc.  It is normally equivalent to InsertFields, but may be redefined
to change the diagram content of certain self-energies."

CreateFeynAmpHook::usage =
"CreateFeynAmpHook[diags, opt] is the function called by SelfEnergy and
DSelfEnergy to create the amplitudes for diagrams diags.  It is normally
equivalent to CreateFeynAmp, but may be redefined to modify the
amplitudes."

ClearSE::usage =
"ClearSE[] clears the internal definitions of already calculated
self-energies."

OptPaint::usage =
"OptPaint[ins, True] invokes Paint[ins]. 
OptPaint[ins, pre] invokes Paint[ins], saving the graphics in files
with names constructed from ProcessName[ins], prefixed with the
string pre (which may include a directory). 
OptPaint[ins, pre, suf] further appends suf to the file name. 
OptPaint[ins] executes OptPaint[ins, $PaintSE]."

$PaintSE::usage =
"$PaintSE determines whether SelfEnergy paints the diagrams it generates
to compute the self-energies.  $PaintSE can be True, False, or a string
which indicates that the output should be saved in a PostScript file
instead of being displayed on screen, and is prepended to the filename."

$LongitudinalSE::usage =
"$LongitudinalSE specifies that the longitudinal rather than the
transverse part of the vector-boson self-energies is taken in
SelfEnergy and DSelfEnergy."


(* low-level Fortran output functions *)

ToList::usage =
"ToList[expr] returns a list of summands of expr."

MkDir::usage =
"MkDir[\"dir1\", \"dir2\", ...] makes sure the directory dir1/dir2/...
exists, creating the individual subdirectories dir1, dir2, ... as
necessary."

ToForm::usage =
"ToForm[expr] returns the FORM form of expr as a string."

OpenForm::usage =
"OpenForm[file] opens file for writing in FORM format. 
OpenForm[] opens a temporary FORM file with a unique name for writing."

ToFortran::usage =
"ToFortran[expr] returns the Fortran form of expr as a string."

OpenFortran::usage =
"OpenFortran[file] opens file for writing in Fortran format."

TimeStamp::usage =
"TimeStamp[] returns a string with the current date and time."

BlockSplit::usage =
"BlockSplit[var -> expr] tries to split the calculation of expr into
subexpressions each of which has a leaf count less than $BlockSize."

FileSplit::usage =
"FileSplit[exprlist, mod, writemod, writeall] splits exprlist into
batches with leaf count less than $FileSize.  If there is only one
batch, writemod[batch, mod] is invoked to write it to file.  Otherwise,
writemod[batch, modN] is invoked on each batch, where modN is mod
suffixed by a running number, and in the end writeall[mod, res] is
called, where res is the list of writemod return values.  The optional
writeall function can be used e.g. to write out a master subroutine
which invokes the individual modules."

FortranNames::usage =
"FortranNames[base, ind] constructs two names out of base and ind,
where the latter are typically indices.  The first is the direct
concatenation of the elements, separated by underscores, and is for
use in file names and similar uncritical places.  The second is
limited to a maximum of $MaxFortranName characters and should be
used for Fortran symbols to comply with compiler limits.  Truncation
is done by leaving out delimiting underscores and if that is not
enough, contracting index names down to 2 characters."

$MaxFortranName::usage =
"$MaxFortranName specifies the maximum length of a symbol name in
Fortran."

RuleAdd::usage =
"RuleAdd[var, expr] is equivalent to var -> var + expr."

PrepareExpr::usage =
"PrepareExpr[{var1 -> expr1, var2 -> expr2, ...}] prepares a list of
variable assignments for write-out to a Fortran file. Expressions with a
leaf count larger than $BlockSize are split into several pieces, as in\n
\tvar = part1\n\
\tvar = var + part2\n\
\t...\n
thereby possibly introducing temporary variables for subexpressions. 
The output is a FortranExpr[vars, tmpvars, exprlist] object, where vars
are the original and tmpvars the temporary variables introduced by
PrepareExpr."

WriteExpr::usage =
"WriteExpr[file, exprlist] writes a list of variable assignments in
Fortran format to file.  The exprlist can either be a FortranExpr object
or a list of expressions of the form {var1 -> expr1, var2 -> expr2,
...}, which is first converted to a FortranExpr object using
PrepareExpr.  WriteExpr returns a list of the subexpressions that were
actually written."

HornerStyle::usage =
"HornerStyle is an option of WriteExpr.  It specifies whether
expressions are arranged in Horner form before writing them out as
Fortran code."

FinalCollect::usage =
"FinalCollect is an option of WriteExpr.  It specifies whether common
factors are collected in the final Fortran expression, just before
write-out."

Type::usage =
"Type is an option of WriteExpr.  If a string is given, e.g. Type ->
\"double precision\", WriteExpr writes out declarations of that type for
the given expressions.  Otherwise no declarations are produced."

TmpType::usage =
"TmpType is an option of WriteExpr.  It is the counterpart of Type for
the temporary variables.  TmpType -> Type uses the settings of the Type
option."

IndexType::usage =
"IndexType is an option of WriteExpr.  It is the counterpart of Type
for do-loop indices.  IndexType -> Type uses the settings of the Type
option."

RealArgs::usage =
"RealArgs is an option of WriteExpr.  It specifies a list of functions
whose numerical arguments must be of a guaranteed type (usually real). 
For example, if foo expects a real argument, it must be invoked as
foo(0D0), not foo(0) in Fortran.
RealArgs[foo] := ... defines the actual conversion for foo.  The
default is to map NArgs over all arguments, which turns the integers
into reals."

NArgs::usage =
"NArgs[args] returns args with integers turned into reals.  Note that
NArgs is not quite the same as N: while it changes 1 to 1., it leaves
m[1] intact so that array indices remain integers."

Newline::usage =
"Newline is an option of WriteExpr.  It specifies a string to be printed
after each Fortran statement."

Optimize::usage =
"Optimize is an option of PrepareExpr.  With Optimize -> True, variables
are introduced for subexpressions which are used more than once."

Expensive::usage =
"Expensive is an option of PrepareExpr.  It specifies patterns of objects
whose evaluation is expensive in terms of CPU time and which should be
hoisted from inner do-loops if possible."

MinLeafCount::usage =
"MinLeafCount is an option of PrepareExpr and Abbreviate.  It specifies
the minimum LeafCount a common subexpression must have in order that a
variable is introduced for it."

DebugLines::usage =
"DebugLines is an option of PrepareExpr.  It specifies whether debugging
statements are written out for each variable.  If instead of True a
string is given, the debugging messages are prefixed by this string. 
The actual debugging statements are constructed from the items in
$DebugCmd."

FinalTouch::usage =
"FinalTouch is an option of PrepareExpr.  It specifies a function which
is applied to each final subexpression, such as will then be written out
to the Fortran file."

FortranExpr::usage =
"FortranExpr[vars, tmpvars, exprlist] is the output of PrepareExpr and
contains a list of expressions ready to be written to a Fortran file,
where vars are the original variables and tmpvars are temporary
variables introduced in order to shrink individual expressions to a size
small enough for Fortran."

DebugLine::usage =
"DebugLine[var] emits a debugging statement (print-out of variable var)
when written to a Fortran file with WriteExpr.  DebugLine[var, tag]
prefixes the debugging message with the string \"tag\"."

$DebugCmd::usage =
"$DebugCmd is a list of three strings used to construct debugging
statements in Fortran.  The first string is the name of a preprocessor
variable to shield the debugging statement, the second the actual
command used to print debugging output, and the third a possible
text to be printed after the debugging line (e.g. to restore
terminal properties)."

SplitSums::usage =
"SplitSums[expr] splits expr into a list of expressions such that index
sums (marked by SumOver) always apply to the whole of each part. 
SplitSums[expr, wrap] applies wrap to the coefficients of the SumOver."

ToDoLoops::usage =
"ToDoLoops[list, ifunc] splits list into patches which must be summed
over the same set of indices.  ifunc is an optional argument, where
ifunc[expr] must return the indices occurring in expr."

DoLoop::usage =
"DoLoop[expr, ind] is a symbol introduced by ToDoLoops indicating that
expr is to be summed over the indices ind."

Dim::usage =
"Dim[i] returns the highest value the index i takes on, as determined
from the amplitude currently being processed."

MoveDepsRight::usage =
"MoveDepsRight[r1, ..., rn] shuffles variable definitions (var -> value)
among the lists of rules ri such that the definitions in each list do
not depend on definitions in ri further to the left.  For example,
MoveDepsRight[{a -> b}, {b -> 5}] produces {{}, {b -> 5, a -> b}}, i.e.
it moves a -> b to the right list because that depends on b."

MoveDepsLeft::usage =
"MoveDepsLeft[r1, ..., rn] shuffles variable definitions (var -> value)
among the lists of rules ri such that the definitions in each list do
not depend on definitions in ri further to the right.  For example,
MoveDepsLeft[{a -> b}, {b -> 5}] produces {{b -> 5, a -> b}, {}}, i.e.
it moves b -> 5 to the left list because that depends on b."

OnePassOrder::usage =
"OnePassOrder[r] orders a list of interdependent rules such that the
definition of each item (item -> ...) comes before its use in the
right-hand sides of other rules."

$OnePassDebug::usage =
"When OnePassOrder detects a recursion among the definitions of a list,
it deposits the offending rules in an internal format in $OnePassDebug
as debugging hints."

Tag::usage =
"Tag[t, expr] tags expr with t (possibly empty).  This tag is 
transparent to the functions MoveDepsLeft, MoveDepsRight, OnePassOrder."

SubroutineDecl::usage =
"SubroutineDecl[name] returns a string with the declaration of the
Fortran subroutine name.  SubroutineDecl[name[args], decl] declares
the subroutine name with arguments args, where decl is a string with
the declaration of args."

VarDecl::usage =
"VarDecl[v, t] returns a string with the declaration of v as variables
of type t in Fortran.  VarDecl[v, {t, c}] additionally makes v members
of common block c.  VarDecl[v, t, f] wraps the declaration in
#ifndef f ... #endif."

CommonDecl::usage =
"The functionality of CommonDecl[v, t, c] has been integrated into
VarDecl[v, {t, c}]."

DoDecl::usage =
"DoDecl[v, m] returns two strings with the do/enddo declarations of a
Fortran loop over v from 1 to m.  DoDecl[v, {a, b}] returns the same for
a loop from a to b.  DoDecl[v] invokes Dim[v] to determine the upper
bound on v."

CallDecl::usage =
"CallDecl[names] returns a string with the invocations of the
subroutines names in Fortran, taking into account possible loops
indicated by DoLoop."

$SymbolPrefix::usage =
"$SymbolPrefix is a string prepended to all externally visible symbols
in the generated Fortran code to avoid symbol collisions."


(* symbols used in the Fortran code *)

Ctree::usage =
"Ctree[Fi] is the ith form factor (the coefficient of Fi) of the
tree-level amplitude."

Cloop::usage =
"Cloop[Fi] is the ith form factor (the coefficient of Fi) of the
one-loop amplitude."

SInvariant::usage =
"SInvariant[ki, kj] represents the s-type (invariant-mass type)
invariant formed from the momenta ki and kj, i.e. s_{ij} = (ki + kj)^2."

TInvariant::usage =
"TInvariant[ki, kj] represents the t-type (momentum-transfer type)
invariant formed from the momenta ki and kj, i.e. t_{ij} = (ki - kj)^2."

DCONJG::usage =
"DCONJG[z] takes the complex conjugate of z in Fortran."

DBLE::usage =
"DBLE[z] takes the real part of z in Fortran."

DIMAG::usage =
"DIMAG[z] takes the imaginary part of z in Fortran."

exp::usage =
"exp[x] is the exponential function in Fortran."

cI::usage =
"cI represents the imaginary unit in Fortran."

SxS::usage =
"SxS[s1, s2] is the Fortran function which computes s1.s2, the direct
product of the two Weyl spinors s1 and s2."

SeS::usage =
"SeS[s1, s2] is the Fortran function which computes s1.eps.s2, the SU(2)
product of the two Weyl spinors s1 and s2."

VxS::usage =
"VxS[v, s] is the Fortran function which computes sigma[v].s, the direct
product of the vector v contracted with sigma and the Weyl spinor s."

VeS::usage =
"VeS[v, s] is the Fortran function which computes sigma[v].eps.s, the
SU(2) product of the vector v contracted with sigma and the Weyl spinor
s."

BxS::usage =
"BxS[v, s] is the Fortran function which computes sigmabar[v].s, the
direct product of the vector v contracted with sigma-bar and the Weyl
spinor s."

BeS::usage =
"BeS[v, s] is the Fortran function which computes sigmabar[v].eps.s,
the SU(2) product of the vector v contracted with sigma-bar and the
Weyl spinor s."

SplitChain::usage =
"SplitChain[w] splits WeylChain[w] into primitive operations SxS, SeS,
VxS, VeS, BxS, BeS for computation."

K::usage =
"K[f, i] is the encoded version of momentum i with (integer) prefactor f."


(* system variables *)

$Editor::usage =
"$Editor specifies the editor command line used in debugging FORM code."

$FormCalc::usage =
"$FormCalc contains the version number of FormCalc."

$FormCalcDir::usage =
"$FormCalcDir points to the directory from which FormCalc was loaded."

$FormCalcProgramDir::usage =
"$FormCalcProgramDir points to the directory which contains the FormCalc
program files."

$ReadForm::usage =
"$ReadForm contains the location of the ReadForm executable."

$FormCmd::usage =
"$FormCmd gives the name of the actual FORM executable.  It may contain
a path."

$DriversDir::usage =
"$DriversDir is the path where the driver programs for the generated
Fortran code are located."

$BlockSize::usage =
"$BlockSize is the maximum LeafCount a single Fortran statement written
out by WriteExpr may have.  Any expression with LeafCount > $BlockSize
will be chopped up before being written to the Fortran file."

$FileSize::usage =
"$FileSize gives the maximum LeafCount the expressions in a single
Fortran file may have.  If the expressions grow larger than $FileSize,
the file is split into several pieces."


Begin["`Private`"]

$FormCalc = 6.1

$FormCalcDir = DirectoryName[ File /.
  FileInformation[System`Private`FindFile[$Input]] ]

$FormCalcProgramDir = ToFileName[{$FormCalcDir, "FormCalc"}]


$ReadForm = ToFileName[{$FormCalcDir, $SystemID}, "ReadForm"];

Check[
  Install[$ReadForm],
  ReadForm::notcompiled = "The ReadForm executable `` could not be \
installed.  Did you run the compile script first?";
  Message[ReadForm::notcompiled, $ReadForm];
  Abort[] ]


$NumberMarks = False

Off[General::spell1, General::spell]


If[ $VersionNumber < 6,
  Needs["Algebra`Horner`"];

  (* actually load the Horner package so that the Off works: *)
  Algebra`Horner[1];
  Off[Algebra`Horner::fail];

  System`HornerForm = Algebra`Horner`Horner;

  Unprotect[StringMatchQ];
  StringMatchQ[s_, l_List] := !VectorQ[l, !StringMatchQ[s, #]&];
  Protect[StringMatchQ]
]


SetOptions[ToString, CharacterEncoding -> "ASCII"]


(* generic functions *)

ParseOpt[func_, opt___] :=
Block[ {names = First/@ Options[func]},
  Message[func::optx, #, func]&/@
    Complement[First/@ {opt}, names];
  names /. {opt} /. Options[func]
]


ToSymbol[x__] := ToExpression[ StringJoin[ToString/@ Flatten[{x}]] ]


Attributes[ToArray] = {Listable}

ToArray[sym_Symbol] :=
Block[ {t = ToString[sym], p},
  p = Position[DigitQ/@ Characters[t], False][[-1, 1]];
  ToExpression[StringTake[t, p]] @ ToExpression[StringDrop[t, p]]
]

ToArray[other_] = other

ToArray[expr_, vars__] :=
Block[ {v = ToString/@ Flatten[{vars}], rules, t, p, h},
  rules = (
    t = ToString[#];
    p = Position[DigitQ/@ Characters[t], False][[-1, 1]];
    h = StringTake[t, p];
    If[ MemberQ[v, h],
      # -> ToExpression[h] @ ToExpression[StringDrop[t, p]],
      {} ]
  )&/@ Symbols[expr];
  expr /. Dispatch[Flatten[rules]]
]


Renumber[expr_, vars__] :=
Block[ {v = Flatten[{vars}], old, new},
  old = Union[Cases[expr, #[__], Infinity]]&/@ v;
  new = MapThread[Array, {v, Length/@ old}];
  expr /. Dispatch[Thread[Flatten[old] -> Flatten[new]]]
]


ExpandSums[x_, ___] := x /; FreeQ[x, SumOver]

ExpandSums[a_. IndexDelta[i_, j_] SumOver[i_, _], h___] :=
  ExpandSums[a /. i -> j, h]

ExpandSums[a_. s__SumOver, h_:Sum] := h[a, Sequence@@ List@@@ {s}]

ExpandSums[other_, h___] := ExpandSums[#, h]&/@ other


Kind[Tag[___, x_]] := Kind[x]

Kind[h_ -> _] := kind[h]

_Kind = Sequence[]

kind[h_[___], ___] = h

kind[h_, ___] = h


KindPattern[Tag[___, x_]] := KindPattern[x]

KindPattern[h_ -> _] := kindPattern[h]

_KindPattern = Sequence[]

kindPattern[h_[___], ___] = _h

kindPattern[h_, ___] = h


Alt[l_List] := Alt@@ Flatten[l]

Alt[s_] = s

Alt[s__] := Alternatives[s]


MaxDims[args__] :=
  listmax/@ Split[Union[Flatten[{args}]], Head[#1] === Head[#2]&]

listmax[{s__Symbol}] := s

listmax[{s__String}] := s

listmax[l_] := l[[1, 0]]@@ MapThread[Max, List@@@ l]


Attributes[Keep] = {HoldFirst}

Keep[lhs_ = rhs_] := lhs = Keep[rhs, Block[{lhs}, ToString[lhs]]]

Keep[lhs_ = rhs_, other__] := lhs = Keep[rhs, other]

Keep[cmd_, name_String, prefix_String:$KeepDir] :=
Block[ {file = ChkExist[prefix, name <> ".m"]},
  If[ FileType[file] === File,
    Print["loading ", file];
    (RegisterAbbr[#1]; RegisterSubexpr[#2]; #3)&@@ Get[file],
  (* else *)
    (Put[{Abbr[], Subexpr[], #}, file]; #)& @ cmd ]
]

$KeepDir = "keep/"


Attributes[ToForm] = {Listable}

ToForm[x_String] = x

ToForm[x_] := ToString[x, InputForm]


ToSeq[li_List] := StringTake[ToString[li, InputForm], {2, -2}]

ToSeq[x_] := ToForm[x]


ToBool[True] = "1"

ToBool[___] = "0"


ToFortran[x_String] = x

ToFortran[x_List] := StringTake[ToString[x, FortranForm], {6, -2}]

ToFortran[x_] := ToString[x, FortranForm]


ToCat[n_, {}] := Table[{}, {n}]

ToCat[_, li_] := Flatten/@ Transpose[li]


Symbols[expr_] := Union[Cases[expr, _Symbol, {-1}]]


FromPlus[h_, p_Plus] := h@@ p

FromPlus[_, other_] = other


numadd[term_] := numadd[
  Numerator[term],
  Denominator[term] //Simplify ]

numadd[n_, x_?NumberQ d_] := numer[d] += n/x

numadd[n_, d_] := numer[d] += n

DenCollect[p_Plus, wrap_:Identity] :=
Block[ {numer},
  _numer = 0;
  numadd/@ p;
  _numer =.;
  Plus@@ (wrap[#2]/#1[[1, 1]] &)@@@ DownValues[numer]
]

DenCollect[x_, wrap_:Identity] := wrap[x] /; FreeQ[x, Plus]

DenCollect[x_, wrap___] := DenCollect[#, wrap]&/@ x


Pool[expr_, wrap_:Identity] := expr /.
  p_Plus :> ploos[wrap]@@ Cases[p, _Times] +
    DeleteCases[p, _Times] /; LeafCount[p] > 10

ploos[wrap_][a_, r__] :=
Block[ {pos = 0, lcmin, lcmax, lc, ov, ovmax, ovpos, ploos},
  lcmin = Floor[LeafCount[a]/3];
  lcmax = -Infinity;
  Scan[ (
    ++pos;
    lc = LeafCount[ov = Intersection[a, #]];
    If[ lc > lcmax,
      lcmax = lc; ovpos = pos; ovmax = ov;
      If[ lc > lcmin, Return[] ] ] )&, {r} ];
  If[ lcmax < 5,
    a + ploos[wrap][r],
  (* else *)
    ovmax wrap[a/ovmax + {r}[[ovpos]]/ovmax] +
      Drop[ploos[wrap][r], {ovpos}] ]
]

ploos[_][other___] := Plus[other]


Attributes[usq] = {Orderless}

usq/: usq[a__][ik__] + usq[b__][ik__] := usq[a, b][ik]

usq/: usq[a__][h_, i_, i_] + usq[i_][h_, k_, k_] := usq[a, k][h, i, i]


usum[i___, _, _, _, _] := {i} /; uexpr[[i, 0]] === Plus

_usum = Sequence[]


ApplyUnitarity[expr_, U_, dim_Integer, simp_:FullSimplify] :=
Block[ {ux, uexpr, pos},
  Evaluate[ux@@ Range[dim]] = KroneckerDelta[##2] &;
  ux[a__] := With[ {cpl = usq@@ Complement[Range[dim], {a}]},
    ux[a] = KroneckerDelta[##2] - cpl[##] &
  ] /; Length[{a}] > dim/2;
  ux[other__] := usq[other];

  uexpr = expr //. {
    (u:U[i_, j_])^n_. (uc:Conjugate[U[k_, j_]])^nc_. :>
      (usq[j][1, i, k]^# u^(n - #) uc^(nc - #) &) @ Min[n, nc],
    (u:U[j_, k_])^n_. (uc:Conjugate[U[j_, i_]])^nc_. :>
      (usq[j][2, i, k]^# u^(n - #) uc^(nc - #) &) @ Min[n, nc]
  } /. usq -> ux;

  pos = usum@@@ Position[uexpr, usq];
  pos = Select[pos, FreeQ[pos, Append[#, __]]&];

  MapAt[simp, uexpr, pos] /. usq -> ux /. {
    usq[a__][1, i_, k_] -> Plus@@ (U[i, #] Conjugate[U[k, #]] &)/@ {a},
    usq[a__][2, i_, k_] -> Plus@@ (U[#, k] Conjugate[U[#, i]] &)/@ {a} }
]


foodef[Map] := Hold[foo[y_] := foo/@ y]

foodef[f_:Identity] := Hold[foo[y_] := f[y]]

foodef[n_, f_, r___] := {
  Hold[foo[y_] := f[y] /; LeafCount[y] < n],
  foodef[r]
}

Attributes[blkdef] = {HoldAll}

blkdef[defs___] := Block[{foo}, defs; foo[#]]&

OnSize[args__] := Level[Flatten[{foodef[args]}], {2}, blkdef]


DiagramType[a_FeynAmp] := Count[a[[3]] /. _FeynAmpDenominator -> 1,
  _PropagatorDenominator, Infinity]


FermionicQ[a_FeynAmp] := !FreeQ[a[[3]], FermionChain | MatrixTrace]


Attributes[IndexDiff] = Attributes[IndexDelta] = {Orderless}

IndexDelta[i_, i_] = 1

IndexDelta[i_Integer, _Integer] = 0

IndexDiff[i_, i_] = 0

IndexDiff[i_Integer, _Integer] = 1


IndexIf[] = 0

IndexIf[a_] = a

IndexIf[True, a_, ___] = a

IndexIf[False, _, b___] := IndexIf[b]

IndexIf[a__] := IndexIf[a, 0] /; EvenQ[Length[{a}]]

IndexIf[cond_, 0, a_] := IndexIf[!cond, a, 0]

IndexIf[cond1_, IndexIf[cond2_, a_, 0], 0] := IndexIf[cond1 && cond2, a, 0]

IndexIf[a___, i_IndexIf] := IndexIf[a, Sequence@@ i]


MapIf[foo_, i_] := MapAt[foo, i,
  List/@ Append[Range[2, # - 1, 2], #]]& @ Length[i]


Off[Optional::opdef]

ToIndexIf[expr_, s_String] := ToIndexIf[expr, Alt[Names[s]]]

ToIndexIf[expr_, patt_:_] :=
  Fold[ singleIf, expr, Union @ Cases[expr,
    (IndexDelta | IndexDiff)[i:patt..] -> {i}, Infinity] ]

singleIf[expr_, {i__}] := expr /.
  { IndexDelta[i] -> suck[1, 1],
    IndexDiff[i] -> suck[2, 1] } /.
  a_. suck[1, x_] + a_. suck[2, y_] -> a IndexIf[Equal[i], x, y] /.
  suck[h_, x_] :> IndexIf[{Equal, Unequal}[[h]][i], x]

suck/: r_ suck[h_, x_] := suck[h, r x] /; FreeQ[r, Mat (*| SumOver*)]

suck/: suck[h_, x_] + suck[h_, y_] := suck[h, x + y]

(* suck/: suck[h_, x_] + y_ := suck[3 - h, y] /; x + y == 0 *)


(* preparations for FORM *)

ExtWF = {e, ec, z, zc, eT, eTc}

ConjWF = ({#1 -> #2, #2 -> #1}&)@@@ Partition[ExtWF, 2] //Flatten

KinVecs = {eta, e, ec, z, zc, k}

KinTens = {eT, eTc}

KinObjs = Join[KinVecs, KinTens]
	(* note: KinObjs determines the order in which vectors
	   and tensors appear in FORM and hence ultimately in
	   functions like Pair and Eps *)


Attributes[KinFunc] = {HoldAll}

KinFunc[args__] := Function[Evaluate[KinObjs], args]


FromFormRules = Outer[ ToSymbol["FormCalc`", #2, #1] -> #2[#1] &,
  Range[8], KinObjs ]

FormKins = Apply[#1&, FromFormRules, {2}]

FromFormRules = Flatten[FromFormRules]


MomThread[i_Index, foo_] := foo[i]

MomThread[p_Symbol, foo_] := foo[p]

MomThread[p_, foo_] := Replace[MomReduce[p], k_Symbol :> foo[k], {-1}]


MomReduce[p_Plus] := Fewest[p, p + MomSum, p - MomSum]

MomReduce[p_] = p


Fewest[a_, b_, r___] := Fewest[a, r] /; Length[a] <= Length[b]

Fewest[_, b__] := Fewest[b]

Fewest[a_] = a


fvec[p_] := (vecs = {vecs, Symbols[p]}; MomReduce[p])

fvec[p_, mu_] := (vecs = {vecs, Symbols[p]}; MomThread[p, #[mu]&])


iname[type_, n_] := iname[type, n] =
  ToSymbol[StringTake[ToString[type], 3], n]


Attributes[idelta] = {Orderless}

idelta[c1:Index[Colour, _], c2_] := SUNT[c1, c2]

idelta[g1:Index[Gluon, _], g2_] := 2 SUNT[g1, g2, 0, 0]

idelta[x__] := IndexDelta[x]


ieps[c__] := Block[{eps = SUNEps[c]}, eps /; !FreeQ[eps, Index[Colour, _]]]

ieps[x__] := IndexEps[x]


psum/: psum[n_][x__]^k_ := psum[n k][x]

isum[expr_, {i_, f_:1, t_}] := (t - f + 1) expr /; FreeQ[expr, i]

isum[expr_, {i_, r__}] :=
Block[ {dummy = Unique[ToString[i]]},
  indices = {indices, dummy};
  (expr /. i -> dummy) SumOver[dummy, r, Renumber]
]

isum[expr_, i_] :=
Block[ {dummy = Unique[ToString[i]]},
  indices = {indices, dummy};
  ranges = {dummy -> dummy, ranges};
  (expr /. i -> dummy)
]


sumover[i_, 1, r_] := sumover[i, r]

sumover[i:Index[Colour | Gluon, _], r__] := SUNSum[i, r]

sumover[other__] := SumOver[other]

SUNSum[_, _, External] = 1


KinFunc[
  pol[_, k, mu:Index[EpsilonScalar, _]] = z[mu];
  polc[_, k, mu:Index[EpsilonScalar, _]] = zc[mu];
  pol[_, k, mu_] = e[mu];
  polc[_, k, mu_] = ec[mu];
  pol[_, k, mu__] = eT[mu];
  polc[_, k, mu__] = eTc[mu] ]@@@ FormKins

Conjugate[pol] ^= polc


Attributes[scalar] = {Orderless}

scalar[0, _] = 0

scalar[a_Symbol, b_Symbol] := a . b

scalar[a_, p:_[__]] := MomThread[p, scalar[a, #]&]


prop[0, m_Symbol, d___] := -m^(-2 d)

prop[p_, m__] := prop[-p, m] /; !FreeQ[p, -q1]

prop[p_, m_, d___] := Den[p, m^2, d]


Attributes[loop] = {Orderless}

loop[a___, d_[p_, m1_], _[p_, m2_], b___] :=
  (loop[a, d[p, m1], b] - loop[a, d[p, m2], b]) inv[Factor[m1 - m2]] /;
  m1 =!= m2

loop[d__] := I Pi^2 intM[d]


inv[t_Times] := inv/@ t

inv[p_Plus] := inv[p, -p]

inv[other_] := 1/other

inv/: inv[p_, _] p_ = 1

inv/: inv[_, p_] p_ = -1


noncomm[p_Plus] := noncomm/@ p

noncomm[g_] := g ga[] /; FreeQ[g, ga | Spinor]

noncomm[g_] = g

noncomm[g__] := NonCommutativeMultiply[g]


Neglect[m_] = m

FormPatt[_[_, m]] = {}

FormPatt[_?NumberQ, _] = {}

FormPatt[_[_[_[lhs_]], rhs_]] := FormPatt[lhs, rhs]

FormPatt[lhs_Alternatives, rhs_] := FormPatt[#, rhs]&/@ List@@ lhs

FormPatt[lhs_, rhs_] :=
Block[ {c = 96, newlhs, newrhs = rhs},
  newlhs = lhs /.
    {Blank -> FormArg, BlankSequence | BlankNullSequence -> FormArgs} /.
    Pattern -> ((newrhs = newrhs /. #1 -> #2; #2)&);
  (newlhs /. patt[x_] :> x <> "?") -> (newrhs /. patt[x_] -> x)
]

FormArg[h_:Identity] := h[patt["ARG" <> FromCharacterCode[++c]]]

FormArgs[h_:Identity] := h["?" <> FromCharacterCode[++c]]


OrdSq[r:_[_, rhs_]] := {{}, r} /; VectorQ[lhs, FreeQ[rhs, #]&]

OrdSq[r_] := {r, {}}

SortSq[dv_] :=
Block[ {lhs = #[[1, 1, 1]]&/@ dv},
  Flatten[Transpose[OrdSq/@ dv]]
]


Attributes[Inv] = {Listable}

Inv[i_, j_] :=
Block[ {s = signs[[i]] signs[[j]], ki = Moms[[i]], kj = Moms[[j]], inv},
  inv = If[ legs === 3, dot[#, #]& @ Moms[[3]], Invariant[s, i, j] ];
  dot[ki, kj] = s/2 inv - s/2 dot[ki, ki] - s/2 dot[kj, kj];
  inv
]


Invariant[1, 1, 2] = S

Invariant[-1, 1, 3] = T

Invariant[-1, 2, 3] = U

Invariant[s_, i_, j_] := Invariant[s, i, j] =
  ToSymbol["FormCalc`", FromCharacterCode[(167 - s)/2], i, j]


OtherProd[{k1___, k_}, {s1___, s_}, zero_] :=
  MapThread[
    (dot[#1, k] = -s (Plus@@ ({s1} dot[#1, {k1}]) - #2 zero))&,
    {{k1}, {s1}} ]


(* global variables set here:
   CurrentProcess, CurrentOptions, CurrentScale,
   FormProcs, FormSymbols,
   Moms, MomSum, InvSum, PairRules, LastAmps *)

CurrentProcess = CurrentOptions = {}

Options[DeclareProcess] = {
  OnShell -> True,
  Invariants -> True,
  Transverse -> True,
  Normalized -> True,
  InvSimplify -> True,
  Antisymmetrize -> True,
  MomElim -> Automatic,
  AbbrScale -> 1 }

Attributes[DeclareProcess] = {Listable}

DeclareProcess::incomp =
"Calculation of incompatible process(es) attempted.  If you want to \
calculate a new process, run ClearProcess[] first."

DeclareProcess::syntax =
"Wrong syntax: DeclareProcess expects FeynAmpList objects as arguments."

DeclareProcess[fal:FeynAmpList[__][___].., opt___Rule] :=
Block[ {proc},
  proc = Process /. List@@@ Head/@ {fal};
  (DeclP@@ Flatten[{Apply[#3 &, proc, {3}], CurrentProcess}])[
    ParseOpt[DeclareProcess, opt], CurrentOptions ];
  LastAmps = {fal};
  Level[ {fal}, {2}, FormAmp @ proc[[1]] ]
]

DeclareProcess[l__List, opt___Rule] := DeclareProcess@@ Flatten[{l, opt}]

DeclareProcess[___] := (Message[DeclareProcess::syntax]; Abort[])


DeclP[(proc_)..] := DeclO[proc]

_DeclP := (Message[DeclareProcess::incomp]; Abort[])


_DeclO[opt_, opt_] = Null

DeclO[proc_][opt:{onshell_, inv_, transv_, norm_,
                  invsimp_, antisymm_, momel_, scale_}, ___] :=
Block[ {momelim = momel, signs, masses, legs, kins, dot, n,
neglect, square, invproc = {}, invs = {},
kikj = {}, eiki = {}, eiei = {}},
  CurrentProcess = proc;
  CurrentOptions = opt;
  CurrentScale = scale;
  MomSum = InvSum = 0;
  FormProcs = {};
  FormSymbols = {scale,
    neglect = FormPatt/@ DownValues[Neglect],
    square = FormPatt/@ SortSq[DownValues[Square]]};

  signs = Flatten[{1&/@ proc[[1]], -1&/@ proc[[2]]}];
  masses = Level[proc, {2}]^2;
  legs = Length[masses];
  kins = Take[FormKins, legs];

  Moms = KinFunc[k]@@@ kins;
  FormVectors =
    Level[{q1, Transpose[KinFunc[Evaluate[KinVecs]]@@@ kins]}, {-1}];
  FormTensors =
    Level[KinFunc[Evaluate[KinTens]]@@@ kins, {-1}];

  Attributes[dot] = {Orderless, Listable};

  kikj = dot[Moms, Moms];
  If[ onshell, kikj = MapThread[Set, {kikj, masses}] ];

  If[ transv, eiki = KinFunc[{e.k -> 0, ec.k -> 0}]@@@ kins ];

  If[ norm, eiei = KinFunc[e.ec -> -1]@@@ kins ];

  Switch[ Length[kins],
    0 | 1,
      Null,
    2,
      (*dot[ Moms[[2]], Moms[[2]] ] =.;*)
      eiki = eiki /. Moms[[2]] -> Moms[[1]];
      FormProcs = {FormProcs, "\
id Spinor(k2, ?m) = Spinor(e2, ?m);\n\
multiply replace_(k2, k1);
id Spinor(e2, ?m) = Spinor(k2, ?m);\n\n"};
      FormSymbols = {FormSymbols, Spinor[]};
      momelim = False,
    _,
      If[ inv,
        invs = Flatten[Array[Inv[Range[# - 1], #]&, legs - 2, 2]];
        Scan[(RealQ[#] = True)&, invs];
        InvSum = kikj[[-1]] + Plus@@ ((legs - 3) Drop[kikj, -1]);

	(* The number of invariants is ninv = (legs - 1)(legs - 2)/2 in
	   total and ndot = (legs - 2) in a dot product (of the terms in
	   pi.plegs = Sum[pi.pj, {j, legs - 1}], pi.pi is a mass^2).
	   Thus, the number of invariants in pi.plegs can be reduced by
	   using the Mandelstam relation only if ndot > ninv/2, that is
	   legs < 5.  The case legs = 3 is handled specially already by
	   Inv, so only legs = 4 remains. *)
        OtherProd[ Moms, signs,
          If[legs === 4, Distribute[(Plus@@ invs - InvSum)/2], 0] ];

        If[ invsimp && onshell && legs =!= 3,
          n = Length[invs] + 1;
          invproc = ToForm[MapIndexed[
            {"id `foo'(ARG?) = `foo'(ARG, replace_(",
              #1, ", ",  InvSum - Plus@@ Drop[invs, #2],
              ")*ARG);\n#call Fewest(`foo')\n"}&, invs ]]
        ];

	  (* not used anywhere in FormCalc, but useful for debugging: *)
        InvSum = Plus@@ invs - InvSum
      ];

      n = signs Moms;
      MomSum = Plus@@ n;
      FormProcs = {FormProcs,
        "#define MomSum \"", MomSum, "\"\n",
        Array[{"#define ", Moms[[#]], " \"",
          -signs[[#]] Plus@@ Drop[n, {#}], "\"\n"}&, legs]};
  ];

  kikj = ReleaseHold[DownValues[dot] /. dot -> Dot];
  FormSymbols = Flatten[{FormSymbols, Last/@ kikj}];

  FormProcs = "\
#define Legs \"" <> ToString[legs] <> "\"\n\
#define MomElim \"" <> ToString[momelim] <> "\"\n\
#define OnShell \"" <> ToBool[onshell] <> "\"\n\
#define Antisymmetrize \"" <> ToBool[antisymm] <> "\"\n\
#define Scale \"" <> ToForm[scale] <> "\"\n\n" <>
    ToForm[FormProcs] <> "\
#procedure Neglect\n" <> FormId[neglect] <> "#endprocedure\n\n\
#procedure Square\n\
id powM(ARG?, ARG1?)^ARG2? = powM(ARG, ARG1*ARG2);\n\
repeat id powM(ARG?, ARG1?) * powM(ARG?, ARG2?) = powM(ARG, ARG1 + ARG2);\n\
id powM(ARG?, ARG1?int_) = ARG^ARG1;\n" <>
    FormSq[square] <> "#endprocedure\n\n\
#procedure eiei\n" <> FormId[eiei] <> "#endprocedure\n\n\
#procedure eiki\n" <> FormId[eiki] <> "#endprocedure\n\n\
#procedure kikj\n" <> FormId[kikj] <> "#endprocedure\n\n\
#procedure InvSimplify(foo)\n" <> invproc <> "#endprocedure\n\n";

	(* not used anywhere in FormCalc, but useful for debugging: *)
  PairRules = Flatten[{eiei, eiki, kikj}] /. Dot -> Pair /. FromFormRules;
]


LevelSelect[_][id_, _, amp_] :=
Block[ {name = AmpName[Select[id, FreeQ[#, Number]&]]},
  {name -> TrivialSums[GenNames[amp]], name}
]

LevelSelect[Automatic][r__, gm_ -> ins_] := LevelSelect[
  Which[
    !FreeQ[ins, Particles], Particles,
    !FreeQ[ins, Classes], Classes,
    True, Generic ]
][r, gm -> ins]

LevelSelect[Generic][id_, _, gen_, ___] :=
Block[ {name = AmpName[Select[id, FreeQ[#, Classes | Particles | Number]&]]},
  {name -> GenNames[gen], name}
]

LevelSelect[lev_][id_, _, gen_, gm_ -> ins_] :=
Block[ {old, new, amp, name = AmpName[id], pc},
  _pc = 0;
  old = TrivialSums/@ Cases[{ins}, Insertions[lev][r__] :> r, Infinity];
  new = Thread[Flatten[ReduceIns[gm, Transpose[old]]], Rule];
  amp = gen /. new[[1]] /.
    {p_PropagatorDenominator :> (p /. small[m_] -> m), _small -> 0};
  { name -> amp,
    If[ Length[ new[[2]] ] === 0,
      Length[old] name,
    (* else *)
      "i" <> name -> name Plus@@ (Level[#, {2}, "replace_"]&)/@
        Transpose[ Thread/@ new[[2]] ] ] }
]


ReduceIns[{g_, rg___},
  {{ins_ /; (* LeafCount[ins] < 10 && *)
		(* Commuting functions in FORM are only commuting in
		   nonnegative powers.  The following FreeQ makes sure
		   that all expressions containing such are always
		   inserted *after* the kinematical simplification so
		   that FORM won't screw up the spinor chains. *)
            FreeQ[ins, _[__]^_?Negative], r___} /; ins === r, rins___}] :=
Block[ {smallg},
  smallg = If[ Neglect[ins] === 0, small[ins], ins ];
  { (g -> smallg) -> Sequence[], ReduceIns[{rg}, {rins}] }
]

ReduceIns[{g_, rg___}, {ins_List, rins___}] :=
Block[ {instype = InsType[g], newg, smallg},
  If[ FreeQ[ins, DenyHide],
    newg = ToSymbol["Form`p", instype, ++pc[instype]],
  (* else *)
    patt = {patt, newg = ToSymbol["FormCalc`c", instype, ++pc[instype]]} ];
  smallg = If[ Union[Neglect/@ ins] === {0}, small[newg], newg ];
  { (g -> smallg) -> (newg -> inssym/@ ins),
    ReduceIns[{rg}, Replace[{rins}, {ins -> smallg, -ins -> -smallg}, 1]] }
]

ReduceIns[{g_, rg___}, {newg_, rins___}] :=
  { (g -> newg) -> Sequence[], ReduceIns[{rg}, {rins}] }

ReduceIns[{}, {}] = {}


InsType[_Mass] = "M"

InsType[RelativeCF] = "R"

InsType[_GaugeXi] = "X"

InsType[_] = "G"


AmpName[g_] :=
Block[ {t, c},
  t = StringJoin@@ ({StringTake[ToString[#1], 1], ToString[#2]}&)@@@ g;
  If[ (c = ++uniq[t]) === 1, t, t <> "v" <> ToString[c] ]
]


TrivialSums[ins_ -> _] := TrivialSums[ins]

TrivialSums[ins_] := ins /; FreeQ[ins, SumOver]

TrivialSums[ins_] :=
Block[ {test = ins /. _SumOver -> 1},
  ins /. SumOver -> CarryOut
]

CarryOut[i_, n_, r___] := (
  ranges = {i -> i == n, ranges};
  sumover[i, n, r]
) /; !FreeQ[test, i]

CarryOut[_, v_, External] := Sqrt[v]

CarryOut[_, v_, ___] = v


Attributes[OffShell] = {Listable}

OffShell[fal:FeynAmpList[__][___].., extmass__Rule] :=
Block[ {subst},
  (subst[#1][_] = #2)&@@@ {extmass};
  subst[_][m_] = m;
  ##&@@ ({fal} /. (Process -> p_) :>
    Block[{c = 0}, Process -> Map[MapAt[subst[++c], #, 3]&, p, {2}]])
]


Combine::incomp = "Warning: combining incompatible processes."

Combine[fal:FeynAmpList[__][___]..] :=
Block[ {amps = {fal}},
  If[ !SameQ@@ Cases[Head/@ amps, _[Process, p_] -> p, {2}],
    Message[Combine::incomp] ];
  Level[ amps, {1}, amps[[1, 0]] ]
]

Combine[amp:Amp[_][___]..] :=
Block[ {amps = {amp}, comp},
  If[ !Level[Head/@ amps, {1}, SameQ],
    Message[Combine::incomp] ];
  _comp = 0;
  Map[Component, amps, {2}];
  _comp =.;
  amps[[1, 0]]@@ Cases[DownValues[comp], _[_[_[x___]], r_] -> r x]
]

Component[r_ s__SumOver] := comp[s] += r

Component[other_] := comp[] += other


m_MultiplyDiagrams[fal:FeynAmpList[__][___]] := m/@ fal

m_MultiplyDiagrams[fal:FormAmp[__][___]] := m/@ fal

MultiplyDiagrams[f_][FeynAmp[id_, q_, gen_, gm_ -> ins_]] :=
  FeynAmp[id, q, gen, gm -> MultiplyIns[f, gen, gm]/@ ins]

MultiplyDiagrams[f_][FeynAmp[id_, q_, amp_]] :=
  FeynAmp[id, q, f[amp] amp]

i_MultiplyIns[ins_ -> more_] := i[ins] -> i/@ more

MultiplyIns[f_, gen_, gm_][{r__, fac_}] := {r, f[gen, gm -> {r, fac}] fac}


TagDiagrams[diags_, tag_:Diagram] :=
Block[ {diag = 0},
  MultiplyDiagrams[tag[++diag]&][diags]
]


IndexDim[Lorentz, i_] := i -> i

IndexDim[Lorentz4, i_] := i -> i == 4

IndexDim[EpsilonScalar, i_] := i -> i == Dminus4

_IndexDim = {}


(* the main function CalcFeynAmp *)

SUNObjs = SUNSum | SUNT | SUNTSum | SUNF | SUNEps

DenyNoExp = Level[{ga, Spinor, Den, A0i, B0i, C0i, D0i, E0i, F0i,
  SumOver, SUNObjs, "replace_"}, {-1}]

DenyHide = Level[{SumOver, IndexDelta, IndexEps, SUNObjs},
  {-1}, Alternatives]

FinalFormRules = {
  (f:_[__])^(n_?Negative) -> powM[f, n],
  x_^n_ :> powM[x, n] /; !IntegerQ[n],
  Complex[a_, b_] -> a + "i_" b }


Attributes[CalcFeynAmp] = {Listable}

Options[CalcFeynAmp] = {
  CalcLevel -> Automatic,
  Dimension -> D,
  FermionChains -> Weyl,
  FermionOrder -> Automatic,
  InsertionPolicy -> Default,
  CutTools -> False,
  NoExpand -> {},
  NoBracket -> {},
  EditCode -> False,
  RetainFile -> False }

CalcFeynAmp::ddim = "Warning: `` \
implies the use of Fierz identities which are in general not \
applicable in D dimensions.  You know what you are doing."

CalcFeynAmp[FormAmp[proc_][amp___], opt___Rule] :=
Block[ {lev, dim, fchain, forder, inspol, cuttools, noexp, nobrk,
edit, retain, uniq, vecs, legs = 0, ic, inssym, mmains, 
indices = {}, ranges = {}, indsym, vars, patt = {}, 
hh, amps, res, traces = 0},

  {lev, dim, fchain, forder, inspol, cuttools, noexp, nobrk,
    edit, retain} = ParseOpt[CalcFeynAmp, opt];

  If[ dim === 0,
    If[ Length[#] > 0, Message[CalcFeynAmp::ddim, #] ]&[{
      If[ forder =!= None, FermionOrder -> forder, {} ],
      If[ fchain === Weyl, FermionChains -> fchain, {} ]
    } //Flatten] ];

(*
  NoExp[p_] := NoExp[p] = Unique["noexp"];
  NoExpandRule = (p_Plus :> NoExp[p] /; !FreeQ[p, #1] && FreeQ[p, #2])&[
    Alt[noexp],
    Alt[{FormVectors, FormTensors, DenyNoExp}] ];
*)

  NoExpandRule = (p_Plus :> (p /. Plus -> addM) /;
    !FreeQ[p, #1] && FreeQ[p, #2])&[
      Alt[noexp],
      Alt[{FormVectors, FormTensors, DenyNoExp}] ];

  If[ NumberQ[inspol],
    _ic = 0;
    inssym[ins_] := ins /; LeafCount[ins] < inspol || !FreeQ[ins, DenyHide];
    inssym[ins_] := inssym[ins] =
      ToSymbol["FormCalc`i", instype, ++ic[instype]],
  (* else *)
    inssym = Identity ];

  _uniq = 0;
  vecs = {};
  amps = LevelSelect[lev]@@@ ({amp} /.
      {g:G[_][_][__][_] -> g,
       IndexDelta -> idelta,
       IndexEps -> ieps,
       IndexSum -> psum[1]} //.
      psum[n_] -> (Product[Fold[isum, #1, {##2}], {n}]&)) /.
    Cases[Thread[(#2&@@@ Level[proc, {2}]) -> Moms],
      _[_FourMomentum, _]] /.
    FourMomentum[Internal, 1] -> q1 /.
    PolarizationVector | PolarizationTensor -> pol /.
    { _MetricTensor^2 -> dim,
      MetricTensor -> "d_",
      LeviCivita -> Eps,
      ScalarProduct -> scalar,
      PropagatorDenominator -> prop,
      FeynAmpDenominator -> loop,
      (DiracSpinor | MajoranaSpinor)[s_. p_Symbol, m_] :>
        Spinor[p, Neglect[m], s],
      ChiralityProjector[c_] :> ga[(13 - c)/2],
      (DiracSlash | DiracMatrix)[p_] :> MomThread[p, ga],
      NonCommutative -> noncomm,
      FourVector -> fvec } /.
    inv[p_, _] -> 1/p;

  FormVectors = Join[FormVectors,
    Complement[Flatten[vecs], FormVectors]];

  amps = DeleteCases[amps, {___, _ -> 0, __}];
  If[ Length[amps] === 0, Return[Amp[CurrentProcess][0]] ];

  mmains = Cases[DownValues[inssym], _[_[_[ins_]], s_Symbol] :> s == ins];

  amps = ToCat[2, amps] /. NoExpandRule /. FinalFormRules;

  indsym[type_, n_] := (
    indices = {indices, #};
    ranges = {IndexDim[type, #], ranges};
    indsym[type, n] = #
  )& @ iname[type, n];

  amps = amps /. Index -> indsym;
  mmains = mmains /. Index -> indsym;
  proc /. Index -> indsym;
  ranges = Append[Union[Flatten[ranges]] /. Index -> indsym, i_ -> i == 0];

  indices = Union[Flatten[ {indices, sM,
    Cases[DownValues[indsym], _[_, s_Symbol] -> s],
	(* possibly some colour or gluon indices have been fixed by
	   the user in FeynArts; these would by default be declared
	   as symbols and not be seen by FORM: *)
    Cases[amps, SUNObjs[o__] :> Cases[{o}, _Symbol], Infinity]} ]];

  If[ Head[forder] === Colour,
    indices = Join[#, Complement[indices, #]]&[
      iname[Colour, #]&/@ List@@ forder ];
    forder = Head[forder] ];

  vars = DeclareVars[{amps, SUNN, CurrentScale, dirM[]},
    indices, ranges, {D, Dminus4}];

  hh = OpenForm[];
  WriteString[hh,
    If[ dim === 4,
      "s Dminus4;\nd 4;\n",
      "s D, Dminus4;\nd D:Dminus4;\n" ] <>
    vars[[1]] <>
    "\n.global\n\n"];

  FormWrite[ hh, amps[[1]] ];

  WriteString[hh,
    Array[{"trace4,", ToString[#], ";\n"}&, traces] <> "\n\
#define Dim \"" <> ToString[dim] <> "\"\n\
#define InsertionPolicy \"" <> ToForm[inspol] <> "\"\n\
#define CutTools \"" <> ToForm[cuttools] <> "\"\n\
#define HaveFermions \"" <> ToBool[!FreeQ[amps, FermionChain]] <> "\"\n\
#define FermionOrder \"" <> ToSeq[forder] <> "\"\n\
#define FermionChains \"" <> ToForm[fchain] <> "\"\n\
#define HaveSUN \"" <> ToBool[!FreeQ[amps, SUNObjs]] <> "\"\n\
#define SUNN \"" <> ToForm[SUNN] <> "\"\n\
#define Patterns \"" <> ToSeq[Union[Flatten[patt]]] <> "\"\n" <>
    FormProcs <> "\
#procedure Const\n" <>
    FormDecl["ab ", DeleteCases[
      Complement[Flatten @ vars[[{3, 4}]], Flatten[{nobrk}]],
      (x_ /; MemberQ[{"FeynArts`", "FormCalc`", "LoopTools`", "Form`"},
        Context[x]]) | SUNObjs | Conjugate ]] <> "\
#endprocedure\n\n\
#procedure MapAll(foo)\n\
#ifdef `Inserted'
#call `foo'(Result)
#else\n" <>
    ({"#call `foo'(", First[#], ")\n"}&)/@ amps[[1]] <> "\
#endif\n\
#endprocedure\n\n\
#procedure Insertions\n"];
  res = Plus@@ FormWrite[ hh, amps[[2]] ];
  WriteString[hh, ".store\n\n"];
  FormWrite[hh, FormCalc`Result -> res];
  WriteString[hh, "\
#endprocedure\n\n\
#include CalcFeynAmp.frm\n"];
  Close[hh];

  Amp[CurrentProcess]@@ RunForm[mmains][edit, retain][[1]]
]

CalcFeynAmp[amps__, opt___Rule] := CalcFeynAmp[
  DeclareProcess[amps, FilterOptions[DeclareProcess, opt]],
  FilterOptions[CalcFeynAmp, opt] ]


(* FORM interfacing *)

trace[expr__] := (traces = Max[traces, ++fline];
  NonCommutativeMultiply[expr] /. {
    a_. ga[li_] ** ga[1] + a_. ga[li_] ** ga[-1] :> a "g_"[fline, li],
    a_. ga[1] + a_. ga[-1] :> a "g_"[fline],
    ga -> om })

chain[expr__] := (fline = Max[++fline, 10];
  NonCommutativeMultiply[expr] /. ga -> om)

dobj[expr__][di__] :=
Block[ {fline = sM},
  dirM[NonCommutativeMultiply[expr] /. ga -> om, di]
]


om[1] := "g_"[fline]

om[5] := "g5_"[fline]

om[6] := "g6_"[fline]/2

om[7] := "g7_"[fline]/2

om[li___] := "g_"[fline, li]


Attributes[FormWrite] = {Listable}

FormWrite[hh_, lhs_ -> rhs_] :=
Block[ {fline = 0},
  Write[hh, "G ", lhs == (rhs /.
    MatrixTrace -> trace /.
    FermionChain -> chain /.
    DiracObject -> dobj /.
    NonCommutativeMultiply[a_] -> a), ";"];
  WriteString[hh, "\n"];
  lhs
]

FormWrite[_, lhs_] = lhs


FormDecl[_, _[]] = {}

FormDecl[type_, _[f_, v___]] :=
Block[ {l, ll, t},
  ll = StringLength[t = type <> ToForm[f]];
  { Fold[
      ( l = StringLength[t = ToForm[#2]] + 2;
        {#1, If[(ll += l) > 70, ll = l; ",\n  ", ", "], t} )&,
      t, {v} ],
    ";\n" }
]


Attributes[FormId] = Attributes[FormSq] = {Listable}

FormId[_[lhs_, rhs_], endl_:";\n"] :=
  {"id ", ToForm[lhs], " = ", ToForm[rhs], endl}

FormSq[_[lhs_, h_[x___]], endl_:";\n"] :=
  {"id ", #1, "^2 = ", #2, endl,
   "id ", #1, "^-2 = 1/(", #2, ")", endl}&[ ToForm[lhs], ToForm[h[x]] ] /;
  Context[h] === "System`"

FormSq[_[lhs_, rhs_], endl_:";\n"] :=
  {"id ", #1, "^2 = ", #2, endl,
   "id ", #1, "^-2 = ", #2, "^-1", endl}&[ ToForm[lhs], ToForm[rhs] ]


NCFuncs = Spinor | DottedSpinor

DeclareVars[expr_, inds_:{}, ranges_:{}, cvars_:{},
  vecs_:FormVectors, tens_:FormTensors] :=
Block[ {theexpr, vars, func},
  theexpr = {expr, FormSymbols};

  vars = Complement[Symbols[theexpr],
    cvars, inds, FormVectors, FormTensors, {i$$}];

  func = Complement[
    Cases[Head/@ Level[theexpr, {1, -2}, Heads -> True], _Symbol],
    cvars, FormVectors, FormTensors, ExtWF,
    {ga, MatrixTrace, FermionChain, NonCommutativeMultiply,
     Rule, Equal, Plus, Times, Power, Dot, Rational} ];

  { { FormDecl["i ", Replace[inds, ranges, 1]],
      FormDecl["s ", vars],
      FormDecl["cf ", DeleteCases[func, NCFuncs]],
      FormDecl["f ", Cases[func, NCFuncs]],
      FormDecl["v ", vecs],
      FormDecl["t ", tens] },
    inds, vars, func }
]


tempnum = 1

OpenForm[] :=
Block[ {hh},
  While[
    tempfile = ToFileName[Directory[],
      "fc" <> ToString[tempnum] <> ".frm"];
    FileType[tempfile] =!= None, ++tempnum ];
  Print[""];
  Print["preparing FORM code in ", tempfile];
  hh = OpenForm[tempfile];
  WriteString[hh, FormSetup];
  hh
]

OpenForm[file_] := OpenWrite[toform <> "\"" <> file <> "\"",
  FormatType -> InputForm, PageWidth -> 73]

toform =
  "!\"" <> ToFileName[{$FormCalcDir, $SystemID}, "ToForm"] <> "\" > "


Attributes[FormExpr] = {HoldAll}

FormExpr[x__] := Block[{addM = Plus, mulM = Times, powM = Power}, {x}]


FormEval[s_] := ToString[ToExpression[s], InputForm]

FormEvalDecl[s_] :=
  (DeclareVars[#,
      Complement[ Cases[#, SumOver[i_, ___] -> i, Infinity], vars[[2]] ], {},
      Flatten @ {Rest[vars], List@@ allint, D, Dminus4},
      {}, {}][[1]] <>
    ToString[# /. FinalFormRules, InputForm])& @ ToExpression[s]


Attributes[FormExec] = {HoldAll}

FormExec[s__Set] := Block[{s},
  ReadForm[$FormCmd, $FormCalcProgramDir, tempfile] ]


RunForm[r___][edit_, retain_] :=
Block[ {res},
  If[ edit, Pause[1]; Run[StringForm[$Editor, tempfile]]; Pause[3] ];
  WriteString["stdout", "running FORM... "];
  If[ Check[
        res = Set@@@
          Level[{r, FromFormRules, {Dot -> pSS}}, {2}, FormExec]; True,
        False] &&
    !retain, DeleteFile[tempfile] ];
  Print["ok"];
  res
]


(* things to do when the amplitude comes back from FORM *)

d$$ = MetricTensor

_dummy$$ = 1

FormSimplify = Simplify

mulM[args__] := FormSimplify[Times[args]] /. Plus -> addM


paveM[A0i[0], args__] := A0[args];
paveM[A0i[0, 0], args__] := A00[args];
paveM[h_[i__], args__] := h[paveid[h, i], args]

paveid[h_, i__] := paveid[h, i] =
Block[ {t = ToLowerCase[StringTake[ToString[h], 1]]},
  ToSymbol["LoopTools`", t, t, i]
]

A0[0] = 0

B0i[id:bb0 | dbb0, p_, m1_, m2_] :=
  B0i[id, p, m2, m1] /; !OrderedQ[{m1, m2}]

(Derivative[0, 1, 0, 0][B0i][#1, args__] := B0i[#2, args])&@@@
  {{bb0, dbb0}, {bb1, dbb1}, {bb00, dbb00}, {bb11, dbb11}}


(* abbreviationing business *)

Mat[0] = 0

fmeM[x_] := ferm[x, CurrentScale]

ferm[x__] := ferm[x] = Unique["F"]


pSS[a_, x_^n_] := pSS[a, x]^n
	(* different operator priority in FORM *)

pSS[x__] := pair[Pair[x], CurrentScale]

pair[x__] := pair[x] = Unique["Pair"]


e$$[x__] := eps[Eps[x], CurrentScale]

eps[x__] := eps[x] = Unique["Eps"]


abbM[p_Plus] := abbsum[abbM/@ p]

abbM[n_?NumberQ x_] := n abbM[x]

abbM[x_?AtomQ] = x

abbM[x_] := abbM[x] = Unique["Abb"]

abbsum[x_] := abbsum[x] = Unique["AbbSum"]


sunM[x_] := sunM[x] = Unique["SUN"]


numM[x_] := Sequence[qcount[x], num[Num[x]]]

num[x_] := num[x] = Unique["Num"]


qcount[p_Plus] := Max[qcount/@ List@@ p]

qcount[x_] := Count[x, q1, {-1}]


qcM[n_?NumberQ x_] := n qcM[x]

qcM[p_Plus] := -qcM[-p] /; MatchQ[p[[1]], _?Negative _.]

qcM[x_] := qcM[x] = Unique["QC"]


Attributes[dv] = {HoldAll}

dv[x_, _] := {} /; !FreeQ[x, Pattern]

dv[_[_[x_]], s_Symbol] := s -> x

dv[_[_[x_, 1]], s_Symbol] := s -> x

dv[_[_[x_, scale_]], s_Symbol] := s -> x/scale^(Count[4711 x, _k, {2}]/2)


Abbr[] := Flatten @ Apply[dv,
  DownValues/@ {ferm, pair, eps, abbM, abbsum, sunM, num, qcM}, {2}]

Abbr[patt__] :=
Block[ {all = Abbr[], need = Flatten[{patt}], omit},
  omit = Cases[need, !x_ -> x];
  need = DeleteCases[need, !_];
  FixedPoint[ abbsel[First/@ #]&, abbsel[need, omit] ]
]

abbsel[{need__}, {omit__}] := Select[all,
  !FreeQ[#, Alternatives[need]] && FreeQ[#, Alternatives[omit]]&]

abbsel[{need__}, ___] := Select[all, !FreeQ[#, Alternatives[need]]&]


IndexHeader[h_, expr_ /; Depth[expr] > 2] :=
Block[ {ind = Symbols[Level[expr, {-2}]]},
  ind = Select[ Union[ind], Head[Dim[#]] =!= Dim & ];
  h@@ ind /; Length[ind] =!= 0
]

IndexHeader[h_, _] = h


MomEncode[other_] := other /; FreeQ[other, k]

MomEncode[p_Plus] := MomEncode/@ p

MomEncode[f_. k[i_]] := K[f, i]


(* intabb introduces abbreviations for the loop integrals.
   They fall into two categories:
   1. A0, A00 (cint..),
   2. B0i, C0i, D0i, E0i, F0i (iint..).
   For the latter the LoopTools functions [BCDEF]get can be used to
   compute all tensor coefficients at once (which is much more efficient).
   Unlike the other integrals, whose results are double complex numbers,
   Cget and Dget return an integer pointing into a cache array.  In the
   conventions of LoopTools 2, the actual tensor coefficients are
   retrieved from the arrays [BCDEF]val. *)

( intabb[#1[i_, args__]] :=
  Block[ {abb = IndexHeader[ToSymbol["iint", ++ic], {args}]},
    iint = {iint, abb -> #2[args]};
    abbint[#1[id_, args]] = #3[id, abb];
    #3[i, abb]
  ] )&@@@
{ {B0i, Bget, Bval},
  {C0i, Cget, Cval},
  {D0i, Dget, Dval},
  {E0i, Eget, Eval},
  {F0i, Fget, Fval} }

intabb[func_] :=
Block[ {abb = IndexHeader[ToSymbol["cint", ++cc], func]},
  cint = {cint, abb -> MomEncode/@ func};
  abbint[func] = abb
]


ExtractInt[expr_] :=
Block[ {abbint, new, cint = {}, cc = 0, iint = {}, ic = 0},
  abbint[f_] := intabb[f];
  new = expr /. int:allint[__] :> abbint[int];
  {Flatten[cint], Flatten[iint], new}
]


GenNames[amp_] := amp /. {
  g:G[_][_][__][__] :> coup[g],
  m_Mass :> mass[m],
  x_GaugeXi :> xi[x],
  v:VertexFunction[_][__] :> vf[v] }

coup[g_] := coup[g] = Unique["Coupling"]

mass[m_] := mass[m] = Unique["Mass"]

xi[x_] := xi[x] = Unique["GaugeXi"]

vf[v_] := vf[v] = Unique["VertexFunction"]


GenericList[] := Flatten @
  Apply[dv, DownValues/@ {coup, mass, xi, vertex}, {2}]


Options[ClearProcess] = {
  ZapFunction -> Remove
}

ClearProcess[opt___Rule] :=
Block[ {zapf},
  {zapf} = ParseOpt[ClearProcess, opt];
  CurrentProcess = CurrentOptions = {};
  ReadFormClear[];
  Apply[zap, DownValues/@ {ferm, pair, eps, abbM, abbsum, sunM,
    num, qcM, coup, mass, xi, vertex}, {2}];
]


Attributes[zap] = {HoldAll}

zap[p_, s_Symbol] := (Unset[p]; zapf[s]) /; FreeQ[p, Pattern]

zap[p_, s_Symbol[__]] := (Unset[p]; zapf[s]) /; FreeQ[p, Pattern]


DenyFunc[leaf_, _[]] := LeafCount[#] < leaf &

DenyFunc[leaf_, deny_] := LeafCount[#] < leaf || !FreeQ[#, deny] &


Options[Abbreviate] = {
  MinLeafCount -> 10,
  Deny -> {k, q1},
  Preprocess -> Identity }

Abbreviate[expr_, obj_:2, opt___Rule] :=
Block[ {minleaf, deny, pre},
  {minleaf, deny, pre} = ParseOpt[Abbreviate, opt];
  deny = DenyFunc[minleaf, Alt[deny]];
  abbrev[expr, obj]
]

abbrev[expr_, lev_Integer] :=
Block[ {abbprep, sums, ind, pos},
  prepdef[pre];
  sums = Position[expr, SumOver];
  If[ Length[sums] === 0,
    Replace[expr, Plus -> abbprep[], {lev, Infinity}, Heads -> True],
  (* else *)
    ind[{0}] = {};
    ind[{r___, _, 0}] := {ind[{r, 0}],
      Extract[expr, Cases[sums, {r, i_, 0} -> {r, i, 1}]]};
    pos = Position[expr, Plus, {lev, Infinity}, Heads -> True];
    ReplacePart[expr, (abbprep@@ Union[ind[#]//Flatten])&/@ pos,
      pos, Array[List, Length[pos]]]
  ]
]

abbrev[expr_, patt_] := iso[patt][expr]


prepdef[Identity] := abbprep = abbplus

prepdef[f_] :=
  abbprep[i___][args__] := f[Plus[args]] /. Plus -> abbplus[i]


abbplus[___][args__] := Plus[args] /; deny[{args}]

abbplus[i___][x:(_Integer | _Rational)?Negative _., r__] :=
  -FromPlus[abbplus[i], -(x + r)]

abbplus[i___][a_^2, -b_^2] := abbplus[i][a, -b] abbplus[i][a, b]

abbplus[i___][args__] :=
  (abbsub@@ Select[{i}, !FreeQ[{args}, #]&])[Plus[args]]


iso[s_, ___][s_] = s

iso[s_, i___][x_] :=
  (abbsub@@ Select[{i}, !FreeQ[x, #]&])[pre[x]] /; FreeQ[x, s]

iso[s__][SumOver[i_, r___] x_] := SumOver[i, r] iso[s, i][x]

iso[s__][h_[x_]] := h[iso[s][x]]

iso[s_, i___][h_[x__]] := h[
  iso[s, i][Select[h[x], !FreeQ[#, s]&], h],
  abbplus[i][Select[h[x], FreeQ[#, s]&]]
] /; Length[Intersection[Attributes[h], {Flat, Orderless}]] === 2

iso[s__][x_] := iso[s]/@ x

iso[s__][h_[x___], h_] := iso[s]/@ h[x]

iso[s__][x_, _] := iso[s][x]


abbsub[][x_] := abbsub[][x] = Unique[$SubPrefix]

abbsub[i__][x_] := abbsub[i][x] = Unique[$SubPrefix][i]

$SubPrefix = "Sub"


Subexpr[] := Cases[ SubValues[abbsub],
  _[_[_[x_]], s_] :> (s -> x) /; FreeQ[x, Pattern] ]


Options[ClearSubexpr] = {
  ZapFunction -> Remove
}

ClearSubexpr[opt___Rule] :=
Block[ {zapf},
  {zapf} = ParseOpt[ClearProcess, opt];
  zap@@@ SubValues[abbsub];
]


likeQ[patt_][s_] := StringMatchQ[ToString[s], patt]


General::unknown = "Don't know how to register ``."

General::conflict = "Definition of `` conflicts with ``."

Attributes[RegisterAbbr] = {Listable}

RegisterAbbr[s_ -> x_Num] := regset[num, x, s]

RegisterAbbr[s_ -> other_] := regabb[s -> 4711 other]

RegisterAbbr[other_] := (Message[RegisterAbbr::unknown, other]; other)

regabb[s_ -> x_. p___Pair e___Eps d___DiracChain w___WeylChain t___eT tc___eTc] :=
  regabb[s, x/4711, Times[p], Times[e], d w, t tc]

regabb[s_, x_, 1, 1, 1, t_] := regabb[s, x t]

  regabb[s_?(likeQ["SUN*"]), x_] := regset[sunM, x, s]

  regabb[s_?(likeQ["AbbSum*"]), x_] := regset[abbsum, x, s]

  regabb[s_?(likeQ["Abb*"]), x_] := regset[abbM, x, s]

  regabb[s_?(likeQ["QC*"]), x_] := regset[qcM, x, s]

regabb[s_, x_, p_, 1, 1, 1] := regabb[s, x, p, pair]

regabb[s_, x_, 1, e_, 1, 1] := regabb[s, x, e, eps]

regabb[s_, x_, p_, e_, f_, t_] := regabb[s, x, p e f t, ferm]

  regabb[s_, 1, t_, h_] := regset[h, t, 1, s]

  regabb[s_, x_, t_, h_] := regset[h, t, x^(-2 Count[4711 t, _k, {2}]), s]

regabb[s_, x__] :=
  (Message[RegisterAbbr::unknown, #]; #)&[ s -> Times[x] ]


regset[h_, args__, s_] :=
  newset[h[args], Cases[DownValues[h], _[_[_[args]], x_] -> x], s]

Attributes[newset] = {HoldFirst}

newset[x_, {}, s_] := x = s

newset[x_, {s_}, s_] = s

newset[x_, {s1_}, s2_] :=
  (Message[RegisterAbbr::conflict, s1, s2]; s1)


Attributes[RegisterSubexpr] = {Listable}

RegisterSubexpr[s_Symbol[i__] -> x_] :=
  regset[abbsub[i], x, s]

RegisterSubexpr[s_Symbol -> x_] :=
  regset[abbsub[], x, s]

RegisterSubexpr[other_] :=
  (Message[RegisterSubexpr::unknown, other]; other)


$OptPrefix = "Opt"

newvar[] := Unique[$OptPrefix]

newvar[i__] := Unique[$OptPrefix][i]


Attributes[set] = {Flat, Orderless}


ptdef[0][x_, t_] := lotdef[x, lot@@ t]

ptdef[1][x_, p_] := (lotdef[-x, FromPlus[lot, -p]]; lotdef[x, lot@@ p])

lotdef[x_, l_lot] := (l = x; {})

lotdef[x_, _Integer x_] := x -> 0

lotdef[x_, other_] := x -> other


lotget[_, _Times] := {{}, {}}

lotget[_[_[terms__]], x_] := {x, set[terms]}


Attributes[ptcom] = {Listable}

ptcom[a_set, b_set, 0] := Intersection[a, b]

ptcom[a_set, b_set, 1] :=
Block[ {cp, cm},
  cp = Intersection[a, b];
  If[ Length[a] > 2 Length[cp],
    cm = Intersection[a, Thread[-b, set]];
    If[ Length[cm] > Length[cp], cp = cm ] ];
  cp
]

_ptcom = set[]


ptrul[a_, b_, 0] := a -> b

ptrul[a_, b_, 1] := {a -> b, Thread[-a, set] -> -b}


Overlap[] = set[]

Overlap[x__] := Intersection[x]


CSE[_, {}] = {}

CSE[pt_, abbr_] :=
Block[ {lot, alias, var, def, com, tmp, new = {}},
  Attributes[lot] = {Flat, Orderless};
  alias = Flatten[ ptdef[pt]@@@
    Sort[abbr, Length[ #1[[2]] ] <= Length[ #2[[2]] ] &] ];
  {var, def} = ToCat[2, Apply[lotget, DownValues[lot], 1]];
  def = def /. alias;

  Do[
    While[ Length[com = ptcom[def[[i]], def[[i + 1]], pt]] > 3,
      tmp = Ceiling[Length[com]/2];
      tmp = Overlap@@ Select[
        ptcom[com, Drop[def, {i, i + 1}], pt],
        Length[#] > tmp & ];
      If[ Length[tmp] > 3, com = tmp ];
      tmp = Position[def, com, 1, 1, Heads -> False];
      If[ Length[tmp] > 0,
        tmp = tmp[[1, 1]];
        def[[tmp, 0]] = List;
        def = def /. ptrul[com, var[[tmp]], pt];
        def[[tmp, 0]] = set,
      (* else *)
        tmp = newvar@@ Select[
          Union[Level[{var[[i]], var[[i + 1]]}, {2}]],
          !FreeQ[com, #]& ];
        new = {new, tmp -> com};
        def = def /. ptrul[com, tmp, pt] ]
    ],
  {i, Length[def] - 1}];

  Flatten[{new, Thread[var -> def], alias}]
]


AbbrCat[rul:_[_, _Plus]] := {{}, {}, rul}

AbbrCat[rul:_[_, t_Times]] := {{}, rul, {}} /; FreeQ[t, DiracChain]

AbbrCat[rul_] := {rul, {}, {}}


OptimizeAbbr[{}] = {}

OptimizeAbbr[rul:{__Rule}] := Flatten[
  {#1, CSE[0, #2] /. set -> Times, CSE[1, #3] /. set -> Plus}&@@
    ToCat[3, AbbrCat/@ rul] ]


Attributes[simple] = {Listable}

simple[r:(_ -> _?NumberQ)] := subst[r] = Sequence[]

simple[r:(_ -> n_. _Symbol)] := (subst[r] = Sequence[]) /; NumberQ[n]

simple[FortranExpr[vars__, expr_]] := FortranExpr[vars, simple[expr]]

simple[other_] = other

SubstSimpleAbbr[arg_] :=
Block[ {subst, new = arg, rul},
  While[ ( new = simple[new];
           Length[rul = (#[[1, 1, 1]]&)/@ DownValues[subst]] > 0 ),
    new = new /. FortranExpr[var_, tmpvar_, expr_] :>
      (FortranExpr[DeleteCases[var, #], DeleteCases[tmpvar, #], expr]&[
        First/@ rul ]) //. rul;
    Clear[subst]
  ];
  new
]


(* UV and IR finiteness checks *)

loopint = A0 | A00 | B0 | B1 | B00 | B11 | B001 | B111 |
  DB0 | DB1 | DB00 | DB11 |
  A0i | B0i | C0i | D0i | E0i | F0i

cutint = Acut | Bcut | Ccut | Dcut | Ecut | Fcut

allint = Join[loopint, cutint]

ToNewBRules = {
  B0[args__] -> B0i[bb0, args],
  B1[args__] -> B0i[bb1, args],
  B00[args__] -> B0i[bb00, args],
  B11[args__] -> B0i[bb11, args],
  B001[args__] -> B0i[bb001, args],
  B111[args__] -> B0i[bb111, args],
  DB0[args__] -> B0i[dbb0, args],
  DB1[args__] -> B0i[dbb1, args],
  DB00[args__] -> B0i[dbb00, args],
  DB11[args__] -> B0i[dbb11, args],
  C0[args__] -> C0i[cc0, args],
  D0[args__] -> D0i[dd0, args],
  E0[args__] -> E0i[ee0, args],
  F0[args__] -> F0i[ff0, args] }

ToOldBRules = {
  B0i[bb0, args__] -> B0[args],
  B0i[bb1, args__] -> B1[args],
  B0i[bb00, args__] -> B00[args],
  B0i[bb11, args__] -> B11[args],
  B0i[bb001, args__] -> B001[args],
  B0i[bb111, args__] -> B111[args],
  B0i[dbb0, args__] -> DB0[args],
  B0i[dbb1, args__] -> DB1[args],
  B0i[dbb00, args__] -> DB00[args],
  B0i[dbb11, args__] -> DB11[args] }


Attributes[FindDiv] = {Listable}

FindDiv[rc_ -> x_] := rc -> FindDiv[x]

FindDiv[i_IndexIf] := MapIf[FindDiv, i]

FindDiv[a:Amp[_][__]] := FindDiv/@ a

FindDiv[other_] := foo[other]


SubDiv[p_Plus] := SubDiv/@ p

SubDiv[x_ r_] := x SubDiv[r] /; FreeQ[x, div]

SubDiv[x_] := x - (x /. div -> 0)


SerDiv[x_] := Series[x /. Divergence -> -2/Dminus4, {Dminus4, 0, 0}]


MapDiv[f_, expr_, fin_] :=
Block[ {div = RCPattern[Divergence], foo = f},
  FindDiv[expr /. ToNewBRules /.
    {int:loopint[__] :> UVDiv[int] + fin int, D -> Dminus4 + 4} /.
    Dminus4^(n_?Negative) -> (-2/Divergence)^n]
]


UVDivergentPart[expr_] := MapDiv[SubDiv, expr, 0]


UVSeries[expr_, pow_] := Coefficient[UVSeries[expr], Dminus4, pow]

UVSeries[expr_] := MapDiv[SerDiv, expr, 1]


UVDiv[A0[m_]] = m Divergence

UVDiv[A00[m_]] = m^2/4 Divergence

UVDiv[B0i[bb0, __]] = Divergence

UVDiv[B0i[bb1, __]] = -1/2 Divergence

UVDiv[B0i[bb00, p_, m1_, m2_]] := -(p - 3 m1 - 3 m2)/12 Divergence

UVDiv[B0i[bb11, __]] = 1/3 Divergence

UVDiv[B0i[bb001, p_, m1_, m2_]] := (p - 2 m1 - 4 m2)/24 Divergence

UVDiv[B0i[bb111, __]] = -1/4 Divergence

UVDiv[B0i[dbb00, __]] = -1/12 Divergence

UVDiv[C0i[cc00, __]] = 1/4 Divergence

UVDiv[C0i[cc001 | cc002, __]] = -1/12 Divergence

UVDiv[C0i[cc0000, p1_, p2_, p1p2_, m1_, m2_, m3_]] :=
  -(p1 + p2 + p1p2 - 4 (m1 + m2 + m3))/96 Divergence

UVDiv[C0i[cc0011 | cc0022, __]] = 1/24 Divergence

UVDiv[C0i[cc0012, __]] = 1/48 Divergence

UVDiv[D0i[dd0000, __]] = 1/24 Divergence

UVDiv[D0i[dd00001 | dd00002 | dd00003, __]] = -1/96 Divergence

_UVDiv = 0


(* FeynCalc compatibility functions *)

PaVe[i__Integer, {p__}, {m__}] :=
  ToExpression[#1 <> "0i"][ ToSymbol[#2, #2, Sort[{i}]], p, m ]&[
    FromCharacterCode[Length[{m}] + 64],
    FromCharacterCode[Length[{m}] + 96] ]


FeynCalcGet[mask___] :=
Block[ {Global`GraphName, Global`Momentum = Identity},
  _Global`GraphName = 0;
  Plus@@ Get/@ FileNames[mask] /. ep_Eps :> I ep
]


FeynCalcPut[expr_, file_] :=
Block[ {C0i, D0i, E0i, F0i, PaVe},
  C0i[cc0, args___] := C0[args];
  D0i[dd0, args___] := D0[args];
  E0i[ee0, args___] := E0[args];
  F0i[ff0, args___] := F0[args];
  C0i[i_, p__, m1_, m2_, m3_] := PaVe[
    Sequence@@ (ToExpression/@ Drop[Characters[ToString[i]], 2]),
    {p}, {m1, m2, m3} ];
  D0i[i_, p__, m1_, m2_, m3_, m4_] := PaVe[
    Sequence@@ (ToExpression/@ Drop[Characters[ToString[i]], 2]),
    {p}, {m1, m2, m3, m4} ];
  E0i[i_, p__, m1_, m2_, m3_, m4_, m5_] := PaVe[
    Sequence@@ (ToExpression/@ Drop[Characters[ToString[i]], 2]),
    {p}, {m1, m2, m3, m4, m5} ];
  F0i[i_, p__, m1_, m2_, m3_, m4_, m5_, m6_] := PaVe[
    Sequence@@ (ToExpression/@ Drop[Characters[ToString[i]], 2]),
    {p}, {m1, m2, m3, m4, m5, m6} ];
  Put[expr /. ToOldBRules, file]
]


(* helicity matrix elements *)

ToTrace[fi_ -> plain_, fj_ -> conj_] := Mat[fi, fj] ->
  plain (conj /. DiracChain -> ConjChain /. ep_Eps -> -ep /. ConjWF)

ConjChain[s1_Spinor, om_Integer, g___, s2_Spinor] :=
  (#1 Reverse[DiracChain[s1, g, #2, s2]])&@@
     ConjOm[EvenQ[Length[{g}]], om]

ConjOm[False, om_] := {1, om};
ConjOm[_, 6] = {1, 7};
ConjOm[_, -6] = {1, -7};
ConjOm[_, 7] = {1, 6};
ConjOm[_, -7] = {1, -6};
ConjOm[_, 5] = {-1, 5};
ConjOm[_, -5] = {-1, -5};
ConjOm[_, 1] = {1, 1}
ConjOm[_, -1] = {1, -1}


HelTab[Hel[i_]] := helM[i, e]

HelTab[h_] := h + e


SelectAbbr[abbr_, All] = abbr

SelectAbbr[abbr_, expr_] := Select[abbr, !FreeQ[ expr, #[[1]] ]&]

SelectAbbr[abbr_, plain_, _[0]] :=
  {#, #}& @ SelectAbbr[abbr, plain]

SelectAbbr[abbr_, plain_, conj_] :=
  {Union[SelectAbbr[abbr, plain], #], #}& @ SelectAbbr[abbr, conj]

SelectAbbr[abbr_, patt_, plain_, conj_] :=
  SelectAbbr[Select[abbr, !FreeQ[#, patt]&], plain, conj]


Rep[Spinor[k[i_], _, s_]] :=
  {FromCharacterCode[(235 - s)/2], i}

Rep[DottedSpinor[k[i_], _, s_]] :=
  {Overscript[FromCharacterCode[(235 - s)/2], "."], i}

Format[DiracChain[s1_Spinor, om_, g___, s2_Spinor]] :=
  SequenceForm@@ Flatten[{"<", Rep[s1], "|", om,
    {",", #}&/@ {g}, "|", Rep[s2], ">"}]

Format[WeylChain[s1:(_Spinor | _DottedSpinor), om_, g___,
         s2:(_Spinor | _DottedSpinor)]] :=
  SequenceForm@@ Flatten[{"(", Rep[s1], "|", om,
    {",", #}&/@ {g}, "|", Rep[s2], ")"}]


Options[HelicityME] = {
  Source :> Abbr[],
  EditCode -> False,
  RetainFile -> False }

HelicityME::noprocess = "No process defined so far.  \
HelicityME works only after DeclareProcess or CalcFeynAmp."

HelicityME::weyl = "Warning: HelicityME does not work on WeylChains.  \
CalcFeynAmp uses DiracChains with the option FermionChains -> Chiral or VA."

General::nomat = "Warning: No matrix elements to compute."

HelicityME[plain_, opt___?OptionQ] := HelicityME[plain, plain, opt]

HelicityME[plain_, conj_, opt___?OptionQ] :=
Block[ {abbr, edit, retain,
part, ind, c = 0, vars, hels, helM, hh, e, mat},

  If[ CurrentProcess === {},
    Message[HelicityME::noprocess];
    Abort[] ];

  {abbr, edit, retain} = ParseOpt[HelicityME, opt];

  If[ !FreeQ[abbr, WeylChain], Message[HelicityME::weyl] ];

  abbr = SelectAbbr[abbr, DiracChain[_Spinor, ___], plain, conj];
  If[ Times@@ Length/@ abbr === 0,
    Message[HelicityME::nomat];
    Return[{}] ];

  part = Cases[abbr, Spinor[k[i_], __] -> i, Infinity] //Union;

  ind = Map[# -> "N" <> ToString[++c] <> "_?" &,
    Union[Cases[#, _Lor, Infinity]]&/@ abbr, {2}];

  mat = Flatten[Outer[ ToTrace,
      abbr[[1]] /. ind[[1]],
      abbr[[2]] /. ind[[2]] ]] /.
    Reverse/@ FromFormRules /. Eps -> "e_" /.
    FinalFormRules;

  Print["> ", Length[mat], " helicity matrix elements"];

  hels = HelTab/@ Hel/@ part;
  vars = DeclareVars[{Last/@ mat, hels}];

  hh = OpenForm[];
  WriteString[hh,
    vars[[1]] <> "\n\
table HEL(" <> ToString[Min[part]] <> ":" <>
               ToString[Max[part]] <> ", e?);\n" <>
    MapThread[{"fill HEL(", ToForm[#1], ") = ", ToForm[#2], ";\n"}&,
      {part, hels}] <> "\n" <>
    FormProcs <> "\
#include HelicityME.frm\n\n"];

  ( Write[hh, "L " <> "Mat" <> ToString/@ List@@ #1 <> " = ", #2, ";"];
    WriteString[hh, "\n#call Emit\n\n"]
  )&@@@ mat;

  WriteString[hh, ".end\n"];
  Close[hh];

  (e[#] = s[#])&/@ part;
  helM[i_, s_] := Hel[i] + s;

  Thread[First/@ mat -> Plus@@@ RunForm[][edit, retain]]
]


(* colour matrix elements *)

sunT[i_Symbol, i_] := SUNN

sunT[_, i_Symbol, i_] = 0

sunT[a_Integer, a_, i_Symbol, i_] = 1/2

sunT[_Integer, _Integer, i_Symbol, i_] = 0

sunT[t1___, a_Symbol, t2___, a_, t3___, i_, j_] :=
  (sunT[t1, t3, i, j] sunTr[t2] -
    sunT[t1, t2, t3, i, j]/SUNN)/2

sunT[t1___, a_Symbol, t2___, i_, j_] sunT[t3___, a_, t4___, k_, l_] ^:=
  (sunT[t1, t4, i, l] sunT[t3, t2, k, j] -
    sunT[t1, t2, i, j] sunT[t3, t4, k, l]/SUNN)/2

sunT[a___, i_, j_Symbol] sunT[b___, j_, k_] ^:=
  sunT[a, b, i, k]

sunT[a___, i_, j_Symbol] sunT[b___, k_, j_] ^:=
  Level[{{a}, Reverse[{b}], {i, k}}, {2}, sunT]

sunT[a___, j_Symbol, i_] sunT[b___, j_, k_] ^:=
  Level[{Reverse[{a}], {b}, {i, k}}, {2}, sunT]

sunT/: sunT[a___, i_, j_Symbol]^2 :=
  Level[{{a}, Reverse[{a}], {i, i}}, {2}, sunT]

sunT/: sunT[a___, i_Symbol, j_]^2 :=
  Level[{Reverse[{a}], {a, j, j}}, {2}, sunT]


sunTr[] := SUNN

sunTr[_] = 0

sunTr[a__] := sunT[a, #, #]& @ Unique["col"]


(* we assume that structures of the form delta[a, a] indicate
   summations over external colour/gluon indices *)

sunText[i_Symbol, i_] := Sqrt[SUNN]

sunText[a_Symbol, a_, 0, 0] := Sqrt[SUNN^2 - 1]/2

sunText[a___, 0, 0] := sunTrace[a]

sunText[other__] := sunT[other]


sunF[a_, b_, c_] := 2 I (sunTrace[c, b, a] - sunTrace[a, b, c])

sunF[a__, b_, c_] := (sunF[a, #] sunF[#, b, c])& @ Unique["glu"]


sunEps/: sunEps[a___, _Symbol, b___]^2 :=
  (sunT[#1, #1] sunT[#2, #2] - sunT[#1, #2]^2)&[a, b]

sunEps[a___, i_Symbol, b___] sunEps[c___, i_, d___] ^:=
  (-1)^Length[{a, c}] *
  ((sunT[#1, #3] sunT[#2, #4] -
    sunT[#1, #2] sunT[#3, #4])&[a, b, c, d])


Conjugate[t_SUNT] ^:= RotateLeft[Reverse[t], 2]

Conjugate[f_SUNF] ^= f


ColourSimplify[plain_, conj_:1] :=
Block[ {res, ind},
  {res, ind} = Reap[ plain Conjugate[conj] /.
    {IndexDelta -> idelta, SumOver -> sumover} /.
    SUNSum[i_, _] :> (Sow[i]; 1) ];
  Simplify[ Expand[ res /.
    Cases[ind, i_Index :> i -> iname@@ i, {2}] /.
    {SUNT -> sunText, SUNF -> sunF, SUNEps -> sunEps} /.
    sunTrace[a__]^n_. :> Times@@ Table[sunTr[a], {n}]
  ] /. {sunT[i_, j_] :> Sort[SUNT[i, j]],
        sunT -> SUNT, sunEps -> SUNEps} ]
]


ColourGrouping[tops_] := DiagramGrouping[ tops,
  Replace[
    ColourSimplify[Times@@
      FeynAmpCases[_[Index[Colour | Gluon, _], ___]][##]],
    _?NumberQ r_ -> r]& ]


ColourFactor[fi_ -> plain_, fj_ -> conj_] :=
  Mat[fi, fj] -> ColourSimplify[plain, conj]


Options[ColourME] = {Source :> Abbr[]}

ColourME[plain_, opt___?OptionQ] := ColourME[plain, plain, opt]

ColourME[plain_, conj_, opt___?OptionQ] :=
Block[ {abbr},
  {abbr} = ParseOpt[ColourME, opt];

  abbr = SelectAbbr[abbr, SUNObjs, plain, conj];
  If[ Times@@ Length/@ abbr === 0,
    Message[ColourME::nomat];
    Return[{}] ];

  Outer[ ColourFactor, abbr[[1]], abbr[[2]] ] //Flatten
]


(* squaring the matrix element *)

UniquefyIndices[conj_, plain__] :=
Block[ {ind},
  ind = Intersection@@
    (Union[Cases[#, SumOver[x_, ___] -> x, Infinity]]&)/@ {conj, plain};
  conj /. Thread[ind -> (ToSymbol[#, "c"]&)/@ ind]
]


SquaredME[amp_] := SquaredME[amp, amp]

SquaredME[Amp[_][0], _] = 0

SquaredME[_, Amp[_][0]] = 0

SquaredME[Amp[_][plain__], Amp[_][conj__]] :=
  Plus[plain] Conjugate[UniquefyIndices[Plus[conj], plain]] /;
  FreeQ[{plain, conj}, Mat]

SquaredME[Amp[_][plain__], Amp[_][conj__]] :=
  Apply[Plus,
    Outer[ ToSquared,
      MatList[Plus[plain]], MatList[UniquefyIndices[Plus[conj], plain]] ],
    {0, 1}]


MatList[expr_] := ToList[Collect[expr, _Mat, Hold]]

ToList[p_Plus] := List@@ p

ToList[other_] := {other}


ToSquared[Mat[m1_] x1_., Mat[m2_] x2_.] :=
  ReleaseHold[x1 Conjugate[x2]] ToMat[m1, m2]

ToMat[m1_Symbol, m2_Symbol] := Mat[m1, m2]

ToMat[m1_, m2_] := Inner[Mat, m1, m2, Times]


Unprotect[Conjugate]

Format[ Conjugate[x_] ] := SequenceForm[x, Superscript["*"]]

(*
Format[ Conjugate[t_Times] ] :=
  SequenceForm["(", t, ")", Superscript["*"]]
*)

Conjugate[D] = D

Conjugate[Dminus4] = Dminus4

Re[Dminus4] ^= Dminus4

Re[Divergence] ^= Divergence

Conjugate[Divergence] ^= Divergence

Conjugate[o_SumOver] = o

Conjugate[d_IndexDelta] = d

Conjugate[p_Plus] := Conjugate/@ p

Conjugate[t_Times] := Conjugate/@ t

Conjugate[d_Den] := Conjugate/@ d

Conjugate[p_Pair] := Conjugate/@ p

Conjugate[ep_Eps] := -Conjugate/@ ep

Conjugate[e[n_]] := ec[n]

Conjugate[ec[n_]] := e[n]

Conjugate[z[n_]] := zc[n]

Conjugate[zc[n_]] := z[n]

Conjugate[eT[n_][i__]] := eTc[n][i]

Conjugate[eTc[n_][i__]] := eT[n][i]

Conjugate[k[n_]] := k[n]

Conjugate[s[n_]] := s[n]

Conjugate[x_?RealQ] = x

Conjugate[(x_?RealQ)^n_] = x^n

Protect[Conjugate]


Unprotect[Re]

Re[D] = D

Re[Dminus4] ^= Dminus4

Re[Divergence] ^= Divergence

Re[o_SumOver] = o

Re[d_IndexDelta] = d

Re[p_Plus] := Re/@ p

Re[d_Den] := Re/@ d

Re[x_?RealQ y_.] := x Re[y]

Re[(x_?RealQ)^n_ y_.] := x^n Re[y]

Protect[Re]


(* performing the polarization sum analytically *)

Options[PolarizationSum] = {
  Source :> Abbr[],
  GaugeTerms -> True,
  EditCode -> False,
  RetainFile -> False }

PolarizationSum::noprocess = "No process defined so far.  \
PolarizationSum works only after DeclareProcess or CalcFeynAmp."

PolarizationSum::incomp = "PolarizationSum used on an amplitude \
other than the last one set up by DeclareProcess or CalcFeynAmp."

PolarizationSum[amp:Amp[_][__].., ___] :=
  Message[PolarizationSum::incomp] /;
    {CurrentProcess} =!= Union[#[[0, 1]]&/@ {amp}]

PolarizationSum[amp:Amp[_][__].., opt___?OptionQ] :=
Block[ {Hel},
  _Hel = 0;
  PolarizationSum[
    SquaredME[amp] /.
      HelicityME[amp, FilterOptions[HelicityME, opt]] /.
      ColourME[amp, FilterOptions[ColourME, opt]],
    opt ]
]

PolarizationSum[expr_, opt___?OptionQ] :=
Block[ {abbr, gauge, edit, retain,
fullexpr, lor, indices, legs, masses, vars, hh},

  If[ CurrentProcess === {},
    Message[PolarizationSum::noprocess];
    Abort[] ];

  {abbr, gauge, edit, retain} = ParseOpt[PolarizationSum, opt];

  fullexpr = expr //. abbr /. FinalFormRules;
  lor = Cases[fullexpr, _Lor, Infinity] //Union;
  indices = FormIndices[[ Level[lor, {2}] ]];
  fullexpr = fullexpr /. Thread[lor -> indices];

  legs = Cases[fullexpr, Alt[ExtWF][i_] -> i,
    Infinity, Heads -> True] //Union;
  masses = Level[CurrentProcess, {2}][[legs]];

  fullexpr = fullexpr /. Reverse/@ FromFormRules /.
    {Eps -> "e_", Pair -> Dot} /.
    NoExpandRule /.
    FinalFormRules;

  vars = DeclareVars[{fullexpr, masses}, indices];

  hh = OpenForm[];
  WriteString[hh,
    vars[[1]] <> "\n\
#define GaugeTerms \"" <> ToBool[gauge] <> "\"\n" <>
    FormProcs <> "\
#include PolarizationSum.frm\n\n"];

  Write[hh, "L SquaredME = ", fullexpr, ";"];

  WriteString[hh,
    MapThread[{"\n#call PolSum(", ToString[#1], ", ", ToForm[#2], ")"}&,
      {legs, masses}] <>
    "\n\n#call Emit\n"];
  Close[hh];

  Plus@@ RunForm[][edit, retain][[1]]
]


FormIndices = Array[ToSymbol["FormCalc`Lor", #]&, 10]


(* set up a directory for the Fortran code *)

ChkExist[file__] := (MkDir[DirectoryName[#]]; #)& @ StringJoin[file]

(*ChkExist[dir_, file_] := (MkDir[dir]; ToFileName[dir, file])*)


MkDir[""] = "."

MkDir[dir_String] := dir /; FileType[dir] === Directory

MkDir[dir_String] := Check[CreateDirectory[dir], Abort[]]

MkDir[dirs_List] := Fold[MkDir[ToFileName[##]]&, {}, dirs]

MkDir[dirs__String] := MkDir[{dirs}]


Off[CopyFile::filex]

Options[SetupCodeDir] = {Drivers -> "drivers"}

SetupCodeDir[dir_, opt___Rule] :=
Block[ {drivers, path, files = {}},
  {drivers} = ParseOpt[SetupCodeDir, opt];
  path = SetDirectory[MkDir[dir]];
  ResetDirectory[];

  If[ FileType[drivers] === Directory,
    SetDirectory[drivers];
    CopyFile[#, ToFileName[path, #]]&/@ (files = FileNames[]);
    ResetDirectory[]
  ];

  If[ FileType[$DriversDir] === Directory,
    SetDirectory[$DriversDir];
    CopyFile[#, ToFileName[path, #]]&/@ Complement[FileNames[], files];
    ResetDirectory[]
  ];

  CopyFile[
    ToFileName[{$FormCalcDir, $SystemID}, "util.a"],
    ToFileName[path, "util.a"] ];

  path
]


(* Fortran code generation *)

Dim[i_Integer] = i


FortranNames[base_, ind___] := {
  #1 <> abridge[0, #2, Infinity],
  #1 <> abridge[0, #2, $MaxFortranName -
          StringLength[$SymbolPrefix <> #1 <> #2]]
}&[ ToString[base], ToString/@ Flatten[{ind}] ]

abridge[0, ind_, n_] :=
  Insert[ind, "_", Array[List, Min[Length[ind], n]]] /; n > 0

abridge[0, ind_, n_] := {"_", abridge[0, ind, {}, n - 1]}

abridge[_, {i___}, {j___}, 0] := {i, j}

abridge[f_, {i___, s_}, {j___}, n_] :=
  If[ StringLength[s] > 2 && !DigitQ[StringTake[s, {2}]],
    abridge[1, {i}, {StringDrop[s, {2}], j}, n + 1],
    abridge[f, {i}, {s, j}, n] ]

abridge[1, {}, ind_, n_] := abridge[0, ind, {}, n]

abridge[0, {}, ind_, _] = ind


Attributes[WriteFF] = {HoldFirst, Listable}

WriteFF[x_Symbol, array_] :=
  FFPut[x, array, Block[{x}, ToString[x]]]

WriteFF[amp_, array_] :=
  FFPut[amp, array, ToString[array] <> ToString[++modnum]]


InvDef[h_[i_], h_[j_]] := Invariant[1, i, j] -> SInvariant[i, j]

InvDef[_[i_], _[j_]] := Invariant[-1, i, j] -> TInvariant[i, j]


InvList[n_, r__] := {InvDef[n, #]&/@ {r}, InvList[r]}

InvList[_] = {}


ProcessCheck[p_] := (
  proc = p;
  name = ToString[p];
  legs = Plus@@ Length/@ p;
  invs = If[ legs < 4, {},
    Flatten[InvList@@
      MapIndexed[ Apply, Drop[Flatten[{1&/@ p[[1]], -1&/@ p[[2]]}], -1] ]] ];
)

ProcessCheck[p_, p_] = 0

_ProcessCheck := Message[WriteSquaredME::incomp]


FFPut[Amp[p_][amp__], array_, file_] := (
  ProcessCheck[p, proc];
  FFWrite[#, array, file]&/@ {amp}
)

FFPut[_[], __] = {}

FFPut[other_, __] := (Message[WriteSquaredME::noamp, other]; Abort[])


FFWrite[0, __] = {}

FFWrite[amp_, array_, file_] :=
Block[ {ind, ff, mods},
  ind = Cases[amp, SumOver[i_, r_] :> (Dim[i] = r; i), Infinity];
  ff = FFList[amp /. unused[array] -> 0 /. fcs /. xrules /. {
    _SumOver -> 1,
    int:allint[__] :> abbint[int] }, array];
  mods = FileSplit[ff, FortranNames[file, ind], FFMod];
  (Indices[#] = ind)&/@ mods;
  mods
]

FFMod[ff_, {fmod_, mod_}] :=
Block[ {file = prefix <> fmod <> ".F", hh},
  hh = OpenFortran[ModName[file]];
  WriteString[hh,
    Hdr[file, "form factors for " <> name] <>
    fincl <>
    SubroutineDecl[mod] <>
    "#include \"" <> prefix <> "vars.h\"\n\n"];
  WriteExpr[hh, ff, Optimize -> True, DebugLines -> fmod];
  WriteString[hh, "\tend\n"];
  Close[hh];
  mod
]


FFList[0, _] = {}

FFList[amp_, array_] :=
  ( maxmat[array] = {Mat[1]};
    RuleAdd[array[1], amp] ) /; FreeQ[amp, Mat]

FFList[amp_Plus, array_] := FFList[#, array]&/@ List@@ amp

FFList[Mat[m_] x_., array_] :=
  ( maxmat[array] = MaxDims[maxmat[array], Level[m, {-2}]];
    RuleAdd[Level[m, {-1}, array], x] )


(* Calculating the abbreviations in a clever way is key to a decent
   performance of the generated code.  Therefore, the abbreviations are
   split into three categories:
   1. objects that depend only on model constants and S
      -> subroutine abbr_s,
   2. objects that depend on other phase-space variables (angles etc.)
      -> subroutine abbr_angle,
   3. objects that depend on the helicities
      -> subroutine abbr_hel.
   The master subroutine SquaredME takes care to invoke these abbr_nnn
   subroutines only when necessary. *)

Category[rul_] := {{}, {}, rul} /;
  !FreeQ[rul, Hel | e | ec | Spinor | DottedSpinor | q1]

Category[rul_] := {{}, rul, {}} /; !FreeQ[rul, angledep]

Category[rul_] := {rul, {}, {}}


setdep[Tag[___, r_]] := setdep[r]

setdep[v_ -> _] := v = TAG[++c]

markdep[Tag[___, r_]] := markdep[r]

markdep[v_ -> x_] := (v = TAG[++c]; {}) /; !FreeQ[x, TAG]

markdep[r_] := {} /; !FreeQ[r, TAG]

MoveDepsRight[li_List] := {li}

MoveDepsRight[f__List, li_List] :=
Block[ {pos = {f}, c = 0, cc = 0},
  Block[ #,
    setdep/@ li;
    pos = Map[markdep, pos, {2}];
    While[c != cc, cc = c; pos = pos]
  ]&[ Union@@ Map[Kind, {f, li}, {2}] ];
  pos = Position[pos, {}, {2}, Heads -> False];
  Append[
    MoveDepsRight@@ Delete[{f}, pos],
    Flatten[{li, Extract[{f}, pos]}] ]
]


depcat[expr_][r:Tag[___, v_ -> _]] := (cpos = 1; {r, {}}) /; FreeQ[expr, v]

depcat[expr_][r:(v_ -> _)] := (cpos = 1; {r, {}}) /; FreeQ[expr, v]

depcat[_][r_] := (cnew = 1; {{}, r})

MoveDepsLeft[li_List] := {li}

MoveDepsLeft[li_List, f__List] :=
Block[ {pos = {f}, old = {}, new = li, cpos = 1, cnew = 1},
  While[ cpos + cnew > 1,
    old = {new, old};
    cpos = cnew = 0;
    {pos, new} = Transpose[ ToCat[2, depcat[new]/@ #1]&/@ pos ];
    new = Flatten[new] ];
  Prepend[MoveDepsLeft@@ pos, Flatten[{new, old}]]
]


toTicket[Tag[___, x_]] := toTicket[x]

toTicket[v_ -> x_] :=
  Ticket[v = Dep[v] /. HoldPattern[Pattern][a_, _] -> a, x]

toTicket[x_] := Ticket[x]

fromTicket[v_, x_] := v -> x

fromTicket[x_] = x

OnePassOrder::recurs = "Recursive definition among `` \
(full list in $OnePassDebug).  Returning list unordered."

OnePassOrder[li_List] :=
Block[ {c = 0, l = Length[li], Dep, Ticket, posmap, prev},
  Attributes[Dep] = {HoldFirst};
  Block[#, posmap = Hold@@ Evaluate[toTicket/@ li]]& @ Union[Kind/@ li];
  Ticket[v_, x_] := (v = Random[]; ++c) /; FreeQ[x, Dep];
  Ticket[x_] := ++c /; FreeQ[x, Dep];
  While[ c < l,
    prev = c;
    posmap = Evaluate/@ posmap;
    If[ c === prev,
      $OnePassDebug = Cases[posmap, t_Ticket :> fromTicket@@ t] /.
        Dep -> Identity;
      Message[OnePassOrder::recurs, $OnePassDebug];
      Return[li] ]
  ];
  li[[ Level[Sort[MapIndexed[List, posmap]], {3}] ]]
]

$OnePassDebug = {}


AbbrMod[{}, _] := Sequence[]

AbbrMod[abbr_, mod_] :=
Block[ {file = prefix <> mod <> ".F", hh},
  hh = OpenFortran[ModName[file]];
  WriteString[hh,
    Hdr[file, "abbreviations for " <> name] <>
    fincl <>
    SubroutineDecl[mod] <>
    "#include \"" <> prefix <> "vars.h\"\n\n"];
  WriteExpr[hh, abbr];
  WriteString[hh, "\tend\n"];
  Close[hh];
  mod
]


NumMod[{}, _] := Sequence[]

NumMod[expr_, mod_] :=
Block[ {file = prefix <> mod <> ".F", hh},
  hh = OpenFortran[ModName[file]];
  WriteString[hh,
    Hdr[file, "numerators for " <> name] <>
    fincl];
  WriteNum[hh, #]&/@ expr;
  Close[hh];
  mod
]


WriteNum[hh_, var_ -> Num[expr_]] :=
Block[ {Global`res, Global`q1in},
  WriteString[hh,
    SubroutineDecl[var[Global`q1in, Global`qt2, Global`res],
      "\tdouble complex q1in(0:3), qt2, res\n"] <>
    "#include \"" <> prefix <> "num.h\"\n\n"];
  WriteExpr[hh, Global`res -> (expr /. WeylChain -> SplitChain)];
  WriteString[hh, "\tend\n\n\n"];
]


$SymbolPrefix = ""


VarDecl[_[], __] = ""

VarDecl[vars_, type_, flag_String] :=
  "\n#ifndef " <> flag <> VarDecl[vars, type] <> "\n#endif"

VarDecl[vars_, type_] := StringJoin @ varDecl[
  DeleteCases[
    Replace[vars /. (x_ -> _) -> x, i_Symbol :> Dim[i], {2, Infinity}],
    _[0] | _[] ],
  type ]

varDecl[vars_, {type_String, common_String}] := {
  varDecl[vars, type],
  varDecl[kind/@ vars, "common /" <> $SymbolPrefix <> common <> "/"],
  "\n" }

varDecl[vars_, type_String] :=
Block[ {lmax = 63 - StringLength[type], llen = Infinity, vlen},
  ( llen += (vlen = StringLength[#] + 2);
    {If[llen > lmax, llen = vlen; {"\n\t", type}, ","], " ", #} )&/@
  ToFortran/@ vars
]


CommonDecl[vars_, type_, common_, flag___] :=
  VarDecl[vars, {type, common}, flag]


SubroutineDecl[name_, decl___String] := "\
\tsubroutine " <> $SymbolPrefix <> ToFortran[name] <> "\n\
\timplicit none\n" <> decl <> "\n"


CallDecl[li_List] := StringJoin[CallDecl/@ li]

CallDecl[DoLoop[name_, ind__]] :=
  ("\n" <> #1 <> CallDecl[name] <> #2 &)@@ DoDecl[ind]

CallDecl[name_] := "\tcall " <> $SymbolPrefix <> name <> "\n"


DoDecl[{var_}] := DoDecl[{var, Dim[var]}]

DoDecl[{_, _Dim}] := {{}, {}}

DoDecl[{var_, from_:1, to__}] := {
  "\tdo " <> ToFortran[var] <> " = " <> ToFortran[{from, to}] <> "\n",
  "\tenddo\n" }

DoDecl[var_] := DoDecl[{var, Dim[var]}]

DoDecl[vars__] := StringJoin/@ Transpose[DoDecl/@ Reverse[{vars}]]


	(* LoopComponents gives back e.g.
		1. {F[jFtree], SUN[jSUNtree]}
		2. "Ctree(jFtree,jSUNtree)"
		3. "Ctree, nFtree, nSUNtree"
		4. "\n\tinteger jFtree, nFtree"
		5. "\n\tparameter (nFtree = 5)"
		6. "\n\tdo jFtree = 1, nFtree
		    \n\tdo jSUNtree = 1, nSUNtree"
		7. "\n\tenddo
                    \n\tenddo" *)

LoopVar[h_[-1]] = {h[1], "", "", "", "", ""}

LoopVar[h_[n_]] :=
Block[ {v = ToString[h] <> type},
  { h[ToSymbol["j", v]],
    ", n" <> v,
    "\n\tinteger j" <> v <> ", n" <> v,
    "\n\tparameter (n" <> v <> " = " <> ToString[n] <> ")",
    "\n\tdo j" <> v <> " = 1, n" <> v,
    "\n\tenddo" }
]

LoopComponents[arr_, {}] = 0

LoopComponents[arr_, maxmat_] :=
Block[ {type = StringDrop[ToString[arr], 1]},
  {#1, ToFortran[Level[#1, {2}, arr]], ToString[arr] <> #2, ##3}&@@
  Transpose[ LoopVar/@ maxmat ]
]


Invoke[mod0_, {}] := CallDecl[mod0]

Invoke[mod0_, mod1_] := {
  CallDecl[mod0],
  "\tLOOP_IF\n",
  CallDecl[mod1],
  "\tLOOP_ENDIF\n" }


LoopReduce[m_] := Transpose[MapThread[LoopNeed, m]] /; SameQ@@ Length/@ m

LoopReduce[m_] := m /. 1 -> -1

LoopNeed[h_[1], h_[1]] := {h[-1], h[-1]}

LoopNeed[other__] := {other}


MatType[_Mat, _] = 1

MatType[h_[i_], h_[j_]] := MatType[h][i, j]

MatType[h_] := MatType[h] = ToSymbol["Mat", h]


Assort[m_Mat -> x_] := {m -> x, {}, {}, {}}

Assort[f_ -> r_. w_WeylChain] := {{}, {}, {}, f -> r SplitChain@@ w}

Assort[f_ -> x_] := {{}, f -> ToArray[f], {}, {}} /;
  !FreeQ[x, DiracChain | SUNT]

Assort[x_ -> n_Num] := {{}, {}, x -> n, {}}

Assort[other_] := {{}, {}, {}, other}


SplitChain[h1_[_[i_], _, s1_], om_, g___, h2_[_[j_], _, s2_]] :=
Block[ {om1 = 13 - om, om2 = 6 + Mod[om + Length[{g}], 2], o},
  o = om2;
  afsign[i, s1, om1] afsign[j, s2, om2] SxS[ h1[i, s1, om1],
    Fold[xs[o = 13 - o][#2, #1]&, h2[j, s2, om2], Reverse[{g}]] ]
]


afsign[i_, -1, 6] = -Hel[i]

afsign[i_, -1, 7] = Hel[i]

_afsign = 1


xs[6] = VxS

xs[7] = BxS

SxS[x_, _[-1, y_]] := SeS[x, y]

VxS[x_, _[-1, y_]] := VeS[x, y]

BxS[x_, _[-1, y_]] := BeS[x, y]


DefFilter[h_][m_[i_, j_]] :=
  h[_[m[x_, y_], _]] := x <= i && y <= j

DefCat[h_][m_[i_, j_]] :=
  h[r:_[m[x_, y_], _]] := If[x <= i && y <= j, {r, {}}, {{}, r}]


Attributes[Cond] = {HoldRest}

Cond[True, args__] := {args}

Cond[False, __] = {}


Attributes[WriteSquaredME] = {HoldAll}

Options[WriteSquaredME] = {
  TreeSquare -> True,
  LoopSquare -> False,
  Folder -> "squaredme",
  ExtraRules -> {},
  FilePrefix -> "",
  SymbolPrefix -> "",
  FileHeader -> "* %f\n* %d\n* generated by FormCalc " <>
    ToString[$FormCalc] <> " %t\n\n",
  FileIncludes -> "",
  SubroutineIncludes -> "#include \"decl.h\"\n" }

WriteSquaredME::incomp = "Warning: writing out Fortran code for \
incompatible processes."

WriteSquaredME::noamp = "`` is not an amplitude."

WriteSquaredME::empty = "Warning: no amplitudes were specified."

WriteSquaredME::badmat = "Incompatible matrix elements `` and ``."

WriteSquaredME[tree_, loop_, dir_, opt___Rule] :=
  WriteSquaredME[tree, loop, Abbr[], dir, opt]

WriteSquaredME[tree_, loop_, abbr__, dir_, opt___Rule] :=
Block[ {treesq, loopsq, folder, xrules, prefix, $SymbolPrefix,
header, fincl, sincl,
ModName, Hdr, proc = Sequence[], name, legs, invs,
mat, nums, fcs, abrs, matsel, treecat, angledep,
Dim, abbint, cint = {}, cc = 0, iint = {}, ic = 0,
mc = 0, Global`c,
Indices, Hel, hels, file, files, hh,
unused, maxmat, ntree, nloop, mats, com, loops,
ffmods, nummods, abbrmods},

  {treesq, loopsq, folder, xrules, prefix, $SymbolPrefix,
    header, fincl, sincl} =
    ParseOpt[WriteSquaredME, opt] /. Options[WriteRenConst];

  (ModName[mod_] := ModName[mod] = ToFileName[#, mod])& @
    MkDir[dir, folder];

  header = StringReplace[header, "%t" -> TimeStamp[]];
  Hdr[f_, d_] := StringReplace[header, {"%f" -> f, "%d" -> d}];

  abbint[f_] := intabb[f];

  {mat, fcs, nums, abrs} = ToCat[4, Assort/@ Flatten[{abbr}]];
  mats = First/@ DeleteCases[mat, _ -> 0];
  unused[Ctree] = Alt@@
    Select[Union[#[[1, 2]]&/@ mat], FreeQ[mats, Mat[_, #]]&];
  unused[Cloop] = Alt@@
    Select[Union[#[[1, 1]]&/@ mat], FreeQ[mats, Mat[#, _]]&];

(* Part 1: the form factors *)

  maxmat[_] = {};
  ffmods = Flatten/@ {
    Block[{modnum = 0}, WriteFF[tree, Ctree]],
    Block[{modnum = 0}, WriteFF[loop, Cloop]] };
  If[ Plus@@ Length/@ ffmods === 0,
    Message[WriteSquaredME::empty]; Return[{}] ];

  abrs = abrs /. xrules /. int:allint[__] :> abbint[int];
  cint = Flatten[cint];
  iint = Flatten[iint];

  ntree = maxmat[Ctree];
  nloop = maxmat[Cloop];
  If[ Length[ntree] === 0, ntree = nloop,
    If[ Length[nloop] === 0, nloop = ntree,
      If[ Head/@ ntree =!= Head/@ nloop,
        Message[WriteSquaredME::badmat, ntree, nloop];
        Abort[] ];
      nloop = MaxDims[nloop, ntree];
      If[ loopsq, ntree = nloop ];
  ] ];
  mats = Select[MapThread[MatType, {nloop, ntree}], Length[#] > 0 &];

  com = VarDecl[invs, {"double precision", "kinvars"}] <>
    VarDecl[{Hel[legs]}, {"integer", "kinvars"}] <>
    VarDecl[abrs, {"double complex", "abbrev"}] <>
    VarDecl[cint, {"double complex", "loopint"}] <>
    VarDecl[iint, {"integer", "loopint"}] <>
    VarDecl[
      #[[1, 1, 1]]&/@ DownValues[Dim],
      {"integer", "indices"} ] <>
    VarDecl[
      Flatten[{ mats,
        Level[maxmat[Ctree], {2}, Ctree],
        Level[maxmat[Cloop], {2}, Cloop] }],
      {"double complex", "formfactors"} ];

(* Part 2: the numerators and abbreviations *)

  Scan[DefFilter[matsel], mats];
  mat = Select[mat /. fcs /. Mat -> MatType, matsel];

  abrs = Flatten[{abrs, mat, cint, iint}];

  (* split into tree/loop *)
  Scan[DefCat[treecat][MatType[#, #]]&, maxmat[Ctree]];
  treecat[r:_[v_, _]] := {r, {}} /; !FreeQ[tree, v];
  treecat[r_] := {{}, r};
  abrs = Join[#1, Tag/@ #2]&@@
    MoveDepsLeft@@ ToCat[2, treecat/@ abrs];

  (* split into s/angle/hel *)
  angledep = Alt[(Range[#2] + #1)&@@ Length/@ proc];
  angledep = Alt[{
    Cases[invs, _[x_, r_] :> x /; MemberQ[r, angledep]],
    (k | s)[angledep],
    kcomp[_, angledep] }];
  abrs = MoveDepsRight@@ ToCat[3, Category/@ abrs];
  abrs = ToCat[2, #]&/@ Replace[abrs,
    {Tag[r_] -> {{}, r}, r_ -> {r, {}}}, {2} ];

  nummods = FileSplit[nums, "num", NumMod];
  nums = (#1 -> $SymbolPrefix <> ToFortran[#1] &)@@@ nums;
  abrs = abrs /. nums;

  abbrmods = MapThread[
    FileSplit[ToDoLoops[OnePassOrder[#1]], #2, AbbrMod]&,
    { abrs, {{"abbr0_s",     "abbr1_s"},
             {"abbr0_angle", "abbr1_angle"},
             {"abbr0_hel",   "abbr1_hel"}} }, 2 ];

(* Part 3: the variable declarations *)

  file = prefix <> "vars.h";
  hh = OpenWrite[ModName[file]];

  WriteString[hh,
    Hdr[file, "variable declarations"] <>
    sincl <>
    com <>
    VarDecl[Last/@ nums, "external", "__num_h__"] <> "\n"];

  Close[hh];

(* Part 4: the makefile *)

  hh = OpenWrite[ModName[prefix <> "makefile"]];

  WriteString[hh, "\
OBJS :=" <> ({" \\\n  $(DIR)/", #, ".o"}&)/@
    Flatten[{abbrmods, nummods, ffmods, prefix <> "SquaredME"}] <> "\n\n\
$(LIB): $(LIB)($(OBJS))\n\n\
$(LIB)($(OBJS)): $(DIR)/" <> prefix <> "vars.h $(DECL_H)\n\n\
LIBS += $(LIB)\n\n"];

  Close[hh];

(* Part 5: the master subroutine SquaredME *)

  hels = Array["Hel" <> ToString[#] -> ToString[#]&, legs];
  {maxmat[Ctree], maxmat[Cloop]} = LoopReduce[{maxmat[Ctree], maxmat[Cloop]}];
  ntree = LoopComponents[Ctree, maxmat[Ctree]];
  nloop = LoopComponents[Cloop, maxmat[Cloop]];
  loops = DeleteCases[{ntree, nloop}, 0];

  file = prefix <> "SquaredME.F";
  hh = OpenFortran[ModName[file]];

  (* a) declarations *)
  WriteString[hh, "\
*#define CHECK\n\n" <>
    Hdr[file, "assembly of squared matrix element"] <>
    fincl <> "\
\tsubroutine " <> $SymbolPrefix <> "SquaredME(result, helicities, flags)\n\
\timplicit none\n\
\tdouble precision result(*)\n\
\tinteger*8 helicities\n\
\tinteger flags\n\n\
#include \"" <> prefix <> "vars.h\"\n\n\
\tdouble precision " <> $SymbolPrefix <> "sumup\n\
\texternal " <> $SymbolPrefix <> "sumup\n" <>
    VarDecl[hels, "integer"] <>
    ({"\n\tequivalence (Hel(", #2, "), ", #1, ")"}&)@@@ hels <>
    ({"\n", #4, #5}&)@@@ loops <>
    ({"\n\n\tdata ", ToString[Head[#]],
      " /", ToString[Times@@ #], "*bogus/"}&)/@ mats <>
    "\n\n"];

  (* b) definitions of the invariants *)
  WriteExpr[hh, invs];

  (* c) calculation of the abbreviations *)
  WriteString[hh, "\n\
#define RESET_IF if( btest(flags, 0) ) then\n\
#define RESET_ENDIF endif\n\n\
#define LOOP_IF if( btest(flags, 1) ) then\n\
#define LOOP_ENDIF endif\n\n\
#define HEL_IF(i) if( btest(helicities, " <>
    ToString[5 legs + 2] <> "-5*i+Hel(i)) ) then\n\
#define HEL_ENDIF(i) endif\n\n\
\tRESET_IF\n" <>
    Invoke@@ abbrmods[[1]] <> "\
\tRESET_ENDIF\n\n\
\tcall markcache\n\n" <>
    Invoke@@ abbrmods[[2]] <> "\n\
\tresult(1) = 0\n\
\tresult(2) = 0\n\n" <>
    ({"\tdo ", #1, " = -2, 2\n\tHEL_IF(", #2, ")\n"}&)@@@ hels <> "\n" <>
    Invoke@@ abbrmods[[3]] <>
    ({#6, "\n\t", #2, " = 0", #7, "\n"}&)@@@ loops <>
    "\n"];

  (* d) calculation of the form factors *)
  If[ ntree =!= 0,
    WriteString[hh,
      CallDecl[ToDoLoops[ffmods[[1]], Indices]] <>
      Cond[ TrueQ[treesq],
        "\n\tresult(1) = result(1) + " <> $SymbolPrefix <> "sumup(",
          ntree[[3]] <> ", " <> ntree[[3]] <> ")\n" ]] ];

  If[ nloop =!= 0,
    WriteString[hh, "\n\tLOOP_IF\n" <>
      CallDecl[ToDoLoops[ffmods[[2]], Indices]] <>
      Cond[ ntree =!= 0,
        "\n\tresult(2) = result(2) + 2*", $SymbolPrefix, "sumup(",
          nloop[[3]], ", ", ntree[[3]], ")" ] <>
      Cond[ ntree === 0 || TrueQ[loopsq],
        "\n\tresult(2) = result(2) + ", $SymbolPrefix, "sumup(",
          nloop[[3]], ", ", nloop[[3]], ")" ] <>
      "\n\tLOOP_ENDIF\n"] ];

  (* e) wrapping up *)
  WriteString[hh,
    ({"\n\tHEL_ENDIF(", #2, ")\n\tenddo"}&)@@@ Reverse[hels] <>
    "\n\n\
\tcall restorecache\n\n\
#ifdef CHECK" <>
    ({"\n\tprint *, '", #, " =', ", #}&)/@
      (ToString[#1]&)@@@ invs <> "\n\n\
\tprint *, 'tree =', result(1)\n\
\tprint *, 'loop =', result(2)\n\
\tstop\n\
#endif\n\
\tend\n\n"];

  If[ ntree === 0, ntree = LoopComponents[Ctree, maxmat[Cloop]] ];
  If[ nloop === 0, nloop = LoopComponents[Cloop, maxmat[Ctree]] ];

  (* f) the sumup function *)
  WriteString[hh, "\n\
************************************************************************\n\n\
\tdouble precision function " <> $SymbolPrefix <> "sumup(C" <>
      nloop[[3]] <> ", C" <> ntree[[3]] <> ")\n\
\timplicit none\n\n\
#include \"" <> prefix <> "vars.h\"\n" <>
    nloop[[4]] <>
    ntree[[4]] <> "\n\
\tdouble complex C" <>
    StringReplace[nloop[[2]] <> ", C" <> ntree[[2]], "j" -> "n"] <> "\n\
\tdouble complex m\n\n\
\t" <> $SymbolPrefix <> "sumup = 0\n" <>
    ntree[[6]] <> "\n\
\tm = 0" <>
    nloop[[6]] <> "\n\
\tm = m + C" <> nloop[[2]] <> "*" <>
      ToFortran[Inner[MatType, nloop[[1]], ntree[[1]], Times]] <>
    nloop[[7]] <> "\n\
\t" <> $SymbolPrefix <> "sumup = " <>
      $SymbolPrefix <> "sumup + DBLE(DCONJG(C" <> ntree[[2]] <> ")*m)" <>
    ntree[[7]] <> "\n\
\tend\n"];

  Close[hh];

  Cases[DownValues[ModName], _[_, s_String] -> s]
]


(* renormalization constants *)

RCPattern[other___] := Alt[{other,
  Union[kindPattern[ #[[1, 1, 1]] ]&/@ DownValues[RenConst]]}]


Attributes[WithOpt] = {HoldFirst}

WithOpt[foo_, {}] := (Needs["FeynArts`"]; foo)

WithOpt[foo_, {opt__}] := (
  Needs["FeynArts`"];
  (Options[InsertFields] = #1; #2)&[
    Options[InsertFields],
    SetOptions[InsertFields, opt]; foo ] )

WithOpt[foo_, opt__] := WithOpt[foo, Flatten[{opt}]]


(* These are special versions of Re and Im where the real and
   imaginary part is taken only of the loop integrals,
   see A. Denner, Forts. Phys. 41 (1993) 307, arXiv:0709.1075. *)

ReTilde[expr_] := expr /. int:allint[__] :> Re[int]

ImTilde[expr_] :=
  (expr /. int:allint[__] :> Im[int]) - (expr /. allint[__] -> 0)


	(* Note: it seems weird that the left-handed vector component
	   is taken as the coefficient of DiracChain[6, k]: this is
	   because DiracChain[6, k] = DiracChain[k, 7]. *)

DiracCoeff[expr_, g__] :=
Block[ {e, ec},
  _e = _ec = Sequence[];
  (((# /. DiracChain[g] -> 1) - #) /. _DiracChain -> 0)& @ expr
]

LVectorCoeff[se_] := DiracCoeff[se, 6, k[1]]

RVectorCoeff[se_] := DiracCoeff[se, 7, k[1]]

LScalarCoeff[se_] := DiracCoeff[se, 7]

RScalarCoeff[se_] := DiracCoeff[se, 6]


SEPart[f_, se_] := f[se]


Attributes[OnlyIf] = {HoldRest}

OnlyIf[True, a_, _] = a

OnlyIf[False, _, b_] = b

OnlyIf[other__] := OnlyIfEval[other]

OnlyIfEval[_, a_, a_] = a

OnlyIfEval[cond_, a_, b_] := Thread @ IndexIf[ cond,
  a /. Cases[{cond}, i_ == j_ -> (IndexDelta[i, j] -> 1), Infinity],
  b /. Cases[{cond}, i_ == j_ -> (IndexDelta[i, j] -> 0)] ]


FermionRC[m_, a_, b_, se_] :=
  m (a SEPart[LVectorCoeff, se] + b SEPart[RVectorCoeff, se]) +
     b SEPart[LScalarCoeff, se] + a SEPart[RScalarCoeff, se]

BosonRC[se_] := SEPart[Identity, se]


MassRC[f_F, opt___Rule] :=
  FermionRC[TheMass[f], 1/2, 1/2, ReTilde[SelfEnergy[f, opt]]]

MassRC[f_, opt___Rule] := BosonRC[ReTilde[SelfEnergy[f, opt]]]

MassRC[f_, f_, opt___Rule] := MassRC[f, opt]

MassRC[f1_, f2_, opt___Rule] :=
  1/2 (BosonRC[ReTilde[SelfEnergy[f2 -> f1, TheMass[f1], opt]]] +
       BosonRC[ReTilde[SelfEnergy[f2 -> f1, TheMass[f2], opt]]])


FieldRC[f_F, opt___Rule] :=
Block[ {m = TheMass[f], se},
  se = ReTilde[SelfEnergy[f, opt]];
  -{SEPart[LVectorCoeff, se], SEPart[RVectorCoeff, se]} -
    FermionRC[m, m, m, ReTilde[DSelfEnergy[f]]]
]

FieldRC[f_, opt___Rule] := -BosonRC[ReTilde[DSelfEnergy[f, opt]]]

FieldRC[f_, f_, opt___Rule] := FieldRC[f, opt]

FieldRC[f1_, f2_, opt___Rule] := FieldRC[f1, f2, 0, opt]

FieldRC[f1_, f2_, c_, opt___Rule] := OnlyIf[
  And@@ MapThread[Equal, Flatten[{##}]&@@@ {f1, f2}],
  FieldRC[f1, opt],
  FieldRC2[f1, f2, c, opt]
]

FieldRC2[f1_F, f2_F, c_, opt___Rule] :=
Block[ {m1 = TheMass[f1], m2 = TheMass[f2]},
  2/(m1^2 - m2^2) (FermionRC[m2, {m2, m1}, {m1, m2},
    ReTilde[SelfEnergy[f2 -> f1, m2, opt]]] - c)
]

FieldRC2[f1_, f2_, c_, opt___Rule] :=
Block[ {m1 = TheMass[f1], m2 = TheMass[f2]},
  2/(m1^2 - m2^2) (BosonRC[ReTilde[SelfEnergy[f2 -> f1, m2, opt]]] - c)
]


TadpoleRC[f_, opt___Rule] :=
  -BosonRC[ReTilde[SelfEnergy[f -> {}, Indeterminate, opt]]]


WidthRC[f_F, opt___Rule] :=
  FermionRC[TheMass[f], 1, 1, ImTilde[SelfEnergy[f, opt]]]

WidthRC[f_, opt___Rule] :=
  BosonRC[ImTilde[SelfEnergy[f, opt]]]/TheMass[f]


Attributes[SelfEnergy] = {HoldRest}

SelfEnergy[proc_Rule, m_, opt___Rule] := (# /. K2 -> m^2)& @
  WithOpt[CalcSelfEnergy[proc, Options[InsertFields]], {opt}]

SelfEnergy[f_, opt___Rule] := SelfEnergy[f -> f, TheMass[f], opt]

SelfEnergy[f_, m__] := SelfEnergy[f -> f, m]


Attributes[DSelfEnergy] = {HoldRest}

DSelfEnergy[proc_Rule, m_, opt___Rule] := (# /. K2 -> m^2)& @
  WithOpt[D[CalcSelfEnergy[proc, Options[InsertFields]], K2], {opt}]

DSelfEnergy[f_, opt___Rule] := DSelfEnergy[f -> f, TheMass[f], opt]

DSelfEnergy[f_, m__] := DSelfEnergy[f -> f, m]


CalcSelfEnergy[proc_, opt_] := CalcSelfEnergy[proc, opt] =
Block[ {se, Neglect},
  ClearProcess[];
  se = InsertFieldsHook[
    CreateTopologies[1, Length[Flatten[{#}]]&/@ proc,
      ExcludeTopologies -> Internal],
    proc ];
  OptPaint[se];
  se = CreateFeynAmpHook[se, Truncated -> !FreeQ[proc, F]];
  Plus@@ CalcFeynAmp[se, OnShell -> False, Transverse -> False,
           FermionChains -> Chiral, AbbrScale -> 1, CutTools -> False] //.
    Abbr[] /. {
    Mat -> Identity,
    Pair[_k, _k] -> K2,
    Pair[_e | _ec, _k] -> If[ MatchQ[proc, _V -> _V],
	(* default: take only the transverse part of vector-boson SEs *)
      If[TrueQ[$LongitudinalSE], I Sqrt[K2], 0],
      1 ],
    Pair[_e, _ec] -> -1,
    SUNT[_, _] -> 1,
    SUNT[_, _, 0, 0] -> 1/2 }
]

InsertFieldsHook[args__] := InsertFields[args]

CreateFeynAmpHook[args__] := CreateFeynAmp[args]


ClearSE[] := (DownValues[CalcSelfEnergy] =
  Select[DownValues[CalcSelfEnergy], #[[1, 1, 1, 0]] === Pattern &];)


OptPaint[ins_] := OptPaint[ins, $PaintSE]

OptPaint[ins_, True, ___] := Paint[ins]

OptPaint[ins_, prefix_String, suffix___String] :=
Block[ {file = ChkExist[prefix, ProcessName[ins] <> suffix <> ".ps"]},
  Paint[ins, DisplayFunction -> (Export[file, #]&)]
]


RenConst::nodef =
"Warning: `` might be renormalization constants, but have no definition."

FindRenConst[expr_] :=
Block[ {test = {expr}, orbit, isym, patt, rcs = {}, new,
SelfEnergy, DSelfEnergy},
  Needs["FeynArts`"];
  If[ $Model === "",
    InitializeModel[
      Model /. Cases[test, HoldPattern[Model -> _], Infinity] /.
        Options[InsertFields],
      GenericModel -> (GenericModel /. Options[InsertFields]),
      Reinitialize -> False ] ];

  Apply[ (orbit[#1] = Range[#2])&,
    { Cases[test, SumOver[i_, r_, ___] :> {i, r}, Infinity],
      Cases[test, IndexSum[_, r___] :> r, Infinity] (*,
      Cases[DownValues[Dim], _[_[_[i_]], r_Integer] :> {i, r}]*) }, {2} ];
  orbit[other_] = other;

  Cases[DownValues[Dim], _[_[_[i_]], r_Integer] :> (isym[i] = x)];
  isym[other_] = other;

  patt = RCPattern[];

  While[ Length[new = Complement[
      Flatten[Cases[test, rc:patt :> Distribute[orbit/@ rc, List], Infinity]], 
      rcs, SameTest -> (isym/@ #1 === isym/@ #2 &) ]] =!= 0,
    rcs = Flatten[{new, rcs}];
    test = RenConst/@ new ];

  new = Select[ Names["Global`d*"],
    (FreeQ[rcs, #] && !FreeQ[expr, #]&)[
      ToExpression[#, InputForm, HoldPattern] ]& ];
  If[ Length[new] =!= 0, Message[RenConst::nodef, new] ];

  Cases[rcs, patt]
]


CalcRenConst[expr_, opt___Rule] :=
  (# -> WithOpt[RenConst[#], Options[kind[#]], opt])&/@
    FindRenConst[expr] /. Plus -> IntCollect


IntCollect[p__] := Plus[p] /; FreeQ[{p}, Re]

(*IntCollect[p__] := Collect[Plus[p], _Re, Simplify]*)

IntCollect[p__] :=
Block[ {Simplify},
  Replace[
    Collect[ Plus[p],
      First/@ DeleteCases[Split @ Sort @ Cases[Plus[p], _Re, Infinity], {_}],
      Simplify ],
    Simplify[x_] -> x, {1} ]
]


RCMod[rcs_, mod_] :=
Block[ {file = prefix <> mod <> ".F", hh},
  hh = OpenFortran[ModName[file]];
  WriteString[hh,
    Hdr[file, "renormalization constants"] <>
    fincl <>
    SubroutineDecl[mod] <>
    sincl <>
    VarDecl[Union[Cases[rcs, SumOver[i_, _] -> i, Infinity]], "integer"] <>
    "\n\n"];
  WriteExpr[hh, rcs, Optimize -> True, DebugLines -> True];
  WriteString[hh, "\tend\n"];
  Close[hh];
  mod
]

RCAll[mod_, mods_] :=
Block[ {file = prefix <> mod <> ".F", hh},
  hh = OpenFortran[ModName[file]];
  WriteString[hh,
    Hdr[file, "RC invocations"] <>
    fincl <>
    SubroutineDecl[mod] <>
    ({"\n\tcall ", $SymbolPrefix, #}&/@ mods) <>
    "\n\tend\n"];
  Close[hh];
  {mod, mods}
]


Options[WriteRenConst] = {
  Folder -> "renconst",
  ExtraRules -> ExtraRules (* i.e. from WriteSquaredME *),
  FilePrefix -> FilePrefix,
  SymbolPrefix -> SymbolPrefix,
  FileHeader -> FileHeader,
  FileIncludes -> FileIncludes,
  SubroutineIncludes -> SubroutineIncludes }

WriteRenConst::norcs = "Warning: no renormalization constants found."

WriteRenConst[rcs:{___Rule}, dir_, opt___Rule] :=
Block[ {folder, xrules, prefix, $SymbolPrefix, fincl, sincl, header,
ModName, Hdr, rcmods, file, hh},

  {folder, xrules, prefix, $SymbolPrefix, header, fincl, sincl} =
    ParseOpt[WriteRenConst, opt] /. Options[WriteSquaredME];

  (ModName[mod_] := ModName[mod] = ToFileName[#, mod])& @
    MkDir[dir, folder];

  header = StringReplace[header, "%t" -> TimeStamp[]];
  Hdr[f_, d_] := StringReplace[header, {"%f" -> f, "%d" -> d}];

(* Part 1: CalcRenConst.F *)

  If[ Length[rcs] === 0, Message[WriteRenConst::norcs] ];

  rcmods = FileSplit[ToDoLoops[OnePassOrder[rcs /. xrules]],
    "CalcRenConst", RCMod, RCAll];

(* Part 2: renconst.h *)

  file = prefix <> "renconst.h";
  hh = OpenFortran[ModName[file]];

  WriteString[hh,
    Hdr[file, "RC declarations"] <>
    VarDecl[MaxDims[Map[Dim, First/@ rcs, {2}]],
      {"double complex", "renconst"}] <> "\n"];

  Close[hh];

(* Part 3: the makefile *)

  hh = OpenWrite[ModName[prefix <> "makefile"]];

  WriteString[hh, "\
OBJS :=" <> ({" \\\n  $(DIR)/", #, ".o"}&)/@ Flatten[rcmods] <> "\n\n\
$(LIB): $(LIB)($(OBJS))\n\n\
$(LIB)($(OBJS)): $(DECL_H)\n\n\
LIBS += $(LIB)\n\n"];

  Close[hh];

  Cases[DownValues[ModName], _[_, s_String] -> s]
]

WriteRenConst[expr_, dir_, opt___Rule] :=
  WriteRenConst[CalcRenConst[expr], dir, opt]


(* low-level Fortran output functions *)

OpenFortran[file_, opt___] := OpenWrite[
  tofortran <> "\"" <> file <> "\"",
  FormatType -> FortranForm, opt, PageWidth -> 67 ]

tofortran =
  "!\"" <> ToFileName[{$FormCalcDir, $SystemID}, "ToFortran"] <> "\" > "


TimeStamp[] :=
  ToString[#3] <> " " <>
  {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
   "Sep", "Oct", "Nov", "Dec"}[[#2]] <> " " <>
  ToString[#1] <> " " <>
  ToString[#4] <> ":" <>
  StringTake["0" <> ToString[#5], -2]&@@ Date[]


(* The following routines are concerned with breaking a large
   expression into pieces the Fortran compiler will compile.
   This is controlled by two variables:

   - $BlockSize is the maximum LeafCount a single Fortran statement
     may have.  The cutting-up of expressions into such blocks is
     performed by the function WriteExpr.

   - $FileSize is the maximum LeafCount a whole file may have.
     The function which separates large expressions into file-size
     fragments is SizeSplit. *)

Attributes[batch] = {Flat}

Coalesce[(ru:Rule | RuleAdd)[v_, p_Plus], r___] :=
  Level[
    { ReplacePart[
        RuleAdd[v, Plus@@ #]&/@ Flatten[Coalesce@@ p],
        ru, {1, 0} ],
      {r} }, {2}, Coalesce ] /; LeafCount[p] > size

Coalesce[a_, b_, r___] :=
  Coalesce[batch[a, b], r] /; LeafCount[{a, b}] < size

Coalesce[a_, r___] := {batch[a], Coalesce[r]}

Coalesce[] = {}


Attributes[BlockSplit] = {Listable}

BlockSplit[expr_] := expr /; LeafCount[expr] < $BlockSize

BlockSplit[expr_] :=
Block[ {size = $BlockSize},
  List@@@ Flatten[Coalesce[expr]]
]


Attributes[NumberMod] = {Listable}

NumberMod[s__String] := StringJoin[s]


FileSplit[expr_List, mod_, writemod_, ___] :=
  {writemod[expr, mod]} /; LeafCount[expr] < $FileSize

FileSplit[expr_List, mod_, writemod_, writeall_:(#2&)] :=
Block[ {size = $FileSize},
  writeall[mod, MapIndexed[
    writemod[List@@ #1, NumberMod[mod, ToString@@ #2]]&,
    Flatten[Coalesce@@ expr] ]]
]

FileSplit[FortranExpr[vars__, expr_List], mod_,
  writemod_, writeall_:(#2&)] :=
Block[ {size = $FileSize},
  writeall[mod, MapIndexed[
    writemod[FortranExpr[vars, List@@ #1], NumberMod[mod, ToString@@ #2]]&,
    Flatten[Coalesce@@ expr] ]]
]

FileSplit[other_, r__] := FileSplit[{other}, r]


RemoveDups[expr_] :=
  Fold[RemoveDups, #, -Range[3, Depth[#] - 1]]& @ Flatten[{expr}]

RemoveDups[expr_, lev_] :=
Block[ {obj, elem, tmps, new},
  _obj = {};
  elem[{p_, r___}] := (obj[#] = {obj[#], p})& @ expr[[p, r]];
  elem/@ Position[ expr /. {_DoLoop -> 1, _IndexIf -> 1},
    p_Plus /; LeafCount[N[p]] > minleaf,
    {lev}, Heads -> False ];
  tmps = Cases[ DownValues[obj],
    _[_[_[x_]], p_ /; Depth[p] > 3] :> {x -> dup[++dupc], Min[p]} ];
  new = Block[{Plus}, Apply[Set, tmps, {2}]; expr];
  new = #1 /. Flatten[#2] & @@
    Reap[new /. r:(_dup -> _dup) :> (Sow[r]; {})];
  Fold[ InsertDef, new, Reverse[tmps] ] //Flatten
]

InsertDef[expr_, {var_ -> tmp_, pos_}] :=
  MapAt[{tmp -> var, #}&, expr, pos]


Attributes[TmpList] = {HoldFirst}

TmpList[expr_] := Reverse[Reap[expr]]

ToTmp[expr_] := (Sow[# -> expr]; #)& @ Unique["tmp"]


Attributes[SplitExpr] = {Listable}

SplitExpr[(ru:Rule | RuleAdd)[var_, expr_]] :=
  BlockSplit @ TmpList[ ru[var, Replace[expr,
    Plus -> (If[LeafCount[#] > $BlockSize, ToTmp[#], #]&[Plus[##]]&),
    {2, Infinity}, Heads -> True]] ]

SplitExpr[other_] = other


Attributes[RhsApply] = {Listable}

RhsApply[foo_, (ru:Rule | RuleAdd)[var_, expr_]] := ru[var, foo[expr]]

RhsApply[_, other_] = other


ListReduce[{a_}] = a

ListReduce[other_] = other


Options[PrepareExpr] = {
  Optimize -> False,
  Expensive -> {},
  MinLeafCount -> 10,
  DebugLines -> False,
  FinalTouch -> Identity }

PrepareExpr[expr_, opt___Rule] :=
Block[ {optim, expen, minleaf, debug, final,
process, doloop, new, vars, tmps, dup, dupc = 0},
  {optim, expen, minleaf, debug, final} = ParseOpt[PrepareExpr, opt];
  process = RhsApply[final, ListReduce @ Flatten[SplitExpr @ Prep[#]]] &;
  If[ TrueQ[optim], process = process /. p_Prep :> RemoveDups[p] ];
  doloop = Hoist@@ Flatten[{expen}];
  new = Flatten[{expr}];
  vars = Cases[new, (ru:Rule | RuleAdd)[var_, _] -> var, Infinity];
  If[ debug =!= False, new = AddDebug[new, debug] ];
  new = process[new];
  dup[c_] := dup[c] = Unique["dup"];
  new = new;
  tmps = Cases[new, (ru:Rule | RuleAdd)[var_, _] -> var, Infinity];
  FortranExpr[MaxDims[vars], Complement[tmps, vars], new]
]


Attributes[AddDebug] = {Listable}

AddDebug[s_String, _] = s

AddDebug[DoLoop[expr_, ind__], tag_] :=
  DoLoop[{DebugLine[#1&@@@ {ind}, tag], AddDebug[expr, tag]}, ind]

AddDebug[(ru:Rule | RuleAdd)[var_, expr_], tag_] :=
  {ru[var, expr], DebugLine[var, tag]}

DebugLine[var_, True] := DebugLine[var]


Attributes[Prep] = {Listable}

Prep[DoLoop[expr_, ind__]] := doloop[process[expr], ind]

Prep[(ru:Rule | RuleAdd)[var_, i_IndexIf]] :=
  MapIf[process[ru[var, #]]&, i]

Prep[(ru:Rule | RuleAdd)[var_, expr_List]] := Prep @
  MapIndexed[IniLHS[ru, var], ToDoLoops[expr] /. _SumOver -> 1]

Prep[(ru:Rule | RuleAdd)[var_, expr_]] := Prep @
  ru[var, SplitSums[expr]] /; !FreeQ[expr, SumOver]

Prep[(ru:Rule | RuleAdd)[var_, expr_]] := Prep @
  TmpList[ru[var, expr /. i_IndexIf :> ToTmp[i]]] /;
  !FreeQ[expr, IndexIf]

Prep[other_] = other


IniLHS[Rule, lhs_][DoLoop[rhs_, ind__], {1}] :=
  {lhs -> 0, DoLoop[RuleAdd[lhs, Plus@@ rhs], ind]}

IniLHS[Rule, lhs_][rhs_, {1}] := lhs -> rhs

IniLHS[_, lhs_][DoLoop[rhs_, ind__], _] :=
  DoLoop[RuleAdd[lhs, Plus@@ rhs], ind]

IniLHS[_, lhs_][rhs_, _] := RuleAdd[lhs, rhs]


Hoist[] = DoLoop

Hoist[a_, b__] := Hoist[a | b]

Hoist[patt_][expr_, i__] :=
Block[ {veto, abb, got = {}, c},
  veto = Alt@@ Union[
    DoIndex/@ Cases[expr, DoLoop[_, j__] :> j, Infinity],
    Cases[expr, SumOver[j_, ___] -> j, Infinity] ];
  abb = Cases[expr, x:patt /; FreeQ[x, veto], Infinity] //Union;
  abb = (got = Join[got, c = Complement[#2, got]]; #1 -> c)&@@@
    Sort[hsel[abb]/@ {i}, Length[ #1[[2]] ] < Length[ #2[[2]] ]&];
  abb[[1, 2]] = {};
  abb = (#1 -> (Unique["tmp"] -> # &)/@ #2)&@@@ abb;
  Fold[hdo, expr /. Reverse/@ Flatten[Last/@ abb], abb]
]


hsel[abb_][it:{s_, ___}] := it -> Select[abb, !FreeQ[#, s]&]

hsel[abb_][it_] := it -> Select[abb, !FreeQ[#, it]&]


hdo[DoLoop[expr_, j__], i_ -> {}] := DoLoop[expr, i, j]

hdo[expr_, i_ -> {}] := DoLoop[expr, i]

hdo[expr_, i_ -> {t__}] := DoLoop[{t, expr}, i]


SplitSums[li_List, wrap___] := SplitSums[Plus@@ li, wrap]

SplitSums[x_, wrap_:Identity] := {wrap[x]} /; FreeQ[x, SumOver]

SplitSums[x_, wrap_:Identity] :=
Block[ {term},
  term[_] = 0;
  assign[Expand[x, SumOver]];
  term[_] =.;
  #[[1, 1, 1]] wrap[Plus@@ Flatten[ #[[2]] ]]&/@ DownValues[term]
]

assign[p_Plus] := assign/@ p

assign[t_Times] := (term[#1] = {term[#1], #2})&@@ cull/@ t

assign[other_] := term[1] = {term[1], other}

cull[o_SumOver] := {o, 1}

cull[other_] := {1, other}


FindIndices[var_ -> _] := Union[Cases[var, _Symbol]]

FindIndices[t_Times] := Cases[t, SumOver[i__] -> {i}]

_FindIndices = {}

ToDoLoops[li:{__}, indices_:FindIndices] :=
Block[ {do, si},
  _do = {};
  Scan[(do[#1] = {do[#1], Tag[##]})&[indices[#], #]&, li];
  si = Flatten[Cases[DownValues[do], _[_[_[{___}]], a_] -> a]];
  DoLoop[ Last/@ #, Sequence@@ #[[1, 1]] ]&/@ 
    Split[OnePassOrder[si], #1[[1]] === #2[[1]] &]
]

(* old version:
ToDoLoops[li_List, indices_:FindIndices] :=
Block[ {do},
  do[_] = {};
  Scan[(do[#1] = {do[#1], #2})&[indices[#], #]&, li];
  Cases[ DownValues[do],
    _[_[_[{ind___}]], a_] :> DoLoop[Flatten[a], ind] ]
]
*)

ToDoLoops[other_, ___] = other


DoLoop[{a___}] = a

DoLoop[a_] = a


DoIndex[{i_, ___}] = i

DoIndex[i_] = i


Options[WriteExpr] = {
  HornerStyle -> True,
  FinalCollect -> False,
  Type -> False,
  TmpType -> (*Type*) "double complex",
  IndexType -> False,
  RealArgs -> Level[{allint, Bget, Cget, Dget, Eget, Fget, Log, Sqrt}, {-1}],
  Newline -> "\n" }

WriteExpr[_, _[], ___] = {}

WriteExpr[hh_, FortranExpr[vars_, tmpvars_, expr_], opt___Rule] :=
Block[ {horner, fcoll, type, tmptype, indtype, rargs, newline,
var, block = 0},
  {horner, fcoll, type, tmptype, indtype, rargs, newline} =
    ParseOpt[WriteExpr, opt];
  rargs = Alt[rargs];
  horner = If[ horner =!= True, {},
    p_Plus :> (HornerForm[p] /. HornerForm -> Identity) /; Depth[p] < 6 ];
  fcoll = If[ fcoll =!= True, {},
    a_ b_. + a_ c_ -> a (b + c) ];
  _var = {};
  VarType[vars, type];
  VarType[tmpvars, tmptype /. Type -> type];
  VarType[Union[DoIndex/@
    Cases[expr, DoLoop[_, i__] -> i, Infinity]], indtype];
  _var =.;
  WriteString[hh, StringDrop[
    VarDecl[ Flatten[#2], #1[[1, 1]] ]&@@@ DownValues[var] <>
    "\n\n", 1 ]];
  Flatten[{WriteBlock[hh, expr]}]
]

WriteExpr[hh_, expr_, opt___Rule] := WriteExpr[hh,
  PrepareExpr[expr, FilterOptions[PrepareExpr, opt]],
  FilterOptions[WriteExpr, opt]]


VarType[_, False] = 0

VarType[v:{__}, s:_String | {__String}] := var[s] = {var[s], v}

VarType[v_List, r:_Rule | {__Rule}] :=
  MapThread[VarSet, {v, Replace[v, r, {1}]}]

VarType[v_List, f_] := MapThread[VarSet, {v, f/@ v}]

VarSet[v_, v_] = 0

VarSet[v_, s:_String | {__String}] := var[s] = {var[s], v}


$DebugCmd = {"DEBUG", "DEB", ""}

DebugStatement[var_, tag_, {"", cmd_, end_}] :=
  "\t" <> cmd <> " '" <> ({ToString[#], ": "}&)/@ tag <>
    var <> " =', " <> var <> end <> "\n"

DebugStatement[v__, {pre_, c__}] :=
  "#ifdef " <> pre <> "\n" <>
    DebugStatement[v, {"", c}] <>
  "#endif\n"


Attributes[WriteBlock] = {Listable}

WriteBlock[hh_, s_String] := (WriteString[hh, s <> newline]; s)

WriteBlock[hh_, DebugLine[var_, tag___]] :=
  WriteBlock[hh, DebugStatement[ToFortran[var], {tag}, $DebugCmd]]

WriteBlock[hh_, DoLoop[expr_, ind__]] :=
  WriteBlock[hh, {#1, expr, #2}]&@@ DoDecl[ind]

WriteBlock[hh_, IndexIf[cond_, a_, r___]] :=
  WriteBlock[hh, {
    "\tif( " <> ToFortran[cond] <> " ) then\n", a,
    ElseIf[r], "\tendif\n" }]

ElseIf[] = ElseIf[RuleAdd[_, 0]] = {}

ElseIf[a_] = {"\telse\n", a}

ElseIf[cond_, a_, r___] :=
  {"\telse if( " <> ToFortran[cond] <> " ) then\n", a, ElseIf[r]}

WriteBlock[hh_, ru_[var_, {sub__, expr_}]] :=
  WriteBlock[hh, {sub, ru[var, expr]}]

WriteBlock[_, RuleAdd[_, 0]] := Sequence[]

WriteBlock[hh_, RuleAdd[var_, expr_]] := (
  Write[hh, var -> var + FExpr[expr]];
  WriteString[hh, newline];
  var -> var + expr
)

WriteBlock[hh_, var_ -> expr_] := (
  Write[hh, var -> FExpr[expr]];
  WriteString[hh, newline];
  var -> expr
)


FExpr[expr_] := expr /.
  {Conjugate -> DCONJG, Re -> DBLE, Im -> DIMAG,
    Complex[a_, b_] -> a + cI b,
    Dminus4 -> -2/Divergence} /.
  E^x_ :> exp[x] /.
  f:rargs[__] :> RealArgs[f] /.
  p_Integer^r_Rational :> (HoldForm[#]^r &)[ N[p] ] /.
  (* p:_Integer^_Rational :> N[p] /. *)
  Den[p_, m_, d___] :> (p - m)^-d /.
  horner //. fcoll /.
  Times -> OptTimes


Scan[(RealArgs[#[i_, r__]] := #[i, NArgs[r]])&, cutint]

RealArgs[f_] := NArgs/@ f


NArgs[i_Integer] := N[i]

NArgs[other_] := other

NArgs[i__] := Sequence@@ NArgs/@ {i}


Attributes[NZap] = {HoldAll}

NZap[sym_] := (NValues[sym] = {}; sym) /;
  Length[NValues[sym]] > 0

_NZap = {}

NClear[patt_String:"Global`*"] := Flatten[
  ToExpression[#, InputForm, NZap]&/@ Names[patt]
] //Union


Unprotect[Rule, Rational, Power]

Format[a_ -> b_, FortranForm] := SequenceForm[a, " = ", b]

Format[Rational[a_, b_], FortranForm] := HoldForm[a]/b

Format[a_^n_Integer /; n < -1, FortranForm] := (1/HoldForm[#] &)[ a^-n ]

Protect[Rule, Rational, Power]


OptTimes[t__] :=
Block[ {p = Position[N[{t}], _Real, 1, Heads -> False]},
  OptNum[Times@@ Extract[{t}, p], Times@@ Delete[{t}, p]]
]

OptNum[const_, 1] = const

OptNum[const_Integer, var_] := const var

OptNum[n_?Negative r_., var_] := -OptNum[-n r, var]

OptNum[const_, var_] := HoldForm[HoldForm[const] var]

End[]


Format[ _Continuation ] = "    "
  (* eliminate those `>' in front of continuation lines so one can cut
     and paste more easily *)

$FormCmd = ToFileName[{$FormCalcDir, "FORM"}, "form_" <> $SystemID]
  (* the filename of the actual FORM executable; may contain a path *)

FormSetup = "\
#-\n\
#:SmallSize 5000000\n\
#:LargeSize 20000000\n\
#:WorkSpace 20000000\n\
#:MaxTermSize 300000\n\
#:TermsInSmall 30000\n\
off stats;\n\n"

$Editor = "${VISUAL:-xterm -e nano} `` &"
  (* editor to use when debugging FORM code *)

$BlockSize = 700

$FileSize = 30 $BlockSize

$MaxFortranName = 30

$RecursionLimit = 1024

$DriversDir = ToFileName[{$FormCalcDir, "drivers"}]

$PaintSE = False

EndPackage[]


(* global definitions for specific models *)

powQ[1] = powQ[-1] = False

powQ[other_] := IntegerQ[other]

Sq/: (Sq[v_] = v2_) := (
  v^(n_?EvenQ) ^:= v2^(n/2);
  (v^(n_?powQ) ^:= v2^((n - Sign[n])/2) #^Sign[n])&[ v /. Pattern -> (#1&) ];
  Square[v] = v2
)

Sq/: (Sq[v_] =.) := (
  v/: v^(n_?EvenQ) =.;
  v/: v^(n_?powQ) =.;
  Square[v] =.
)


(* definitions for the Standard Model *)

Sq[EL] = 4 Pi Alfa;
Sq[Alfa] = Alfa2;
Alfa2/: Alfa2/Alfa = Alfa

Sq[GS] = 4 Pi Alfas;
Sq[Alfas] = Alfas2;
Alfas2/: Alfas2/Alfas = Alfas

Sq[SW] = SW2;
Sq[CW] = CW2;
CW2/: CW2 + SW2 = 1

Sq[MZ] = MZ2;
Sq[MW] = MW2;
Sq[MH] = MH2;

Sq[ME] = ME2;  Sq[MM] = MM2;  Sq[ML] = ML2;
Sq[MU] = MU2;  Sq[MC] = MC2;  Sq[MT] = MT2;
Sq[MD] = MD2;  Sq[MS] = MS2;  Sq[MB] = MB2

MLE[a__] := Mf[2,a];
MQU[a__] := Mf[3,a];
MQD[a__] := Mf[4,a];
Sq[Mf[a__]] = Mf2[a]

Mf[2,1] = ME;  Mf[2,2] = MM;  Mf[2,3] = ML;
Mf[3,1] = MU;  Mf[3,2] = MC;  Mf[3,3] = MT;
Mf[4,1] = MD;  Mf[4,2] = MS;  Mf[4,3] = MB

Mf2[2,1] = ME2;  Mf2[2,2] = MM2;  Mf2[2,3] = ML2;
Mf2[3,1] = MU2;  Mf2[3,2] = MC2;  Mf2[3,3] = MT2;
Mf2[4,1] = MD2;  Mf2[4,2] = MS2;  Mf2[4,3] = MB2

Conjugate[CKM[a__]] ^:= CKMC[a];
Conjugate[CKMC[a__]] ^:= CKM[a]

SUNN = 3

(* these symbols represent real quantities, i.e. Conjugate[sym] = sym
   for any of these.  Thinking e.g. of complex masses this looks
   dangerous but then again it's easy to remove any such definition.
   The function that really needs this is SquaredME. *)

Scan[ (RealQ[#] = True)&,
  { EL, Alfa, Alfa2, GS, Alfas, Alfas2,
    SW, CW, SW2, CW2,
    MW, MW2, MZ, MZ2,
    MH, MH2, MG0, MG02, MGp, MGp2,
    ME, ME2, MM, MM2, ML, ML2,
    MU, MU2, MC, MC2, MT, MT2,
    MD, MD2, MS, MS2, MB, MB2, _Mf, _Mf2 } ]

(* Model parameters which are defined using the parameter statement in
   Fortran (i.e. as numeric constants; see model.h) are given some
   numeric value here.  Using this information, the OptTimes function can
   significantly optimize the generated Fortran code.  The idea is to put
   everything that is known as constant at compile time in one place,
   i.e. rearrange products such that they are of the form (const)*(vars),
   then the compiler will usually collect all of these constants into one
   number. *)

Scan[ (N[#] = Random[])&,
  { cI, Alfa, Alfa2, SW2, CW, CW2,
    MW, MW2, MZ, MZ2,
    ME, ME2, MM, MM2, ML, ML2,
    MU, MU2, MC, MC2, MT, MT2,
    MD, MD2, MS, MS2, MB, MB2 } ]

SMReduce[foo_][expr_, r___] := foo[expr /.
  SW2 -> 1 - CW2 /.
  {CW -> MW/MZ, CW2 -> MW2/MZ2}, r]

SMShorten[foo_][x__] := SMReduce[foo][x] /.
  MW2 - MZ2 -> -SW2 MZ2 /.
  {MZ/MW -> 1/CW, MZ2/MW2 -> 1/CW2,
   MW/MZ -> CW, MW2/MZ2 -> CW2}

SMSimplify = SMShorten[Simplify];
SMFullSimplify = SMShorten[FullSimplify]

(* definitions for the MSSM *)

SetOptions[ CalcFeynAmp,
  NoExpand -> {USf, USfC, UASf, UASfC,
    UCha, UChaC, VCha, VChaC, ZNeu, ZNeuC,
    UHiggs, UHiggsC, ZHiggs, ZHiggsC},
  NoBracket -> {CKM, CKMC, USf, USfC, UASf, UASfC,
    UCha, UChaC, VCha, VChaC, ZNeu, ZNeuC,
    UHiggs, UHiggsC, ZHiggs, ZHiggsC} ]

Af[t_, g_] := Af[t, g, g]

USf[t_, g_][a_, b_] := USf[a, b, t, g];
UASf[t_][a_, b_] := UASf[a, b, t]

Conjugate[USf[a__]] ^:= USfC[a];
Conjugate[USfC[a__]] ^:= USf[a]

Conjugate[UASf[a__]] ^:= UASfC[a];
Conjugate[UASfC[a__]] ^:= UASf[a]

Conjugate[UCha[a__]] ^:= UChaC[a];
Conjugate[UChaC[a__]] ^:= UCha[a]

Conjugate[VCha[a__]] ^:= VChaC[a];
Conjugate[VChaC[a__]] ^:= VCha[a]

Conjugate[ZNeu[a__]] ^:= ZNeuC[a];
Conjugate[ZNeuC[a__]] ^:= ZNeu[a]

Conjugate[UHiggs[a__]] ^:= UHiggsC[a];
Conjugate[UHiggsC[a__]] ^:= UHiggs[a]

Conjugate[ZHiggs[a__]] ^:= ZHiggsC[a];
Conjugate[ZHiggsC[a__]] ^:= ZHiggs[a]

Conjugate[Af[a__]] ^:= AfC[a];
Conjugate[AfC[a__]] ^:= Af[a]

Conjugate[MUE] ^:= MUEC;
Conjugate[MUEC] ^:= MUE

Conjugate[Mino3] ^:= Mino3C;
Conjugate[Mino3C] ^:= Mino3

Conjugate[SqrtEGl] ^:= SqrtEGlC;
Conjugate[SqrtEGlC] ^:= SqrtEGl

SqrtEGl/: SqrtEGl SqrtEGlC = 1

Sq[SqrtEGl] = Mino3/MGl;
Sq[SqrtEGlC] = Mino3C/MGl

Sq[SA] = SA2;
Sq[CA] = CA2;
CA2/: CA2 + SA2 = 1

Sq[TB] = TB2;
Sq[SB] = SB2;
Sq[CB] = CB2;
CB2/: CB2 + SB2 = 1;
CB2/: CB2 TB2 = SB2;
CB/: CB TB = SB

Sq[CBA] = CBA2;
Sq[SBA] = SBA2;
CBA2/: CBA2 + SBA2 = 1

SUSYTrigExpand[expr_] := expr /. {
  SB -> sb, CB -> cb, SB2 -> sb^2, CB2 -> cb^2,
  TB -> sb/cb, TB2 -> sb^2/cb^2,
  SA -> sa, CA -> ca, SA2 -> sa^2, CA2 -> ca^2,
  C2A -> ca^2 - sa^2, S2A -> 2 ca sa,
  C2B -> cb^2 - sb^2, S2B -> 2 cb sb,
  CAB -> ca cb - sa sb, SAB -> cb sa + ca sb,
  CBA -> ca cb + sa sb, SBA -> ca sb - cb sa,
  CBA2 -> (ca cb + sa sb)^2, SBA2 -> (ca sb - cb sa)^2 }

Sq[MGl] = MGl2;
Sq[MSf[a__]] = MSf2[a];
Sq[MASf[a__]] = MASf2[a];
Sq[MCha[a__]] = MCha2[a];
Sq[MNeu[a__]] = MNeu2[a]

Sq[MHiggs[a__]] = MHiggs2[a];
Sq[MHiggstree[a__]] = MHiggstree2[a]

Sq[Mh0] = Mh02;  Sq[Mh0tree] = Mh0tree2;
Sq[MHH] = MHH2;  Sq[MHHtree] = MHHtree2;
Sq[MA0] = MA02;  Sq[MA0tree] = MA0tree2;
Sq[MHp] = MHp2;  Sq[MHptree] = MHptree2

Scan[ (RealQ[#] = True)&,
  { TB, CB, SB, TB2, CB2, SB2, C2B, S2B,
    CA, SA, CA2, SA2, C2A, S2A,
    CAB, SAB, CBA, SBA, CBA2, SBA2,
    Mh0, Mh02, MHH, MHH2, MA0, MA02, MHp, MHp2, MGl,
    _MSf, _MSf2, _MCha, _MCha2, _MNeu, _MNeu2 } ]

MSSMReduce[foo_, red_:SMReduce][expr_, r___] :=
  red[foo][expr /. SBA2 -> 1 - CBA2, r]

MSSMShorten[foo_, red_:SMSimplify][x__] :=
  MSSMReduce[foo, red][x] /. CBA2 - 1 -> -SBA2

MSSMSimplify = MSSMShorten[Simplify];
MSSMFullSimplify = MSSMShorten[FullSimplify]

(* make Needs["FeynArts`"] work after loading FormCalc *)

If[ !NameQ["FeynArts`$FeynArts"],
  Unprotect[$Packages];
  $Packages = DeleteCases[$Packages, "FeynArts`"];
  Protect[$Packages] ]

Null

