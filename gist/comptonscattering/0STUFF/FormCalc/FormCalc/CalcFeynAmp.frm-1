* CalcFeynAmp.frm
* the FORM part of the CalcFeynAmp function
* this file is part of FormCalc
* last modified 14 May 10 th


***********************************************************************

#procedure Contract
repeat once e_([i]?, [j]?, [k]?, [LA]?)*e_([I]?, [J]?, [K]?, [LA]?) =
#if "`Dim'" != "4"
  (1 + Dminus4) *
#endif
  ( d_([i], [I]) * (d_([j], [J])*d_([k], [K]) - d_([j], [K])*d_([k], [J])) +
    d_([i], [J]) * (d_([j], [K])*d_([k], [I]) - d_([j], [I])*d_([k], [K])) +
    d_([i], [K]) * (d_([j], [I])*d_([k], [J]) - d_([j], [J])*d_([k], [I])) );
#endprocedure

***********************************************************************

#procedure DotSimplify(cmd)
#if "`MomElim'" == "Automatic"

.sort

#call eiki

ab `Vectors';
.sort

keep brackets;

#$term = 0;
$term = $term + 1;
mul TMP($term, termsinbracket_);
id TMP([x]?, 1) = 1;

b TMP, `Vectors';
.sort

keep brackets;

id TMP([x]?, ?r) = <TMP([x], 0)> + ... + <TMP([x], `Legs')>;

#do i = 1, `Legs'
#ifdef `k`i''
if( match(TMP([x]?, {`Legs'-`i'})) ) id k`i' = `k`i'';
#endif
#enddo

`cmd'
#call eiki
#call kikj

b TMP;
.sort

keep brackets;

id TMP([x]?, [y]?) = TMP([x], termsinbracket_, [y]);

b TMP;
.sort

keep brackets;

#$term = 0;
$prev = $term;
id TMP([x]?$term, ?a) = deltap_([x], $prev);

#else

b `Vectors';
.sort

keep brackets;

#ifdef `k`MomElim''
#call eiki
id k`MomElim' = `k`MomElim'';
#endif

`cmd'
#call eiki
#call kikj

#endif
#endprocedure

***********************************************************************

#procedure DiracOrder
label 1;

#if `OnShell' == 1

* Apply Dirac equation to right spinor

repeat;
  id GA([om]?, ?a, [p1]?, ?b) * Spinor([p1]?, [m1]?, [s1]?) =
    ( 2*GD([om], ?a, [p1]) * distrib_(-1, 1, GD, GD, ?b) +
      sign_(nargs_(?b)) * [s1]*[m1] * GA([om], ?a, ?b) ) *
    Spinor([p1], [m1], [s1]);
  id GD(?a, [p1]?) * GD([mu]?) * GD(?b) =
    d_([p1], [mu]) * GA(?a, ?b) * TAG;
endrepeat;

* Apply Dirac equation to left spinor

repeat;
  id Spinor([p1]?, [m1]?, [s1]?) * GA([om]?{6,7}[[n]], ?a, [p1]?, ?b) =
    Spinor([p1], [m1], [s1]) * sign_(nargs_(?a)) *
    ( [s1]*[m1] * GA({7,6}[[n]], ?a, ?b) -
      2*distrib_(-1, 1, GD, GD, ?a) * GD([p1], [om], ?b) );
  id GD([mu]?) * GD(?a) * GD([p1]?, [om]?, ?b) =
    d_([p1], [mu]) * GA([om], ?a, ?b) * TAG;
endrepeat;
#endif

* Eliminate contractions within each Dirac chain using the
* formulas from M. Veltman's Gammatrica [Nucl Phys B319 (1989) 253]

id GA([om]?, [mu]?, ?a) = GA([om]) * GB([mu], ?a);
while( count(GB, 1) );
  repeat;
    id GB([LA]?, [LA]?, ?a) = d_([LA], [LA]) * GB(?a);
    also GB([LA]?, [mu]?, [LA]?, ?a) = (2 - d_([LA], [LA])) * GB([mu], ?a);
    also GB([LA]?, [mu]?, [nu]?, [LA]?, ?a) =
#if "`Dim'" != "4"
      Dminus4 * GB([mu], [nu], ?a) +
#endif
      4*d_([mu], [nu]) * GB(?a) * TAG;
    also GB([LA]?, [mu]?, [nu]?, [ro]?, ?b, [LA]?, ?a) =
      -sign_(nargs_(?b)) * (
#if "`Dim'" != "4"
        Dminus4 * GB([mu], [nu], [ro], ?b, ?a) +
#endif
        2*GB([ro], [nu], [mu], ?b, ?a) +
        2*GD([mu], [nu], [ro]) * distrib_(-1, 1, GD, GD, ?b) * GD(?a) );
    id GD(?a) * GD([mu]?) * GD(?b) * GD(?c) = GB([mu], ?a, ?b, ?c);
  endrepeat;
  id GB([mu]?, ?a) = GC([mu]) * GB(?a);
  id GB() = 1;
endwhile;

* Order the gamma matrices canonically

repeat;
  id GC([p1]?) * GC([p1]?) = [p1].[p1];
  disorder GC([mu]?) * GC([nu]?) = 2*d_([mu], [nu]) * TAG - GC([nu]) * GC([mu]);
endrepeat;
repeat id GA(?a) * GC(?b) = GA(?a, ?b);

id ifmatch->1 TAG = 1;
#endprocedure

***********************************************************************

#procedure DiracSimplify
#call DiracOrder

#if `OnShell' == 1

#do i = 1, `Legs'
#ifdef `k`i''
b `Tensors';
.sort

keep brackets;

id k`i' = `k`i'';
#call DiracOrder
#endif
#enddo

#endif
#endprocedure

***********************************************************************

#procedure FierzPre(ord)
b `Tensors', ORD;
.sort

keep brackets;

id Spinor([p1]?, ?a) * GA(?g) * Spinor([p2]?, ?b) =
  CH(Spinor([p1], ?a), ?g, Spinor([p2], ?b))
#if `ord'
  * ORD([p1]) * ORD([p2]);

repeat id ORD([p1]?) * ORD(?a, [p1]?, ?b) = ORD(?a, [p1], ?b);

chainin ORD
#endif
  ;

* these relations are obtained by Fierzing twice

id CH(Spinor(?a), 6, [mu]?, [nu]?, Spinor(?b)) *
   CH(Spinor(?c), 7, [mu]?, [nu]?, Spinor(?d)) =
  4*CH(Spinor(?a), 6, Spinor(?b)) *
    CH(Spinor(?c), 7, Spinor(?d));

id CH(Spinor(?a), [om]?, [mu]?, [nu]?, Spinor(?b)) *
   CH(Spinor(?c), [om]?, [ro]?, [mu]?, [nu]?, Spinor(?d)) =
  4*CH(Spinor(?a), [om], Spinor(?b)) *
    CH(Spinor(?c), [om], [ro], Spinor(?d));
also CH(Spinor(?a), [omA]?, [mu]?, [nu]?, Spinor(?b)) *
   CH(Spinor(?c), [omB]?, [ro]?, [mu]?, [nu]?, Spinor(?d)) =
  4*CH(Spinor(?a), [omA], [ro], [mu], Spinor(?b)) *
    CH(Spinor(?c), [omB], [mu], Spinor(?d));

id CH(Spinor(?a), 6, [ro]?, [mu]?, [nu]?, Spinor(?b)) *
   CH(Spinor(?c), 7, [si]?, [mu]?, [nu]?, Spinor(?d)) =
  4*CH(Spinor(?a), 6, [ro], Spinor(?b)) *
    CH(Spinor(?c), 7, [si], Spinor(?d));

id CH(Spinor(?a), [om]?, [mu]?, [nu]?, [ro]?, Spinor(?b)) *
   CH(Spinor(?c), [om]?, [mu]?, [nu]?, [ro]?, Spinor(?d)) =
  16*CH(Spinor(?a), [om], [mu], Spinor(?b)) *
     CH(Spinor(?c), [om], [mu], Spinor(?d));
also CH(Spinor(?a), [omA]?, [mu]?, [nu]?, [ro]?, Spinor(?b)) *
   CH(Spinor(?c), [omB]?, [mu]?, [nu]?, [ro]?, Spinor(?d)) =
  4*CH(Spinor(?a), [omA], [mu], Spinor(?b)) *
    CH(Spinor(?c), [omB], [mu], Spinor(?d));
#endprocedure

***********************************************************************

#procedure FierzPost
id CH([x]?, ?g, [y]?) = [x] * GA(?g) * [y];

repeat;
  once GA([om]?, [mu]?, [nu]?, [ro]?, ?a) =
    sum_([c], 1, 5,
      DUAL([c], MUC, NUC) * CHI([om])*g_(1, [mu], [nu], [ro], ?a) *
      BASIS([c], MUC, NUC) );
  trace4, 1;
  sum MUC, NUC;
endrepeat;

#call Contract

* Chisholm's identity backwards to get rid of all Eps
repeat;
  once GA([om]?, ?a, [LA]?, ?b) * e_([mu]?, [nu]?, [ro]?, [LA]?) =
    sign_([om]) * sign_(nargs_(?a)) * (
      GA([om], ?a, [mu], [nu], [ro], ?b) -
      d_([mu], [nu]) * GA([om], ?a, [ro], ?b) +
      d_([mu], [ro]) * GA([om], ?a, [nu], ?b) -
      d_([nu], [ro]) * GA([om], ?a, [mu], ?b) );
endrepeat;

id D = 4;
#call DiracSimplify

.sort

id D = 4;
#endprocedure

***********************************************************************

#procedure FierzUnordered
#call FierzPre(0)

* The following general Fierz identity is from hep-ph/0412245.
repeat;
  once CH(Spinor(?a), [omA]?, ?A, Spinor(?d)) *
       CH(Spinor(?c), [omB]?, ?B, Spinor(?b)) =
    sum_([c], 1, 5, sum_([d], 1, 5,
      DUAL([d], MUD, NUD) * CHI([omA])*g_(1, ?A) *
      DUAL([c], MUC, NUC) * CHI([omB])*g_(1, ?B) *
      Spinor(?a) * BASIS([d], MUD, NUD) * Spinor(?b) *
      Spinor(?c) * BASIS([c], MUC, NUC) * Spinor(?d) ));
  trace4, 1;
  sum MUC, NUC, MUD, NUD;
endrepeat;

#call FierzPost
#endprocedure

***********************************************************************

#procedure FierzOrdered
while( count(ORD, 1) );
  once ORD([p1]?, [p2]?, ?r) = NOW([p1], [p2]) * ORD(?r);
  id ORD() = 1;

  repeat;
    id NOW([p1]?, [p2]?) *
         CH(Spinor([p1]?, ?m1), ?g, Spinor([p2]?, ?m2)) =
      CH(Spinor([p1], ?m1), ?g, Spinor([p2], ?m2));

* charge conjugation to get first spinor in front
* the rules for this are exactly as in HelicityME
    id NOW([p1]?, [p3]?) * CH(Spinor([p2]?, [m2]?, [s2]?), [x]?, ?g,
                              Spinor([p1]?, [m1]?, [s1]?)) =
      -NOW([p1], [p3]) * sign_(nargs_(?g)) *
      CH(Spinor([p1], [m1], -[s1]),
         (sign_(nargs_(?g))*(2*[x] - 13) + 13)/2, reverse_(?g),
         Spinor([p2], [m2], -[s2]));

* charge conjugation to get second spinor in back
    id NOW([p3]?, [p2]?) * CH(Spinor([p2]?, [m2]?, [s2]?), [x]?, ?g,
                              Spinor([p1]?, [m1]?, [s1]?)) =
      -NOW([p3], [p2]) * sign_(nargs_(?g)) *
      CH(Spinor([p1], [m1], -[s1]),
         (sign_(nargs_(?g))*(2*[x] - 13) + 13)/2, reverse_(?g),
         Spinor([p2], [m2], -[s2]));

* Fierz to get second spinor together with first
    once NOW([p1]?, [p2]?) *
      CH(Spinor([p1]?, ?a), [omA]?, ?A, Spinor([p4]?, ?d)) *
      CH(Spinor([p3]?, ?c), [omB]?, ?B, Spinor([p2]?, ?b)) =
      NOW([p1], [p2]) * sum_([c], 1, 5, sum_([d], 1, 5,
        DUAL([d], MUD, NUD) * CHI([omA])*g_(1, ?A) *
        DUAL([c], MUC, NUC) * CHI([omB])*g_(1, ?B) *
        Spinor([p1], ?a) * BASIS([d], MUD, NUD) * Spinor([p2], ?b) *
        Spinor([p3], ?c) * BASIS([c], MUC, NUC) * Spinor([p4], ?d) ));
    trace4, 1;
    sum MUC, NUC, MUD, NUD;

    id Spinor(?a) * GA(?g) * Spinor(?b) =
      CH(Spinor(?a), ?g, Spinor(?b));
  endrepeat;
endwhile;

#call FierzPost
#endprocedure

***********************************************************************

#procedure DiracFinal

#if `Antisymmetrize' == 1

* introduce antisymmetrized Dirac chains
repeat once GA([om]?, ?g) = GB([om]) *
  sum_(DUMMY, 0, nargs_(?g), 2, distrib_(-1, DUMMY, DD, GD, ?g));

id DD() = 1;
repeat;
  once DD(?a) = g_(1, ?a)/4;
  trace4, 1;
endrepeat;

antisymm GD;
id D = 4;

id GB([om]?) * GD(?g) = GA([om], ?g);

#endif

#if "`FermionChains'" == "VA"
id GA([om]?, ?a) = GA(1, ?a)/2 + sign_([om]) * GA(5, ?a)/2;
#endif

b `Tensors';
.sort

keep brackets;

id Spinor(?a) * GA(?g) * Spinor(?b) =
  ABB(1, DiracChain(Spinor(?a), ?g, Spinor(?b)), ?g);

id GA(?g) = ABB(1, DiracChain(?g), ?g);

#if `Antisymmetrize' == 1
id ABB(1, DiracChain(Spinor(?a), [x]?, [mu]?, [nu]?, ?r), ?g) =
  ABB(1, DiracChain(Spinor(?a), -[x], [mu], [nu], ?r), ?g);
also ABB(1, DiracChain([x]?int_, [mu]?, [nu]?, ?r), ?g) =
  ABB(1, DiracChain(-[x], [mu], [nu], ?r), ?g);
#endif

repeat id ABB(1, [x]?, ?g) * ABB(1, [y]?, ?h) =
  ABB(1, [x]*[y], ?g, ?h);

id ABB(1, DiracChain(?a, -6, [mu]?, [nu]?, ?b) *
          DiracChain(?c, -7, [mu]?, [nu]?, ?d), ?g) = 0;

#call Abbreviate
#endprocedure

***********************************************************************

#procedure Abbreviate
#call DotSimplify()

.sort

*----------------------------------------------------------------------

#if "`CutTools'" != "False"

if( count(cutM, 1) );

id q1.[p1]? = qfM(Pair(q1, [p1]));
id e_(q1, [p1]?, [p2]?, [p3]?) = qfM(Eps(q1, [p1], [p2], [p3]));
id abbM(fmeM(WeylChain(?a, q1, ?b))) = qfM(WeylChain(?a, q1, ?b));

repeat id qfM([x]?) * qfM([y]?) = qfM([x] * [y]);

endif;

.sort

#endif

*----------------------------------------------------------------------

id [p1]?.[p2]? = ABB(0, [p1].[p2], [p1], [p2]);

id e_([mu]?, [nu]?, [ro]?, [si]?) =
  ABB(0, e_([mu], [nu], [ro], [si]), [mu], [nu], [ro], [si]);

id d_([mu]?, [nu]?) = ABB(0, d_([mu], [nu]), [mu], [nu]);

id [t]?(?a) = ABB(0, [t](?a), ?a);

id [p1]?([mu]?) = ABB(0, [p1]([mu]), [p1]);

repeat;
  once ABB([s1]?, [x]?, ?a, [mu]?!fixed_, ?b) *
       ABB([s2]?, [y]?, ?c, [mu]?, ?d) =
    ABB([s1] + [s2], [x]*[y], ?a, ?b, ?c, ?d) * replace_([mu], DUMMY);
  also once ABB([s1]?, [x]?, ?a, [mu]?!fixed_, ?b, [mu]?, ?c) =
    ABB([s1], [x], ?a, ?b, ?c) * replace_([mu], DUMMY);
  sum DUMMY;
endrepeat;

id ABB(0, [x]?, ?a) = abbM([x])
#if "`Scale'" != "1"
  * MOM(?a)
#endif
  ;
id ABB([i]?, [x]?, ?a) = fmeM([x])
#if "`Scale'" != "1"
  * MOM(?a)
#endif
  ;

#if "`FermionChains'" != "Weyl"
repeat id fmeM([x]?) * fmeM([y]?) = fmeM([x] * [y]);
argument fmeM;
id e_([mu]?, [nu]?, [ro]?, [si]?) = Eps([mu], [nu], [ro], [si]);
endargument;
id fmeM([x]?) = Mat(fmeM([x]));
#endif

#if "`Scale'" != "1"
chainout MOM;
id MOM([p1]?MOMS) = SCALE;
id MOM(?a) = 1;
id SCALE^[x]? = powM(`Scale', [x]/2);
id powM([x]?, [y]?pos_) = [x]^[y];
#endif

b addM, mulM;
moduleoption polyfun=abbM;
.sort

normalize abbM;

id abbM(1) = 1;

*----------------------------------------------------------------------

#ifdef `PIPES_'
#include External.frm
#endif
#endprocedure

***********************************************************************

#procedure ToMma(expr)
#toexternal "%E,", `expr'
#endprocedure

***********************************************************************

#procedure FromMma(expr)
G `expr' =
#fromexternal
  ;
#endprocedure

***********************************************************************

#procedure MomSquare
id MOM(0) = 0;

#ifdef `MomSum'
* Apply momentum conservation to generate as few terms as possible

id MOM([p1]?) = MOM([p1], nterms_([p1]),
  [p1] + (`MomSum'), nterms_([p1] + (`MomSum')),
  [p1] - (`MomSum'), nterms_([p1] - (`MomSum')));

symm MOM (2,1) (4,3) (6,5);
#endif

id MOM([p1]?, ?a) = [p1].[p1];
#call kikj
#endprocedure

***********************************************************************

#procedure Fewest(foo)
argument `foo';
#call Neglect
endargument;
id `foo'([x]?, [y]?) = `foo'([x], nterms_([x]), [y], nterms_([y]));
symm `foo' (2,1), (4,3);
id `foo'([x]?, 1, ?a) = [x];
id `foo'([x]?, ?a) = `foo'([x]);
#endprocedure

***********************************************************************

#procedure Factor(foo)
factarg `foo';
chainout `foo';
id `foo'([x]?number_) = [x];
id `foo'([x]?symbol_) = [x];
#endprocedure

***********************************************************************

#procedure DoInsertions
.store

#call Insertions
#call Neglect

.sort
#endprocedure


***********************************************************************
*** the main program starts here

#if "`InsertionPolicy'" == "Begin"
#call DoInsertions
#define Inserted
#else
#call Neglect
#endif

#call Const
.sort

*----------------------------------------------------------------------

#define LoopInt "A0i, B0i, C0i, D0i, E0i, F0i"

#define CutInt "Acut, Bcut, Ccut, Dcut, Ecut, Fcut"

#define SUNObjs "SUNSum, SUNT, SUNTSum, SUNF, SUNEps"

#define Tensors "Spinor, GA, e_, eT1,...,eT`Legs', eTc1,...,eTc`Legs'"

#define Vectors "k1,...,k`Legs', e1,...,e`Legs', ec1,...,ec`Legs'"

* variables appearing in the CalcFeynAmp input and output
s I, D, Renumber;
cf Mat, Den, Pair, Eps, DiracChain, WeylChain, FormSimplify;
cf SumOver, IndexDelta, IndexEps, `SUNObjs';
cf `LoopInt';
cf `CutInt';
f Spinor, DottedSpinor;
i Col1,...,Col`Legs', Ind1,...,Ind10;

* variables that make it into Mma but don't appear in the output
cf addM, mulM, powM, intM, paveM, abbM, fmeM, sunM;
cf cutM, numM, qfM, qcM;

* patterns
s [x], [y], [z], [w], [n], [h];
s [k1], [k2], [k1k2];
s <[m0]>,...,<[m20]>;
s <[s0]>,...,<[s20]>;
v <[p0]>,...,<[p20]>;
i <[i0]>,...,<[i20]>;
i [mu], [nu], [ro], [si], [LA];
i [om], [omA], [omB];
i [i], [j], [k], [l], [I], [J], [K];
i [a], [b], [c], [d];
cf [f];
t [t];

* variables internal to FORM
s TAG, CUTRAT, SCALE;
i DUMMY, MUC, NUC, MUD, NUD;
cf TMP, MOM, ABB, ORD, NOW, CH, SIGN(antisymm);
t NUM, EQ, NEQ, DD, EPS(antisymm);
nt GA, GB, GC, GD;
f CC, WC;
auto s ARG;
set LOOPINT: `LoopInt';
set CUTINT: `CutInt';
set MOMS: k1,...,k`Legs';
set COLS: Col1,...,Col`Legs';

ntable BASIS(1:5, [mu]?, [nu]?);
ntable DUAL(1:5, [mu]?, [nu]?);
ntable CHI(6:7);

fill BASIS(1) = GA(6);
fill BASIS(2) = GA(7);
fill BASIS(3) = GA(6, [mu]);
fill BASIS(4) = GA(7, [mu]);
fill BASIS(5) = i_/2*(GA(6, [mu], [nu]) + GA(7, [mu], [nu]) -
                      (GA(6) + GA(7))*d_([mu], [nu]));

fill DUAL(1) = g_(1, 6_)/4;
fill DUAL(2) = g_(1, 7_)/4;
fill DUAL(3) = g_(1, 7_, [mu])/4;
fill DUAL(4) = g_(1, 6_, [mu])/4;
fill DUAL(5) = i_/4*(g_(1, [mu], [nu]) - d_([mu], [nu]));

fill CHI(6) = g6_(1)/2;
fill CHI(7) = g7_(1)/2;

.global

*----------------------------------------------------------------------

collect mulM;

moduleoption polyfun=mulM;
.sort

#call Square
argument;
#call Square
endargument;

id mulM(0) = 0;

normalize mulM;

id mulM(1) = 1;

id Eps([mu]?, [nu]?, [ro]?, [si]?) = e_([mu], [nu], [ro], [si]);


if( count(dirM, 1) ) ;

repeat;
  repeat id dirM([x]?, [i]?, [j]?) * dirM([y]?, [j]?, [k]?) =
    dirM([x]*[y], [i], [k]);
  id dirM([x]?, [i]?, [j]?) * dirM([y]?, [k]?, [j]?) =
    dirM([x]*CC([y]), [i], [k]);
  id dirM([x]?, [j]?, [i]?) * dirM([y]?, [j]?, [k]?) =
    dirM(CC([x])*[y], [i], [k]);
  argument dirM, 1;
    argument CC;
      id g_(sM, [mu]?) = -CC([mu]);
      chainin CC;
      id CC(?a) = CC(reverse_(?a));
      id g5_(sM) * CC(?a) = g_(sM, ?a, 5_);
      id g6_(sM) * CC(?a) = g_(sM, ?a, 6_);
      id g7_(sM) * CC(?a) = g_(sM, ?a, 7_);
      id CC(?a) = g_(sM, ?a);
    endargument;
    id CC([x]?) = [x];
  endargument;
endrepeat;

id dirM([x]?, [i]?) * dirM([y]?, [i]?, [j]?) *
     dirM(Spinor(?k, [s2]?), [j]?) =
  dirM([x]*[y]*Spinor(?k, -[s2]), [i], [j]);

repeat;
  once dirM([x]?, [i]?, [i]?) = -TMP([x]*replace_(sM, 0));
  id TMP([x]?) = [x];
  trace4, 0;
endrepeat;

#do fline = 1, 4
once dirM([x]?, [i]?, [j]?) = ORD([x]*replace_(sM, `fline'), [i], [j]);
#enddo

if( count(ORD, 1) );
  redefine HaveFermions "1";
  id ORD([x]?, ?a) = [x]*ORD(?a);
  chainin ORD;
  id ORD(?a) = SIGN(?a)*sign_(nargs_(?a)/2);
  mul replace_(SIGN, ORD);
  id ORD(?a) = 1;
endif;

endif;

*----------------------------------------------------------------------

b g_, `Tensors';
.sort

#if `HaveFermions' == 1

keep brackets;

id gi_([i]?) = g6_([i])/2 + g7_([i])/2;
id g5_([i]?) = g6_([i])/2 - g7_([i])/2;
id g_([i]?, [mu]?) = g_([i], 6_, [mu])/2 + g_([i], 7_, [mu])/2;

id g6_([i]?) = 2*GA(6);
id g7_([i]?) = 2*GA(7);
repeat id GA(?g) * g_([i]?, [mu]?) = GA(?g, [mu]);

#call DiracSimplify

#endif

*----------------------------------------------------------------------

b q1, intM;
.sort

keep brackets;

* cancel q^2's in the numerator

while( match(q1.q1) );
  once q1.q1 * intM(?a, Den(q1, [m1]?), ?b) =
    intM(?a, ?b) + [m1] * intM(?a, Den(q1, [m1]), ?b);
  also once intM(?a, Den([p1]?, 0), ?b) =
    replace_(q1, 2*q1 - [p1]) * intM(?a, Den([p1], 0), ?b);
  also once intM(Den([p1]?, [m1]?), ?a) =
    replace_(q1, 2*q1 - [p1]) * intM(Den([p1], [m1]), ?a);
endwhile;

id intM() = 1;
id intM(Den([p1]?, 0)) = 0;

*----------------------------------------------------------------------

b q1, intM;
.sort

keep brackets;

argument intM;
id Den([p1]?, [m1]?) = Den([m1])*MOM([p1]);
endargument;

symm intM;

id intM(Den([m0]?)*MOM([p0]?)) =
  replace_(q1, 2*q1 - [p0]) * intM(A0i(0), [m0]);

#do n = 1, 5
also intM(<Den([m0]?)*MOM([p0]?)>,...,<Den([m`n']?)*MOM([p`n']?)>) =
  replace_(q1, 2*q1 - [p0]) * (
#if "`CutTools'" != "True"
    ORD(<paveM(1)*([p1]-[p0])>+...+<paveM(`n')*([p`n']-[p0])>) *
    intM(LOOPINT[{`n'+1}](0),
#do i = 1, {(`n'+1)/2}
      <MOM([p`i']-[p0])>,...,<MOM([p`n']-[p{`n'-`i'}])>,
#if {2*`i'} <= `n'
      <MOM([p0]-[p{`n'-`i'+1}])>,...,<MOM([p{`i'-1}]-[p`n'])>,
#endif
#enddo
      <[m0]>,...,<[m`n']>)
#if "`CutTools'" == "Rational"
    * CUTRAT
#endif
#endif
#if "`CutTools'" != "False"
  + cutM(CUTINT[{`n'+1}],
      <[p1]-[p0]>,...,<[p`n']-[p0]>,
      <[m0]>,...,<[m`n']>)
#endif
  );
#enddo

*----------------------------------------------------------------------

#call DotSimplify(#call Contract)

*----------------------------------------------------------------------

b q1, NUM, ORD, Den, intM, `Tensors', CUTRAT;
.sort

keep brackets;

if( count(intM, 1) ) totensor q1, NUM;

#if "`Dim'" == "4"
* add local terms for dimred/CDR

id NUM([mu]?, [nu]?, [ro]?, [si]?) * intM(D0i(0), ?p) =
  NUM([mu], [nu], [ro], [si]) * intM(D0i(0), ?p) -
  5/144 * TMP(e_([mu], [nu], [ro], [si]), [mu], [nu], [ro], [si]) +
  1/8 * distrib_(1, 2, EQ, NEQ, [mu], [nu], [ro], [si]);

id EQ([mu]?, [mu]?) = 1;
id NEQ([mu]?, [mu]?) = 0;
id NEQ([mu]?, [nu]?) = d_([mu], [nu]);

id NUM([mu]?, [nu]?, [ro]?) * intM(C0i(0), MOM([p1]?), [x]?, MOM([p2]?), ?m) =
  NUM([mu], [nu], [ro]) * intM(C0i(0), MOM([p1]), [x], MOM([p2]), ?m) +
  1/36 * TMP(e_([mu], [nu], [ro]), [mu], [nu], [ro], [p1] - [p2]);

id TMP(0, ?a) = dd_(?a);
id TMP(?a) = 0;

id CUTRAT * intM(?a) = 0;
id CUTRAT = 1;

#endif


id Den([p1]?, ?m) = Den(MOM([p1]), ?m);

argument intM, Den;
#call MomSquare
endargument;

* decompose into Lorentz-covariant tensors

* The following statement introduces the g_{\mu\nu}'s in a smart way.
* Lifted from: S.A. Larin, T. van Ritbergen, and J.A.M. Vermaseren,
* The optimization of a huge FORM program,
* in: Proceedings Oberammergau 1993, ISBN 9-810-21699-8.

id NUM(?i) = sum_(DUMMY, 0, nargs_(?i), 2,
  paveM(0)^DUMMY * distrib_(1, DUMMY, dd_, NUM, ?i));

id ORD(0) * NUM([mu]?, ?i) = 0;
repeat id ORD([p1]?) * NUM([mu]?, ?i) = ORD([p1]) * d_([p1], [mu]) * NUM(?i);

id ORD(?p) = 1;
id NUM() = 1;

id intM(A0i(0), [m1]?) * NUM(?i) = 0;

*----------------------------------------------------------------------

b paveM, intM;
.sort

keep brackets;

if( count(paveM, 1) );
  chainin paveM;
  id paveM(?i) * intM([f]?(0), ?r) = intM([f](?i), ?r);
else;
  symm intM:4 3, 4;
endif;

*----------------------------------------------------------------------

#if `HaveFermions' == 1
* Dirac algebra on open fermion chains again

#if "`FermionChains'" == "Weyl"

.sort

* Chisholm's identity backwards to get rid of all Eps
*repeat id GA([om]?, ?a, [LA]?, ?b) * e_([mu]?, [nu]?, [ro]?, [LA]?) =
*  1/4 * sign_([om]) * sign_(nargs_(?a)) * (
*    GA([om], ?a, [mu], [nu], [ro], ?b) -
*    GA([om], ?a, [ro], [nu], [mu], ?b) );
repeat;
  once GA([om]?, ?a, [LA]?, ?b) * e_([mu]?, [nu]?, [ro]?, [LA]?) =
    sign_([om]) * sign_(nargs_(?a)) * (
      GA([om], ?a, [mu], [nu], [ro], ?b) -
      d_([mu], [nu]) * GA([om], ?a, [ro], ?b) +
      d_([mu], [ro]) * GA([om], ?a, [nu], ?b) -
      d_([nu], [ro]) * GA([om], ?a, [mu], ?b) );
endrepeat;

#elseif "`Dim'" == "4"

.sort

* this is Chisholm's identity:
repeat;
  once GA([om]?, [mu]?, [nu]?, [ro]?, ?a) =
    sign_([om]) * GA([om], DUMMY, ?a) * e_([mu], [nu], [ro], DUMMY) +
    d_([mu], [nu]) * GA([om], [ro], ?a) -
    d_([mu], [ro]) * GA([om], [nu], ?a) +
    d_([nu], [ro]) * GA([om], [mu], ?a);
  sum DUMMY;
endrepeat;

#call Contract

#endif

b `Tensors';
.sort

keep brackets;

#call DiracSimplify

#endif

*----------------------------------------------------------------------

#if "`Dim'" != 4
b D, Dminus4, intM, cutM, CUTRAT;
.sort

keep brackets;

id D = Dminus4 + 4;

#if "`Dim'" == "D"

id Dminus4 * cutM(?a) = 0;

* add local terms for dimreg
also Dminus4 * intM(A0i(0), [m1]?) = -2*[m1];
also Dminus4 * intM(A0i(0,0), [m1]?) = -[m1]^2/2;
also Dminus4 * intM(B0i(0), ?a) = -2;
also Dminus4 * intM(B0i(1), ?a) = 1;
also Dminus4 * intM(B0i(0,0), [k1]?, [m1]?, [m2]?) =
  1/6*([k1] - 3*[m1] - 3*[m2]);
also Dminus4 * intM(B0i(1,1), ?a) = -2/3;
also Dminus4 * intM(B0i(0,0,1), [k1]?, [m1]?, [m2]?) =
  -1/12*([k1] - 2*[m1] - 4*[m2]);
also Dminus4 * intM(B0i(1,1,1), ?a) = 1/2;
also Dminus4 * intM(C0i(0,0), ?a) = -1/2;
also Dminus4 * intM(C0i(0,0,[i]?), ?a) = 1/6;
also Dminus4 * intM(C0i(0,0,0,0), [k1]?, [k2]?, [k1k2]?, [m1]?, [m2]?, [m3]?) =
  1/48*([k1] + [k2] + [k1k2] - 4*([m1] + [m2] + [m3]));
also Dminus4 * intM(C0i(0,0,[i]?,[i]?), ?a) = -1/12;
also Dminus4 * intM(C0i(0,0,[i]?,[j]?), ?a) = -1/24;
also Dminus4 * intM(D0i(0,0,0,0), ?a) = -1/12;
also Dminus4 * intM(D0i(0,0,0,0,[i]?), ?a) = 1/48;

also Dminus4 = 0;

#redefine Dim "4"

#endif

id CUTRAT * intM(?a) = 0;
id CUTRAT = 1;

#endif

*----------------------------------------------------------------------

#if `HaveFermions' == 0

#call Abbreviate

#else

#if "`FermionChains'" == "Weyl"

b Spinor;
.sort

keep brackets;

id Spinor(?a) * GA(?g) * Spinor(?b) =
  CH(DottedSpinor(?a), ?g, Spinor(?b));

repeat;
  once CH([s1]?, [x]?, ?a, [LA]?, ?b) *
       CH([s2]?, [y]?, ?c, [LA]?, ?d) =
    WC(sign_(nargs_(?a, ?c) + [x] + [y]), [s1], [x], ?a) *
    WC(?b) * WC([s2], [y], ?c) * WC(?d);

* Fierz 1: <A|sig_mu|B> <C|sigbar^mu|D> = 2 <A|D> <C|B>
  id WC(-1, ?a) * WC(?b) * WC(?c) * WC(?d) =
    2*CH(?a, ?d) * CH(?c, ?b);

* Fierz 2: <A|sig(bar)_mu|B> <C|sig(bar)^mu|D> = 2 <A|eps|C> <B|eps|D>
  also WC(1, ?a) * WC(?b, [s1]?) * WC([s2]?, [x]?, ?c) * WC(?d) =
    2 * TMP(?c, ?b) * CH(?a, reverse_(?c), -1, [s2]) *
      CH([s1], 6 + mod_([x] + nargs_(?c, ?b), 2), -1, reverse_(?b), ?d);
  chainout TMP;
  id TMP(-1) = -1;
  id TMP(?a) = 1;

  id CH(?a, -1, -1, ?b) = -CH(?a, ?b);

* due to the canonical ordering of the Dirac chains this
* is the only(?) case we need of Fierz on the same chain:
  repeat id CH(?a, [LA]?, [LA]?, ?b) = 4*CH(?a, ?b);
endrepeat;

id CH([s1]?, ?g, [s2]?) = abbM(fmeM(WeylChain([s1], ?g, [s2])))
#if "`Scale'" != "1"
  * MOM(?g)
#endif
  ;

#call Abbreviate

#else

#switch "`FermionOrder'"

#case "None"
.sort
#call DiracFinal
#break

#case "Fierz"
* Fierz twice for simplification
#call FierzUnordered
#call FierzUnordered
#call DiracFinal
#break

#case "Automatic"
* lexicographical ordering
#call FierzPre(1)
#call FierzOrdered
#call DiracFinal
#break

#case "Colour"
* postponed until after SUNT simplification
#break

#default
#do i = {`FermionOrder'}
#ifndef `order'
#define order "k{`i'}"
#else
#redefine order "`order'\,k{`i'}"
#endif
#enddo
mul ORD(`order');
#call FierzPre(1)
#call FierzOrdered
#call DiracFinal

#endswitch

#endif

#endif

*----------------------------------------------------------------------

#ifdef `Inserted'
.sort
#else
#call DoInsertions
#define Inserted
#endif

id intM(A0i(?a), 0) = 0;

#call Const
.sort

collect mulM;

argument mulM;
#call Neglect
endargument;

argument;
#call Square
endargument;

id mulM(0) = 0;

*----------------------------------------------------------------------
* index handling

.sort

repeat;
  once SumOver([i]?, ?a, Renumber) =
    TMP(DUMMY) * SumOver(DUMMY, ?a) * replace_([i], DUMMY);
  sum DUMMY;
endrepeat;

id IndexEps([i]?, [j]?, [k]?) = EPS([i], [j], [k]);

repeat;
  id EPS([I]?, [J]?, [K]?) * EPS([I]?, [J]?, [K]?) *
    SumOver([I]?, 3) * SumOver([J]?, 3) * SumOver([K]?, 3) = 6;
  id EPS([I]?, [J]?, [k]?) * EPS([I]?, [J]?, [c]?) *
    SumOver([I]?, 3) * SumOver([J]?, 3) = 2*IndexDelta([k], [c]);
  id EPS([I]?, [j]?, [k]?) * EPS([I]?, [b]?, [c]?) *
    SumOver([I]?, 3) =
    IndexDelta([j], [b])*IndexDelta([k], [c]) -
    IndexDelta([j], [k])*IndexDelta([b], [c]);
  repeat;
    id IndexDelta([I]?, [I]?) = 1;
    symm IndexDelta;
    once ifmatch->1 IndexDelta([i]?, [J]?) * SumOver([J]?, [x]?) =
      replace_([J], [i]);
    once IndexDelta([I]?, [j]?) * SumOver([I]?, [x]?) =
      replace_([I], [j]);
    label 1;
  endrepeat;
endrepeat;

id IndexDelta([x]?int_, [y]?int_) = delta_([x], [y]);

id TMP([x]?int_) = 1;
repeat id TMP([I]?)^2 = TMP([I]);

#do i = 1, 10
once TMP([I]?) = replace_([I], Ind`i');
#enddo

id EPS([i]?, [j]?, [k]?) = IndexEps([j], [k], [i]);

*----------------------------------------------------------------------

#if `HaveSUN' == 1
* simplification of SU(N) structures

* The algorithm implemented here is an extension of the one given in
* J.A.M. Vermaseren, The use of computer algebra in QCD,
* in: Proceedings Schladming 1996, ISBN 3-540-62478-3.

* The idea is to transform all SU(N) objects to generators, SUNT.
* In the output, only two types of objects can appear:
* - chains of SUNTs (with external colour indices), or
* - traces of SUNTs.
* A chain of SUNTs is denoted by SUNT(a, b, ..., i, j), where
* a, b, ... are gluon indices and i and j are colour indices.
* SUNT(i, j) is the special case of the identity in colour space.
* A trace over SUNTs is marked by both colour indices being zero,
* i.e. SUNT(a, b, ..., 0, 0).

b `SUNObjs';
.sort

keep brackets;

if( count(SUNF, 1) );

  repeat;
    once SUNF(?a, [a]?, [b]?, [c]?, [d]?) =
      SUNF(?a, [a], [b], DUMMY) * SUNF(DUMMY, [c], [d]) * SUNSum(DUMMY);
    sum DUMMY;
  endrepeat;

* f^{abc} = 2 i Tr(T^c T^b T^a - T^a T^b T^c)

  id SUNF([a]?, [b]?, [c]?) =
    2*i_*(SUNT([c], [b], [a], 0, 0) - SUNT([a], [b], [c], 0, 0));

endif;


repeat;
  once SUNT(?a, 0, 0) = SUNT(?a, DUMMY, DUMMY) * SUNSum(DUMMY);
  sum DUMMY;
endrepeat;

repeat;
  once SUNT(?a, [a]?, [b]?, [i]?, [j]?) =
    SUNT(?a, [a], [i], DUMMY) * SUNT([b], DUMMY, [j]) * SUNSum(DUMMY);
  sum DUMMY;
endrepeat;


* T^a_{ij} T^a_{kl} =
*   1/2 (delta_{il} delta_{jk} - 1/N delta_{ij} delta_{kl})

id SUNT([a]?, [i]?, [j]?) * SUNT([a]?, [k]?, [l]?) * SUNSum([a]?, ?a) =
  1/2 * SUNT([i], [l]) * SUNT([j], [k]) -
  1/2/(`SUNN') * SUNT([i], [j]) * SUNT([k], [l]);

id SUNTSum([i]?, [j]?, [k]?, [l]?) =
  1/2 * SUNT([i], [l]) * SUNT([j], [k]) -
  1/2/(`SUNN') * SUNT([i], [j]) * SUNT([k], [l]);


id SUNEps([i]?, [j]?, [k]?) = EPS([i], [j], [k]);


* cleaning up, step 1: get rid of the deltas

repeat;
  id EPS([I]?, [j]?, [k]?) * EPS([I]?, [b]?, [c]?) *
    SUNSum([I]?, [x]?) = 
    SUNT([j], [b])*SUNT([k], [c]) -
    SUNT([j], [k])*SUNT([b], [c]);
  repeat;
    id SUNT([I]?, [I]?) * SUNSum([I]?, ?a) = `SUNN';
    symm SUNT:2 1, 2;
    once ifmatch->1 SUNT([I]?, [j]?) * SUNSum([I]?, ?a) = replace_([I], [j]);
    once SUNT([i]?, [J]?) * SUNSum([J]?, ?a) = replace_([J], [i]);
    label 1;
  endrepeat;
endrepeat;

id SUNT([x]?int_, [y]?int_) = delta_([x], [y]);
id SUNT([a]?, [i]?, [i]?) * SUNSum([i]?, ?a) = 0;

id EPS([i]?, [j]?, [k]?) = sunM(SUNEps([j], [k], [i]));

* cleaning up, step 2: bead up the SUNTs into chains

repeat;
  once SUNT(?a, [a]?, [i]?, [j]?) = TMP(?a, [a], [i], [j]);
  repeat;
    id TMP(?a, [i]?, [j]?) * SUNT(?b, [j]?, [k]?) * SUNSum([j]?, ?c) =
      TMP(?a, ?b, [i], [k]);
    id SUNT(?a, [i]?, [j]?) * TMP(?b, [j]?, [k]?) * SUNSum([j]?, ?c) =
      TMP(?a, ?b, [i], [k]);
  endrepeat;

  id TMP(?a, [i]?, [i]?) * SUNSum([i]?, ?b) = TMP(?a, 0, 0);

* special case of Tr(T^a T^b) = 1/2 delta_{ab}
*  id TMP([a]?, [a]?, 0, 0) = 1/2;
  id TMP([x]?int_, [y]?int_, 0, 0) = 1/2*delta_([x], [y]);

  id TMP(?a) = sunM(SUNT(?a));
endrepeat;

symm SUNT;

id SUNT(?a) = sunM(SUNT(?a));

#if "`FermionOrder'" == "Colour"
id sunM(SUNT(?a, [i]?COLS[[x]], [j]?COLS[[y]])) =
  sunM(SUNT(?a, [i], [j])) * ORD(MOMS[[x]], MOMS[[y]]);
#endif

repeat id sunM([x]?) * sunM([y]?) = sunM([x] * [y]);

id SUNSum([i]?, [x]?) = [x];

* the Mat(...) are kept at the almost outermost level (only SumOver
* comes before), i.e. the amplitude is of the form Sum[c[i] Mat[i], i];
* this is required for the calculation of the squared amplitude

id sunM([x]?) = Mat(sunM([x]));

#if "`FermionOrder'" == "Colour"
#call FierzPre(1)
#call FierzOrdered
#call DiracFinal
#endif

#endif

*----------------------------------------------------------------------

.sort

mul replace_(i_, I);
id I = mulM(I);

repeat id Mat([x]?) * Mat([y]?) = Mat([x] * [y]);

*----------------------------------------------------------------------

moduleoption polyfun=mulM;
.sort

normalize mulM;

moduleoption polyfun=abbM;
.sort

normalize abbM;
id abbM(1) = 1;

*----------------------------------------------------------------------

b SumOver, Mat, Den, abbM, mulM, intM, cutM, qfM, Dminus4;
.sort

collect FormSimplify;

moduleoption polyfun=FormSimplify;
.sort

normalize FormSimplify;

#call Factor(FormSimplify)
#call InvSimplify(FormSimplify)

id FormSimplify(?a) = TMP(FormSimplify(?a));
id abbM(?a) = TMP(abbM(?a));

moduleoption polyfun=TMP;
.sort

*----------------------------------------------------------------------

normalize TMP;
id TMP(1) = 1;
id TMP(?a) = dum_(?a);

moduleoption polyfun=mulM;
.sort

normalize mulM;
id mulM(1) = 1;
factarg mulM;

*----------------------------------------------------------------------

#if "`CutTools'" != "False"

b cutM, intM, qfM, Dminus4, SumOver, Mat, Den;
.sort

collect FormSimplify;

b cutM, intM, qfM, Dminus4, SumOver, Mat;
.sort

collect qcM;
normalize qcM;

id qfM([x]?) = [x];

b cutM, intM, SumOver, Mat;
.sort

collect numM, numM;
normalize numM;

id cutM([f]?, ?a) * numM([x]?) = [f](numM([x]), ?a);

id numM([x]?) = [x];

#endif

*----------------------------------------------------------------------

mul replace_(intM, paveM);

b SumOver, Mat, Den, paveM;
print;

.end

