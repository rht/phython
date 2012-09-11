* HelicityME.frm
* the FORM part of the HelicityME function
* this file is part of FormCalc
* last modified 6 Apr 11 th


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
repeat id `foo'([x]?, [y]?, ?a) = `foo'([x]) * `foo'([y], ?a);
id `foo'([x]?number_) = [x];
id `foo'([x]?symbol_) = [x];
#endprocedure

***********************************************************************

#procedure Emit

id DiracChain(Spinor(?p), [x]?pos_, ?g, Spinor(?q)) =
  DiracChain(Spinor(?p), 1, [x], ?g, Spinor(?q));

* un-antisymmetrize the Dirac chains if necessary
also DiracChain(Spinor(?p), [x]?, ?g, Spinor(?q)) =
  Spinor(?p) * GM(-[x]) * ( GD(?g) -
    sum_(KK, 2, nargs_(?g), 2, distrib_(-1, KK, DD, GD, ?g)) ) *
  Spinor(?q);

id DD() = 1;
repeat;
  once DD(?a) = g_(1, ?a)/4;
  trace4, 1;
endrepeat;

id Spinor(?p) * GM([om]?) * GD(?g) * Spinor(?q) =
  DiracChain(Spinor(?p), 1, [om], ?g, Spinor(?q));

* The explicit 1 above is make each chirality projector count as
* two, such that sign_(nargs_(.)) effectively ignores the projector.

.sort

repeat;
  id DiracChain(?a, Spinor(?p)) * DiracChain(Spinor(?p), ?b) =
    DiracChain(?a, RHO(?p), ?b);

* If the spinors at the ends don't match directly, i.e.
*   <s2| g1 g2... |s1> <s2| ... |>,
* we use charge conjugation to reverse the first chain to have the
* |s1>'s side by side for substituting the projector (|s2><s2|).
* Inserting 1 = C C^-1 results in
*   <s2|C (C^-1 g1 C) (C^-1 g2 C) ... C^-1|s1>
*   = <anti-s2| (-g1)^T (-g2)^T ... |anti-s1>
*   = (-1)^(# gammas) <anti-s2| (... g2 g1)^T |anti-s1>
*   = (-1)^(# gammas) <anti-s1| ... g2 g1 |anti-s2>
* Thus follow the rules:
*   a) reverse the chain and exchange u <-> v,
*   b) gamma_mu -> -gamma_mu.
*   c) add a global minus sign to compensate for the
*      change in the permutation of the external fermions.
* For more details see the Denner/Eck/Hahn/Kueblbeck paper.
* Note that RHO and RHOC are counted as gamma matrices towards
* the overall sign; this is corrected in the RHOC substitution
* later.

  id DiracChain(Spinor([p2]?, ?p), ?a, Spinor([p1]?, [m1]?, [s1]?)) *
     DiracChain(Spinor([p2]?, ?q), ?b) =
    -sign_(nargs_(?a)) *
    DiracChain(Spinor([p1], [m1], -[s1]),
      reverse_(?a)*replace_(RHO, RHOC, RHOC, RHO),
      RHO([p2], ?q), ?b);
endrepeat;

repeat;
  once DiracChain(Spinor(?p), ?g, Spinor(?p)) = RHO(?p) * GM(?g);
  chainout GM;

  id GM([mu]?) = CHI([mu]);
  id CHI([mu]?) = g_(1, [mu]);
  id GM([x]?) = [x];

  argument RHO, RHOC;
#call Neglect
  endargument;

  id RHO([p1]?MOMS[[x]], 0, [s1]?) =
    [s1]/4*(g6_(1)*HEL([x], [s1]) - g7_(1)*HEL([x], -[s1])) *
      g_(1, [p1]);

  id RHOC([p1]?MOMS[[x]], 0, [s1]?) =
    g_(1, [p1]) *
      [s1]/4*(g6_(1)*HEL([x], [s1]) - g7_(1)*HEL([x], -[s1]));

  id RHO([p1]?MOMS[[x]], [m1]?, [s1]?) =
    (g_(1) + HEL([x], 0)*g_(1, 5_, EPSS[[x]]))/2 *
      (g_(1, [p1]) + [s1]*[m1]*g_(1));

  id RHOC([p1]?MOMS[[x]], [m1]?, [s1]?) =
    (g_(1, [p1]) - [s1]*[m1]*g_(1)) *
      (g_(1) - HEL([x], 0)*g_(1, EPSS[[x]], 5_))/2;

  trace4, 1;
endrepeat;

id D = Dminus4 + 4;

contract;

#call eiei

#do i = 1, `Legs'
#ifdef `k`i''
#call eiki
id k`i' = `k`i'';
.sort;
#endif
#enddo

#call eiki
#call kikj

#call Neglect

id D = Dminus4Eps + 4;

.sort

id [p1]?.[p2]? = abbM([p1].[p2], [p1], [p2]);

id e_([mu]?, [nu]?, [ro]?, [si]?) =
  abbM(e_([mu], [nu], [ro], [si]), [mu], [nu], [ro], [si]);

id d_([mu]?, [nu]?) = abbM(d_([mu], [nu]), [mu], [nu]);

id [t]?(?a) = abbM([t](?a), ?a);

id [p1]?([mu]?) = abbM([p1]([mu]), [p1]);

repeat;
  once abbM([x]?, ?a, [mu]?!fixed_, ?b) *
       abbM([y]?, ?c, [mu]?, ?d) =
    abbM([x]*[y], ?a, ?b, ?c, ?d) * replace_([mu], N100_?);
  also once abbM([x]?, ?a, [mu]?!fixed_, ?b, [mu]?, ?c) =
    abbM([x], ?a, ?b, ?c) * replace_([mu], N100_?);
  renumber;
endrepeat;

id abbM([x]?, ?a) = abbM([x])
#if "`Scale'" != "1"
  * MOM(?a)
#endif
  ;

#if "`Scale'" != "1"
chainout MOM;
id MOM([p1]?MOMS) = SCALE;
id MOM(?a) = 1;
id SCALE^[x]? = powM(`Scale', [x]/2);
#endif

#call Square

moduleoption polyfun=abbM;
.sort

normalize abbM;
id abbM(1) = 1;

b helM, abbM;
.sort

collect FormSimplify, FormSimplify;
normalize FormSimplify;

#call Factor(FormSimplify)
#call InvSimplify(FormSimplify)

b helM;
print;
.store
#endprocedure

***********************************************************************

i [om], [mu], [nu], [ro], [si];
v [p1], [p2];
s [m1], [m2], [s1], [s2], [x], [y];
t [t];

i KK, DUMMY;
t DD;
nt GD;
f RHO, RHOC, GM;
cf MOM;
s SCALE;
auto s ARG;
set MOMS: k1,...,k`Legs';
set EPSS: e1,...,e`Legs';

s D, Dminus4, Dminus4Eps;
cf abbM, helM, powM, DiracChain, FormSimplify;

ntable CHI(0:7);
fill CHI(0) = g_(1);
fill CHI(1) = g_(1);
fill CHI(4) = -g5_(1);
fill CHI(5) = g5_(1);
fill CHI(6) = g6_(1)/2;
fill CHI(7) = g7_(1)/2;

.global

