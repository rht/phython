* PolarizationSum.frm
* the FORM part of the PolarizationSum function
* this file is part of FormCalc
* last modified 19 May 10 th


#procedure PolSum(i, m)
#$dim = 4;
if( count(z`i',1, zc`i',1) ) $dim = Dminus4;

b z`i', zc`i', e`i', ec`i', eT`i', eTc`i';
.sort
d `$dim';

keep brackets;

id e`i' = ET(?);
id ec`i' = ETC(?);
id z`i' = ET(?);
id zc`i' = ETC(?);

#if `m' == "0"
* massless case

multiply 2;
id ET([mu]?) * ETC([nu]?) = 1/2*( -d_([mu], [nu])
    + (d_(eta`i', [mu])*d_(k`i', [nu]) +
       d_(eta`i', [nu])*d_(k`i', [mu]))/(eta`i'.k`i') );
* The eta are gauge-dependent vectors.
* Their appearance in the result is supposed to alert
* the user to the presence of gauge-dependent terms.
* The eta must fulfill eta.eta = e.eta = 0 and k.eta != 0.
* Instead of imposing eta.eta = 0 one can add
* - d_(k`i', [mu])*d_(k`i', [nu])*(eta`i'.eta`i')/(eta`i'.k`i')^2

id eT`i'([mu]?, [nu]?) * eTc`i'([ro]?, [si]?) = 1/2*(
  d_([mu], [ro])*d_([nu], [si]) +
  d_([mu], [si])*d_([nu], [ro]) -
  d_([mu], [nu])*d_([ro], [si]) );

#else
* massive case

multiply 3;
id ET([mu]?) * ETC([nu]?) = 1/3*( -d_([mu], [nu]) +
  k`i'([mu])*k`i'([nu])/(`m')^2 );

#endif

.sort

#call eiei
#call eiki
#call kikj
#call Neglect
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
repeat id `foo'([x]?, [y]?, ?a) = `foo'([x]) * `foo'([y], ?a);
id `foo'([x]?number_) = [x];
id `foo'([x]?symbol_) = [x];
#endprocedure

***********************************************************************

#procedure Emit
contract 0;

#do i = 1, `Legs'
#ifdef `k`i''
id k`i' = `k`i'';
.sort;
#endif
#enddo

#call kikj
#call Neglect

#if `GaugeTerms' == 0
#do i = 1, `Legs'
id eta`i' = 0;
#enddo
#endif

.sort

id [p1]?.[p2]? = abbM([p1].[p2], [p1], [p2]);

id 1/[p1]?.[p2]? = 1/abbM([p1].[p2])
#if "`Scale'" != "1"
  * IMOM([p1], [p2])
#endif
  ;

id e_([mu]?, [nu]?, [ro]?, [si]?) =
  abbM(e_([mu], [nu], [ro], [si]), [mu], [nu], [ro], [si]);

id d_([mu]?, [nu]?) = abbM(d_([mu], [nu]), [mu], [nu]);

id [t]?(?a) = abbM([t](?a), ?a);

id [p1]?([mu]?) = abbM([p1]([mu]), [p1]);

repeat;
  once abbM([x]?, ?a, [mu]?!fixed_, ?b) *
       abbM([y]?, ?c, [mu]?, ?d) =
    abbM([x]*[y], ?a, ?b, ?c, ?d) * replace_([mu], DUMMY);
  also once abbM([x]?, ?a, [mu]?!fixed_, ?b, [mu]?, ?c) =
    abbM([x], ?a, ?b, ?c) * replace_([mu], DUMMY);
  sum DUMMY;
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
chainout IMOM;
id IMOM([p1]?MOMS) = 1/SCALE;
id IMOM(?a) = 1;
id SCALE^[x]? = powM(`Scale', [x]/2);
#endif

#call Square

b abbM, `Bracket';
.sort

collect FormSimplify, FormSimplify;
normalize FormSimplify;

#call Factor(FormSimplify)
#call InvSimplify(FormSimplify)

id FormSimplify(0) = 0;

.sort

moduleoption polyfun=abbM;
.sort

normalize abbM;
id abbM(1) = 1;

b abbM, `Bracket';
print;
.end
#endprocedure

***********************************************************************

#define Bracket "Den, A0, A00, B0, B1, B00, B11, B001, B111, A0i, B0i, C0i, D0i, E0i, F0i"

i DUMMY;
cf MOM, IMOM;
s SCALE;
set MOMS: k1,...,k`Legs';
auto s ARG;

s Dminus4;
i [mu], [nu], [ro], [si];
v [p1], [p2];
s [x], [y];
t [t];

cf abbM, powM, FormSimplify;
cf `Bracket';
t ET, ETC;

.global

