(*
	bfunc.m
		explicit formulas for the one- and two-point
		functions and their derivatives
		this file is part of FormCalc
		last modified 30 Sep 09 th

The formulas are more or less directly from the Passarino-Veltman
paper. The regularization parameters are
	Mudim  -- the mass scale squared
	Lambda -- the photon mass squared
	Delta  -- the infinity -2/(4 - D) - EulerGamma + Log[4 Pi]
*)


F[n_, x_] := -x^n Log[(x - 1)/x] - Sum[x^(n - m)/m, {m, n}]

Block[ {x, theroots},
  theroots = Solve[-p x^2 + (p + m1 - m2) x - m1 == 0, x];
  {x1, x2} = (x /. theroots) + {I eps, -I eps};
]

N[eps, ___] = 10^-50

N[Delta, ___] = 0

N[Mudim, ___] = 1

N[Lambda, ___] = 1


A0[0] = 0

A0[m_] = m (1 - Log[m/Mudim] + Delta)


B0[0, 0, 0] = BAD	(* divergent, must cancel *)

B0[p_, 0, 0] = 2 - Log[-p/Mudim - I eps] + Delta

B0[p_, m1_, 0] := B0[p, 0, m1]

B0[0, m_, m_] = -Log[m/Mudim] + Delta

B0[0, m1_, m2_] = -Log[m2/Mudim] - Simplify[F[1, m1/(m1 - m2)]] + Delta

B0[p_, m1_, m2_] = -Log[m2/Mudim] - Simplify[F[1, x1] + F[1, x2]] + Delta


B1[p_, 0, 0] := -1/2 B0[p, 0, 0];

B1[p_, m1_, 0] := -B1[p, 0, m1] - B0[p, 0, m1]

B1[0, m_, m_] = (Log[m/Mudim] - Delta)/2

B1[0, m1_, m2_] = (Log[m2/Mudim] + Simplify[F[2, m1/(m1 - m2)]] - Delta)/2

B1[p_, m1_, m2_] = (Log[m2/Mudim] + Simplify[F[2, x1] + F[2, x2]] - Delta)/2


B11[p_, m1_, 0] := B11[p, 0, m1] -
  (A0[m1] - (m1 + 2 p) B0[p, 0, m1] - 4 p B1[p, 0, m1])/(3 p)

B11[0, m1_, m2_] = (-Log[m2/Mudim] - F[3, m1/(m1 - m2)] + Delta)/3

B11[p_, m1_, m2_] = (-Log[m2/Mudim] - Simplify[F[3, x1] + F[3, x2]] + Delta)/3


B111[p_, m1_, m2_] = (Log[m2/Mudim] + Simplify[F[4, x1] + F[4, x2]] - Delta)/4


B00[p_, m1_, m2_] :=
  (3 (m1 + m2) - p)/18 + m1 B0[p, m1, m2]/3 +
  (A0[m2] + (m1 - m2 + p) B1[p, m1, m2])/6


Derivative[1, 0, 0][B0] = DB0

DB0[m_, 0, m_] = DB0[m_, m_, 0] = -(1 + Log[Lambda/m]/2)/m

DB0[0, m_, m_] = 1/6/m

DB0[0, m1_, m2_] = (1/2 + m2/(m2 - m1) F[1, m1/(m1 - m2)])/(m1 - m2)

DB0[p_, m1_, m2_] = ((1 - x2) F[1, x2] - (1 - x1) F[1, x1])/p/(x1 - x2)


Derivative[1, 0, 0][B1] = DB1

DB1[m_, m_, 0] = (3 + Log[Lambda/m])/2/m

DB1[0, m_, m_] = -1/12/m

DB1[0, m1_, m2_] = -(1/3 + m2/(m2 - m1) F[2, m1/(m1 - m2)])/(m1 - m2)

DB1[p_, m1_, m2_] = ((1 - x1) F[2, x1] - (1 - x2) F[2, x2])/p/(x1 - x2)


Derivative[1, 0, 0][B11] = DB11

DB11[p_, m1_, m2_] = ((1 - x2) F[3, x2] - (1 - x1) F[3, x1])/p/(x1 - x2)


Derivative[1, 0, 0][B00] = DB00

DB00[p_, m1_, m2_] := -1/18 + m1 DB0[p, m1, m2]/3 + 
 (B1[p, m1, m2] + (p + m1 - m2) DB1[p, m1, m2])/6

(* note: DB00 looks like it could be IR divergent --
   it's not: the lambda-dependence cancels :-) *)


Null

