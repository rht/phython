(*
	VecSet.m
		explicit representation of external vectors
		and derived quantities
		this file is part of FormCalc
		last modified 15 Jun 10 th
*)


BeginPackage["VecSet`", "FormCalc`"]

VecSet::usage = "VecSet[i, m, p, e] defines the vectors for external
particle i with mass m and three-momentum p moving in direction
e = {ex, ey, ez}."

ToComponents::usage = "ToComponents[expr, pol] substitutes the vectors
appearing in expr by their components as set by VecSet.  The argument
pol specifies the polarizations of the external particles, they are
either given as a list of integers (-1 = left-handed, +1 = right-handed,
0 = longitudinal) or as a string of \"+\", \"-\", \"0\"."

{ VecK, VecE, MomSpec, Spi, DotSpi,
  SpecM, SpecK, SpecE, SpecDeltaK, SpecEx, SpecEy, SpecEz,
  SpecPhi, SpecKT, SpecET, SpecKT, SpecPRap, SpecRap }


Begin["`Private`"]

cp[{r_, i_}] := r + I i

cm[{r_, i_}] := r - I i


safelog[n_, d_] := Infinity /; n == 0. || d == 0.

safelog[n_, d_] := Log[n/d]


(* i: the index of the momentum
   m, p: mass and three-momentum of the particle
   ex,ey,ez: the unit three-vector of the momentum *)

VecSet[i_, m_, p_, {ex_,ey_,ez_}, simp_:Simplify] :=
Block[ {p0, deltap, eT, sinth, expIphi, phi, onePez, oneMez},
  If[ TrueQ[m < 10^-7],
    p0 = p;
    deltap = 0,
  (* else *)
    p0 = Sqrt[p^2 + m^2];
    deltap = m^2/(p0 + p) ];

  eT = {ex, ey};
  sinth = eT.eT;
  onePez = 1 + ez;
  oneMez = 1 - ez;

  If[ TrueQ[sinth < 10^-14],
	(* phi is irrelevant when theta = 0 *)
    expIphi = {1, 0};
    phi = 0,
  (* else *)
    sinth = Sqrt[sinth];
    expIphi = eT/sinth;
    phi = ArcTan@@ eT ];

  VecK[i] = {{p0 onePez - deltap ez, p cm[eT]},
             {p cp[eT], p0 oneMez + deltap ez}} //simp;

  If[ m =!= 0,
    VecE[i, 0] = {{(p onePez + deltap ez)/m, p0/m cm[eT]},
                  {p0/m cp[eT], (p oneMez - deltap ez)/m}} //simp ];

  VecE[i, -1] = 1/Sqrt[2] {{-sinth, -oneMez cm[expIphi]},
                           {onePez cp[expIphi], sinth}} //simp;
  VecE[i, +1] = 1/Sqrt[2] {{-sinth, onePez cm[expIphi]},
                           {-oneMez cp[expIphi], sinth}} //simp;

  MomSpec[i] = {
    SpecM -> m,
    SpecK -> p,
    SpecE -> p0,
    SpecDeltaK -> deltap,
    SpecEx -> ex,
    SpecEy -> ey,
    SpecEz -> ez,
    SpecKT -> p sinth,
    SpecET -> Sqrt[(p sinth)^2 + m^2],
    SpecPhi -> phi,
    SpecPRap -> safelog[onePez, sinth],
    SpecRap -> safelog[VecK[i][[1,1]], Sqrt[(p sinth)^2 + m^2]]
  };

  (* this is E^[I phi] cos[th/2] = 1/Sqrt[2] Sqrt[1 + ez] expIphi: *)
  expIphi *= 1/Sqrt[2] Sqrt[onePez];

  (* this is sin[th/2]: *)
  sinth = 1/Sqrt[2] Sqrt[oneMez];

  sump = Sqrt[p0 + p];
  deltap = Sqrt[deltap];

  Spi[i, -1, 6] = deltap {sinth, -cp[expIphi]} //simp;
  DotSpi[i, -1, 6] = deltap {sinth, -cm[expIphi]} //simp;

  Spi[i, -1, 7] = sump {sinth, -cp[expIphi]} //simp;
  DotSpi[i, -1, 7] = sump {sinth, -cm[expIphi]} //simp;

  Spi[i, 1, 6] = sump {cm[expIphi], sinth} //simp;
  DotSpi[i, 1, 6] = sump {cp[expIphi], sinth} //simp;

  Spi[i, 1, 7] = deltap {cm[expIphi], sinth} //simp;
  DotSpi[i, 1, 7] = deltap {cp[expIphi], sinth} //simp;

  True
]


ToComponents[expr_, pol_String] :=
  ToComponents[expr, Characters[pol] /. {"+" -> 1, "-" -> -1, "0" -> 0}]

ToComponents[expr_, pol_List] := expr /.
  WeylChain -> SplitChain /.
  { Spinor[i_, s_, om_] -> Spi[i, s Hel[i], om],
    DottedSpinor[i_, s_, om_] -> DotSpi[i, s Hel[i], om] } /.
  Flatten[MapIndexed[vecrule, pol]] /.
  h:SxS | SeS | VxS | VeS | BxS | BeS | Pair | Eps :> rep[h] /.
  Hel[i_] :> pol[[i]]


vecrule[pol_, {n_}] := {
  k[n] -> VecK[n],
  e[n] -> VecE[n, pol],
  ec[n] -> VecE[n, -pol] }


rep[Pair][A_List, B_List] :=
  1/2 (A[[1,1]] B[[2,2]] + A[[2,2]] B[[1,1]] -
       A[[1,2]] B[[2,1]] - A[[2,1]] B[[1,2]])

diag[i_, j_] := i[[1,1]] j[[2,2]] - i[[2,2]] j[[1,1]]

offdiag[i_, j_] := i[[2,1]] j[[1,2]] - i[[1,2]] j[[2,1]]

rep[Eps][A_List, B_List, C_List, D_List] := 1/4 (
  diag[A, B] offdiag[C, D] +
  diag[C, D] offdiag[A, B] -
  diag[A, C] offdiag[B, D] -
  diag[B, D] offdiag[A, C] +
  diag[A, D] offdiag[B, C] +
  diag[B, C] offdiag[A, D] )


eps[{b1_, b2_}] := {b2, -b1}

bar[{{A11_, A12_}, {A21_, A22_}}] := {{A22, -A12}, {-A21, A11}}

rep[SxS][a_List, b_List] := a.b

rep[SeS][a_List, b_List] := a.eps[b]

rep[VxS][A_List, b_List] := A.b

rep[VeS][A_List, b_List] := A.eps[b]

rep[BxS][A_List, b_List] := bar[A].b

rep[BeS][A_List, b_List] := bar[A].eps[b]


End[]

EndPackage[]

