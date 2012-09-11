(*
	Convert2to3.m
		adapts FeynArts 2.2 GraphInfo files to the
		new, canonically sorted topologies of FA 3
		last modified 23 Jan 03 th
*)


tosort[ Vertex[1][n_] ] := (vc = Max[vc, n]; vert[1][n, 0])

tosort[ Vertex[e_, c_:0][n_] ] := vert[e][1000 + n, c]

tosort[ Propagator[Loop[n_]][v__] ] :=
  prop[Loop, psort[v], {Loop[1000 + n]}]

tosort[ Propagator[type_][v__] ] := prop[type, psort[v]]

fromv[ vert[e_][n_, _] ] := Vertex[e][n]

(* fromv[ vert[e_][n_, c_] ] := Vertex[e, c][n] *)

fromp[ t_, v_ ] := fromv/@ Propagator[t]@@ v

fromp[ _, v_, _[t_] ] := fromv/@ Propagator[t]@@ v


Attributes[psort] = Attributes[plist] = {Orderless}

plist[ t_, r___ ] := { Map[renum, t, {2}], plist[r] }

plist[ ] = {}


renum[ x:_[n_, ___] ] := x /; n < 500

renum[ x_Loop ] := x = Loop[++lc]

renum[ x:h_[n_, r___] ] := (x = h[++vc, r]) /; n > 1000

renum[ x:h_[_, r___] ] := x = h[--pc, r]



xlat[ _[vert[e_][nold_, 0]], vert[e_][nnew_, 0] ] :=
  Vertex[e][nold - 1000] -> Vertex[e][nnew]

xlat[ _[vert[e_][nold_, c_]], vert[e_][nnew_, c_] ] :=
  Vertex[e, c][nold - 1000] -> Vertex[e][nnew]

xlat[ _[Loop[nold_]], Loop[nnew_] ] :=
  Loop[nold - 1000] -> Loop[nnew]


buildsmap[ p_[a_, b_] ] := {-1, p[b, a]} /; !OrderedQ[{a, b}]

buildsmap[ p_ ] := {1, p}


buildomap[a___, {b_, c___}, b_, n_] := {a, n, c}

buildomap[a___, {b_, c___}, x_, n_] := buildomap[a, b, {c}, x, n]


CanonicalOrder[ top_ ] :=
Block[ {ntop, otop, xtab, vmap, smap, omap, cc = 0},
  Block[ {vert, Loop, lc = 0, vc = 0, pc = 0},
    ntop = Apply[fromp, Flatten[plist@@ Map[tosort, top, 2]], 1];
    xtab = Join[ SubValues[vert], DownValues[Loop] ];
  ];
  vmap = Apply[xlat, xtab, 1];
  otop = List@@ top /. vmap;
  {smap, otop} = Transpose[buildsmap/@ otop];
  omap = Fold[buildomap[#1, #2, ++cc] &, ntop, otop];
  smap = smap[[omap]];
  {Head[top]@@ ntop, omap, smap, vmap}
]

CanonicalOrder[_[h_[top_]], ginfo_] :=
Block[ {vert, Loop, lc = 0, vc = 0, pc = 0, omap, smap, vmap,
newtop, newvert, newginfo},
  {newtop, omap, smap, vmap} =
    CanonicalOrder[top /. {Incoming -> AAA, Outgoing -> AAB}];
  newtop = h[newtop /. {AAA -> Incoming, AAB -> Outgoing}];
  If[ dup[newtop], Return[Duplicate] ];
  dup[newtop] = True;
  newvert = ginfo[[1]] /. vmap;
  newginfo = MapThread[ UpdateGraphInfo,
    { Apply[List, newtop[[1]], {0, 1}] /. newvert,
      ginfo[[2, omap]],
      ginfo[[3, omap]],
      smap } ];
  { newtop, Prepend[Transpose[newginfo], newvert] }
]

DefaultRadius = 1.3

NPi = N[Pi]

FlipProp[xy_List, _] := xy

FlipProp[height_, sign_] := sign height

	(* if a straight prop is reversed, we have
	   to flip the label to the other side *)
FlipLabel[0 | 0., sign_] := sign

FlipLabel[__] = 1


UpdateGraphInfo[ x_, height_, 0 | 0., sign_ ] :=
  {FlipProp[height, sign], FlipLabel[height, sign] 1.}

UpdateGraphInfo[ _, height_, {dr_, 0 | 0.}, sign_ ] :=
  {FlipProp[height, sign], FlipLabel[height, sign] dr/DefaultRadius}

UpdateGraphInfo[ {from_, to_}, height_, {dr_, dphi_}, sign_ ] :=
Block[ {cs, ctr, mid, dir, circ, pos, alpha, angle, rad, flip},
  cs := cs = {Cos[ommc], Sin[ommc]};
  If[ to === from,
    mid = If[Head[height] === List, height, from + {0., 2.}];
    ctr = .5 (from + mid);   
    rad = Distance[from, ctr];
    ommc = Orientation[from, ctr];
    dir = ommc + .5 NPi;
    alpha = NPi,
  (* else *)
    lab = .5 Distance[from, to];
    dir = Orientation[from, to];
    ommc = dir - .5 NPi;
    ctr = mid = .5 (from + to);
    If[ TrueQ[height != 0],
      If[ height < 0, ommc += NPi ];
      alpha = 2. ArcTan[h = Abs[height]];
      rad = lab / Sin[alpha];
      mid += h lab cs;
      h = If[h > 1., -1, 1],
    (* else *)
      rad = 2000.;
      alpha = ArcSin[lab / rad];
      h = 1
    ];
    ctr -= h Sqrt[rad^2 - lab^2] cs;
  ];

  pos = ctr + (rad + dr) Through[{Cos, Sin}[ommc + dphi]];
  rad = Distance[pos, mid];
  angle = Orientation[mid, pos] - ommc;
  If[Abs[angle] < .1, angle = 0];
  If[Abs[angle - NPi] < .1, angle = NPi];
  flip = FlipLabel[height, sign];
  { FlipProp[height, sign],
    Which[
      angle == 0, flip rad/DefaultRadius,
      angle == NPi, -flip rad/DefaultRadius,
      True, {rad, If[flip < 0, NPi + angle, angle]} ] }
]


Distance[ p1_, p2_ ] := Block[ {d = p2 - p1}, Sqrt[d . d] ]
    
Orientation[ p1_, p2_ ] := N[ArcTan@@ (p2 - p1)]



Convert[file_] :=
Block[ {fi = StringReplace[file, "backup/" -> ""]},
  Convert[ file, fi, ToExpression[StringReplace[fi, ".m" -> ""]] ]
]

Convert[infile_, outfile_, tag_] :=
Block[ {dv, dup, n},
  Print["converting ", outfile];
  Get[infile];
  dv = DownValues[tag];
  Clear[tag];
  dv = Apply[CanonicalOrder, dv, 1];
  Switch[ n = Count[dv, Duplicate],
    0, Null,
    1, Print["  found 1 duplicate GraphInfo"],
    _, Print["  found ", n, " duplicate GraphInfos"]
  ];
  Apply[Set, DeleteCases[dv, Duplicate], 1];
  Put[ Definition[tag], outfile ];
  Clear[tag];
  outfile
]



If[ FileType["GraphInfo"] === Directory,

  SetDirectory["GraphInfo"];
  If[FileType["backup"] === None,
    Print["backing up GraphInfo/*.m in GraphInfo/backup/*.m"];
    Check[
      (CreateDirectory["backup"];
       RenameFile[#, "backup/" <> #]&/@ FileNames["Inc*.m"]),
      Abort[] ]
  ];

  Check[ Convert/@ FileNames["Inc*.m", "backup"], Abort[] ]

]

