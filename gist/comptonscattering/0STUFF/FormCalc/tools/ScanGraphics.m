(*
	ScanGraphics.m
		special graphing routines for parameter scans
		last modified 7 Oct 02 th
*)


BeginPackage["ScanGraphics`"]

ScanPlot3D::usage = "ScanPlot3D[var1, var2, n, opts] makes a 3D plot of a
scan performed in the variables var1 and var2. A range may be given with
the variable name, as in {var1, min1, max1}. The data to be plotted must
be stored in data sets (e.g. as by ReadData) and n specifies which data
sets are used (5 means 1...5, or {7, 10} means 7...10). The argument opts
may contain additional Graphics3D options."

ScanDensityPlot::usage = "ScanDensityPlot[var1, var2, n, opts] makes a
density plot of a scan performed in the variables var1 and var2. A range
may be given with the variable name, as in {var1, min1, max1}. The data
to be plotted must be stored in data sets (e.g. as by ReadData) and n
specifies which data sets are used (5 means 1...5, or {7, 10} means
7...10). The argument opts may contain additional Graphics options."

ScanContourPlot::usage = "ScanContourPlot[var1, var2, n, opts] makes a
contour plot of a scan performed in the variables var1 and var2. A range
may be given with the variable name, as in {var1, min1, max1}. The data
to be plotted must be stored in data sets (e.g. as by ReadData) and n
specifies which data sets are used (5 means 1...5, or {7, 10} means
7...10). The argument opts may contain additional Graphics options."

$MissingPoints::usage = "$MissingPoints contains a list of points which
have been found to be missing in the scan. $MissingPoints is updated by
ScanPlot3D, ScanDensityPlot, and ScanContourPlot."

NoClipping::usage = "NoClipping is an option of ScanPlot3D,
ScanDensityPlot, and ScanContourPlot. If set to true, it allows the
removal of the outermost row of data if that contains the only missing
points in the plot. This is useful to view preliminary plots of ongoing
calculations because Mathematica refuses to make proper contour plots if
there are missing values."

Assign::usage = "Assign[z, {var1, var2}, c, n] is called by ScanPlot3D,
ScanDensityPlot, and ScanContourPlot for every data set that is plotted
to build a grid of plot points. Assign is supposed to define a value
for z[v1, v2], where v1 and v2 are the numerical values of var1 and
var2 as given by the nth parameter set. An optional column index c may
be given in ScanPlot3D, ScanDensityPlot, or ScanContourPlot.  Its
default value is -1, i.e. by default the assignment
z[v1, v2] = Data[n][[1, -1]] is made. Assign can of course be redefined
to perform more complex tasks, like finding the minimum of all sets
that have the same values of var1 and var2."

z::usage = "z[x, y] represents a missing point at coordinates x and y."

$RoundFactor::usage = "$RoundFactor is a factor with which each data
point is multiplied and then rounded to an integer. This is done to
remove spurious rounding effects which can prevent data points to be
found when building the grid."


Begin["`Private`"]

$RoundFactor = 1000

Assign::novar = "Variable `` not contained in Para[``]."

Assign[z_, {var1_, var2_}, c_, n_] :=
Block[ {v1, v2},
  {v1, v2} = {var1, var2} /. Global`Para[n];
  If[ v1 === var1, Message[Assign::novar, var1, n]; Abort[] ];
  If[ v2 === var2, Message[Assign::novar, var2, n]; Abort[] ];

  (* make sure we have no spurious effects: *)
  {v1, v2} = Round[$RoundFactor {v1, v2}]/$RoundFactor;

  z[v1, v2] = Global`Data[n][[1, c]]
]


FindRange[l_List, min_, max_] :=
Block[ {dif, c, delta, smin, smax, lo, hi},
  dif = Apply[#2 - #1 &, Partition[l, 2, 1], 1];
  c = Count[dif, #]&/@ dif;
  delta = dif[[ Position[c, Max[c], 1, 1][[1, 1]] ]];

  lo = Min[l];
  If[ lo < min, lo += delta Ceiling[(min - lo)/delta] ];
  hi = Max[l];
  If[ hi > max, hi -= delta Floor[(hi - max)/delta] ];

  {lo, hi, delta}
]


Attributes[Clip] = {HoldFirst}

Clip[range_, {miss_}] := (range[[2]] -= range[[3]]; True) /;
  range[[2]] === miss

Clip[__] = False


ToLabel[{var_, __}] := ToString[var]

ToLabel[var_] := ToString[var]


$MissingPoints::missing =
"Warning: there are missing points in your data set. See $MissingPoints."

$MissingPoints::clipped =
"Warning: edge of plotting region clipped (calculation probably still in
progress)."

MakeGrid[{var1_, min1_, max1_}, {var2_, min2_, max2_},
  n1_Integer | {n1_Integer, n2_Integer}, gtype_, c_Integer:-1, gopt___] :=
Block[ {i, xrange, yrange, z = z, grid, miss, graph},
  Do[Assign[z, {var1, var2}, c, i], {i, n1, n2}];
  xrange = FindRange[ Union[#[[1, 1, 1]]&/@ DownValues[z]], min1, max1 ];
  yrange = FindRange[ Union[#[[1, 1, 2]]&/@ DownValues[z]], min2, max2 ];
  grid = Outer[z, Range@@ xrange, Range@@ yrange];

  $MissingPoints = Cases[grid, _z, {2}];
  If[ Length[$MissingPoints] =!= 0,
    If[ (NoClipping /. {gopt}) =!= True &&
        (Clip[xrange, Union[#[[1]]&/@ $MissingPoints]] ||
         Clip[yrange, Union[#[[2]]&/@ $MissingPoints]]),
      grid = Outer[z, Range@@ xrange, Range@@ yrange];
      Message[$MissingPoints::clipped],
    (* else *)
      Message[$MissingPoints::missing] ];
    $MissingPoints = Map[N, $MissingPoints, {2}]
  ];

  graph = gtype[Transpose[grid],
    MeshRange -> {Drop[xrange, -1], Drop[yrange, -1]}, gopt];
  graph = DeleteCases[graph, NoClipping -> _];
  Show[graph];
  graph
]

MakeGrid[var1_, r__] := MakeGrid[{var1, -Infinity, Infinity}, r]

MakeGrid[var1_List, var2_, r__] := MakeGrid[var1, {var2, -Infinity, Infinity}, r]


ScanPlot3D[var1_, var2_, n:_Integer | {_Integer, _Integer}, gopts___] :=
  MakeGrid[var1, var2, n, SurfaceGraphics, gopts,
    Axes -> True,
    AxesLabel -> {ToLabel[var1], ToLabel[var2], ""},
    PlotRange -> All
  ]

ScanDensityPlot[var1_, var2_, n:_Integer | {_Integer, _Integer}, gopts___] :=
  MakeGrid[var1, var2, n, DensityGraphics, gopts,
    Frame -> True,
    FrameLabel -> {ToLabel[var1], ToLabel[var2]},
    PlotRange -> All
  ]

ScanContourPlot[var1_, var2_, n:_Integer | {_Integer, _Integer}, gopts___] :=
  MakeGrid[var1, var2, n, ContourGraphics, gopts,
    Frame -> True,
    FrameLabel -> {ToLabel[var1], ToLabel[var2]},
    PlotRange -> All
  ]

End[]

EndPackage[]

