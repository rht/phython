(*
	Graphics.m
		Graphics routines for FeynArts
		last modified 2 Jun 12 th
*)

Begin["`Graphics`"]

NPi = N[Pi]

	(* for "PS" format: paper size in bp for finding margins *)
PaperSize = {595, 842}	(* A4, use 72 {8.5, 11} for Letter *)

DefaultImageSize = 72 {6, 7}

	(* Dimensions of a single diagram *)
DiagramBorder = 1

DiagramCanvas = 20

DiagramSize = DiagramCanvas + 2 DiagramBorder

LabelFontSize = 2 (* units of the diagram coordinate system *)

LabelFont = "Helvetica"

	(* for ordinary vertices and counter terms: *)
PropagatorThickness = Thickness[.11/DiagramSize]

VertexThickness = PointSize[8 .11/DiagramSize]

CounterThickness = Thickness[3 .11/DiagramSize]

	(* for counter terms: *)
CrossWidth = 1

ArrowLength = 1.2

	(* height of an arrow's thick end over the baseline *)
ArrowHeight = .4

DampingConst = (.65 ArrowLength)^4

ScalarDashing = Dashing[{.66, .66}/DiagramSize]

GhostDashing = Dashing[{.11, .66}/DiagramSize]

	(* # of points per unit arc length used in drawing Sine and
	   Cycles lines *)
NPoints = 10

	(* # of Sine or Cycles wave crests per unit arc length *)
NCrestsSine = .5

NCrestsCycles = .65

SineAmp = .25

CyclesAmp = .6

	(* the breadth of the "spikes" of the cycloid *)
CyclesBreadth = .16

	(* default radius for label positioning, must be the
	   same as DEFAULT_RADIUS in TopologyEditor.java *)
DefaultRadius = 1.3


Options[ Paint ] = {
  PaintLevel -> InsertionLevel,
  ColumnsXRows -> 3,
  AutoEdit -> True,
  SheetHeader -> Automatic,
  Numbering -> Full,
  FieldNumbers -> False,
  If[ $VersionNumber >= 6 && $Notebooks,
    DisplayFunction :> (Print/@ Render[#] &),
    DisplayFunction :> $DisplayFunction ]
}

Paint::nolevel =
"Warning: Level `1` is not contained in this insertion."

Paint::colrows =
"ColumnsXRows is not an integer or a pair of integers."

Paint[ top:P$Topology, opt:P$Options ] :=
  Paint[ TopologyList[top], opt ]

Paint[ top:(P$Topology -> _), opt:P$Options ] :=
  Paint[ TopologyList[][top], opt ]

Paint[ tops_TopologyList, options:P$Options ] :=
Block[ {fnum, ghead, opt = ActualOptions[Paint, options]},
  fnum = FieldNumbers /. opt;
  ghead = Switch[ ghead = SheetHeader /. opt,
    None | False,
      FeynArtsGraphics[],
    Automatic | True, 
      ghead = #[[0, 1]]&/@ tops[[1]];
      FeynArtsGraphics[
        Count[ghead, Incoming | External] -> Count[ghead, Outgoing] ],
    _,
      FeynArtsGraphics[ghead]
  ];

  PaintSheet[tops]
]

Paint[ tops:TopologyList[info___][___], options:P$Options ] :=
Block[ {plevel, ins, ghead,
fnum = False, opt = ActualOptions[Paint, options]},
  If[ (plevel = ResolveLevel[PaintLevel /. opt /. {info} /.
        Options[InsertFields]]) === $Failed,
    Return[$Failed] ];
  ins = PickLevel[plevel][tops];
  If[ FreeQ[plevel, Generic],
    ins = Select[
      ins /. Insertions[Generic][ gr__ ] :> Join@@ TakeIns/@ {gr},
      Length[#] =!= 0 & ] ];
  Scan[
    If[FreeQ[ins, Insertions[#]], Message[Paint::nolevel, #]]&,
    plevel ];

  If[ InitializeModel[ Model /. {info} /. Options[InsertFields],
    GenericModel -> (GenericModel /. {info} /. Options[InsertFields]),
    Reinitialize -> False ] === $Failed, Return[$Failed] ];

  ghead = Switch[ ghead = SheetHeader /. opt,
    None | False,
      FeynArtsGraphics[],
    Automatic | True, 
      If[ (ghead = Process /. {info}) === Process,
        FeynArtsGraphics[],
      (* else *)
        FeynArtsGraphics[Map[TheLabel[#, External]&, ghead, {2}] /.
          _Index -> Null]
      ],
    _,
      FeynArtsGraphics[ghead]
  ];

  PaintSheet[ ins //.
    (x:FeynmanGraph[___][__] -> Insertions[_][gr___]) :> Seq[x, gr] ]
]


PaintSheet[ tops_ ] :=
Block[ {auto, disp, cols, rows, dhead, g, topnr = 0, runnr = 0},
  auto = AutoEdit /. opt /. {False -> 2, _ -> 1};
  Switch[ cols = ColumnsXRows /. opt,
    _Integer,
      rows = cols,
    {_Integer, _Integer},
      {cols, rows} = cols,
    _,
      Message[Paint::colrows]; cols = rows = 3
  ];
  Switch[ Numbering /. opt,
    Full,
      dhead[gr_] := 
      Block[ {specs},
        ++runnr;
        specs = Cases[Head[gr], lev_ == n_ :>
          {" ", StringTake[ToString[lev], 1], ToString[n]}];
        If[ Length[specs] =!= 0,
          specs = {specs, " N", ToString[runnr]} ];
        DiagramGraphics["T" <> ToString[topnr] <> specs]
      ],
    Simple,
      dhead[_] := DiagramGraphics[ ToString[++runnr] ],
    _,
      dhead[_] = DiagramGraphics[]
  ];

  g = Flatten[TopologyGraphics/@ List@@ tops] /. _Index -> Null;
  g = Flatten[{g, Table[Null, {Mod[rows cols - Length[g], rows cols]}]}];
  g = ghead@@ Partition[Partition[g, cols], rows];
  (DisplayFunction /. opt)[g];
  g
]


TopologyGraphics[ top_ -> gr_ ] :=
Block[ {shapedata, gtop, vertexplot},
  shapedata = Shape[gtop = top /. Vertex[e_, _] -> Vertex[e], auto];
  gtop = Transpose[{
    List@@ gtop /. shapedata[[1]],
    shapedata[[2]],
    shapedata[[3]] }];
  vertexplot = VGraphics/@ Vertices[top] /. shapedata[[1]];
  ++topnr;
  FAPrint[2, "> Top. ", topnr, ": ", Pluralize[{Length[gr]}, " diagram"]];
  dhead[#]@@ Flatten[{
    PGraphics@@@ (gtop /. List@@ # /. Field[_] -> 0),
    vertexplot }]&/@ List@@ gr
]

TopologyGraphics[ top_ ] :=
  TopologyGraphics[
    MapIndexed[Append[ Take[#1, 2], #2[[1]] ]&, top] -> {{}} ] /; fnum

TopologyGraphics[ top_ ] :=
  TopologyGraphics[ Append[Take[#, 2], 0]&/@ top -> {{}} ]


VGraphics[ _[e_, c_:0][n_] ] := VertexGraphics[c][ Vertex[e][n] ]

PGraphics[ _[from_, to_, 0, ___], height_, labelpos_ ] :=
  PropagatorGraphics[Straight][from, to, height, labelpos]

PGraphics[ _[type_][from_, to_, particle_, ___], height_, labelpos_ ] :=
  PropagatorGraphics[
    PropagatorType[particle],
    TheLabel[particle, ResolveType[type]],
    PropagatorArrow[particle] ][from, to, height, labelpos]


Format[ DiagramGraphics[h___][__] ] := SequenceForm["[", h, "]"]

Format[ FeynArtsGraphics[h___][l__List] ] :=
  FeynArtsGraphics[h]@@ MatrixForm/@ {l}


SetOptions[OpenWrite, CharacterEncoding -> {}]

Unprotect[Show, Display, Export]

Show[ g:FeynArtsGraphics[___][___] ] := Show/@ Render[g]

Display[ chan_, g:FeynArtsGraphics[___][___], format___String,
  opt:P$Options ] :=
Block[ {rg},
  rg = Render[g, InferFormat[chan, format],
    ImageSize -> (ImageSize /. {opt} /. Options[Display] /.
                   ImageSize -> Automatic)];
  MapThread[Display[##, format, opt]&,
    {FilePerSheet[chan, Length[rg]], rg}]
]

Display[ chan_, s_String, ___ ] := (WriteString[chan, s]; Close[chan])

Export[ chan_, g:FeynArtsGraphics[___][___], format___String,
  opt:P$Options ] :=
Block[ {rg},
  rg = Render[g, InferFormat[chan, format],
    ImageSize -> (ImageSize /. {opt} /. Options[Export] /.
                   ImageSize -> Automatic)];
  MapThread[Export[##, format, opt]&,
    {FilePerSheet[chan, Length[rg]], rg}]
]

Export[ chan_, s_String, ___ ] := (WriteString[chan, s]; Close[chan])

	(* we have to play some tricks to get our definitions
	   for Export before Mma's ones *)
DownValues[Export] =
  Sort[DownValues[Export], !FreeQ[#1, chan] && FreeQ[#2, chan] &]

Protect[Show, Display, Export]


InferFormat[ _, format_ ] = format

InferFormat[ file_String ] := "PS" /; StringMatchQ[file, "*.ps"]

InferFormat[ file_String ] := "EPS" /; StringMatchQ[file, "*.eps"]

InferFormat[ file_String ] := "TeX" /; StringMatchQ[file, "*.tex"]

InferFormat[ _ ] = Sequence[]


Options[Render] = {ImageSize -> Automatic}

getsize[ opt___, def_ ] :=
  If[NumberQ[#], {#, #}, #]& @
    Round[ImageSize /. ActualOptions[Render, opt] /. Automatic -> def]

prologue := prologue =
  ReadList[ToFileName[$FeynArtsProgramDir, "FeynArts.pro"],
    Record, RecordSeparators -> ""][[1]]

epsf = ""

Render[ g:FeynArtsGraphics[___][___], "EPS", opt___Rule ] :=
Block[ {PaperSize = imgsize, epsf = " EPSF-3.0"},
  Flatten[ Render[Head[g][#], "PS", opt]&/@ List@@ g ]
]

Render[ g:FeynArtsGraphics[___][___], "PS", opt___Rule ] :=
Block[ {imgsize = getsize[opt, DefaultImageSize], bbox,
None = 0, Forward = 1, Backward = -1},
  bbox = Round[.5 (PaperSize - imgsize)];
  bbox = {bbox, bbox + imgsize};
  { PSString[ "\
%!PS-Adobe-3.0", epsf, "\n\
%%BoundingBox: 0 0 ", PaperSize, "\n\
%%Pages: ", Length[g], "\n",
    prologue,
    MapIndexed[
      { "\n%%Page: ", #2, #2, "\ngsave\n",
        bbox, #1, "\ngrestore\nshowpage\n" }&,
      PSRender[Head[g]/@ List@@ g] ], "\n\
%%Trailer\n\
end\n\
%%EOF\n" ] }
]

Render[ g:FeynArtsGraphics[___][___], "TeX", opt___Rule ] :=
Block[ {imgsize = getsize[opt, DefaultImageSize],
None = 0, Forward = 1, Backward = -1},
  { "\\unitlength=1bp%\n\n" <> TeXRender[Head[g]/@ List@@ g] }
]

Render[ g:FeynArtsGraphics[___][___], ___String, opt___Rule ] :=
Block[ {imgsize = getsize[opt, {288, 288}],
None = 0, Forward = 1, Backward = -1},
  MmaRender[Head[g]/@ List@@ g]
]


FilePerSheet[ file_, 1 ] := {file}

FilePerSheet[ file:{__}, n_ ] := Thread[FilePerSheet[#, n]&/@ file]

FilePerSheet[ file_String, n_ ] :=
Block[ {p, pre, post},
  p = StringPosition[file, "."];
  If[ Length[p] === 0, pre = file <> "_"; post = "",
    p = p[[-1, 1]] - 1;
    pre = StringTake[file, p] <> "_";
    post = StringDrop[file, p] ];
  Array[pre <> ToString[#] <> post &, n]
] /; file =!= "stdout"

FilePerSheet[ file_, n_ ] := Table[file, {n}]


PSString[ s_String ] = s

PSString[ RGBColor[r_, g_, b_, ___] ] :=
  ToString[r] <> " " <> ToString[g] <> " " <> ToString[b] <>
  " setrgbcolor "

PSString[ CMYKColor[c_, m_, y_, k_, ___] ] :=
  ToString[c] <> " " <> ToString[m] <> " " <> ToString[y] <>
  " " <> ToString[k] <> " setcmykcolor "

PSString[ Hue[h_, s_:1, b_:1, ___] ] :=
  ToString[h] <> " " <> ToString[s] <> " " <> ToString[b] <>
  " sethsbcolor "

PSString[ GrayLevel[g_, ___] ] := ToString[g] <> " setgray "

PSString[ Thickness[t_] ] := ToString[t DiagramSize] <> " setlinewidth "

PSString[ s_ ] := ToString[s] <> " " /; Head[s] =!= List

PSString[ s__ ] := StringJoin[PSString/@ Flatten[{s}]]


TeXString[ {x_, y_} ] := "(" <> ToString[x] <> "," <> ToString[y] <> ")"

TeXString[ s_ ] := ToString[s]

TeXString[ s__ ] := StringJoin[TeXString/@ {s}]


Attributes[ PSRender ] = Attributes[ TeXRender ] =
  Attributes[ MmaRender ] = {Listable}


PSRender[ FeynArtsGraphics[h___][sheet_] ] :=
Block[ {rows, cols, g},
  {rows, cols} = Dimensions[sheet];
  g = PSRender[{Title[h], sheet}];
  PSString[cols, rows, "Layout\n" <> g]
]

Attributes[ TeXJoin ] = {Flat, OneIdentity}

TeXJoin[ n_Integer ] := ToString[n]

TeXJoin[ s1_String, s2_String ] := s1 <> "\\quad " <> s2

TeXRender[ FeynArtsGraphics[in_ -> out_][sheet_] ] :=
  TeXRender[
    FeynArtsGraphics[TeXJoin@@ Flatten[{in, "\\to", out}]][sheet] ]

TeXRender[ FeynArtsGraphics[h___][sheet_] ] :=
Block[ {rows, cols, g},
  {rows, cols} = Dimensions[sheet];
  g = TeXRender[{Title[h], sheet}];
  TeXString["\\begin{feynartspicture}", imgsize, {cols, rows}, "\n" <>
    g <> "\\end{feynartspicture}\n\n"]
]

MmaRender[ FeynArtsGraphics[h___][sheet_] ] :=
Block[ {rows, cols, fsize, g, title,
(* magnify the labels a bit for screen viewing: *)
LabelFontSize = 1.26 LabelFontSize},
  {rows, cols} = Dimensions[sheet];
  g = MapIndexed[DiagramBox, sheet, {2}];
  title = MmaRender[Title[h]];
  fsize = LabelFontSize Min[imgsize/{cols, rows}]/DiagramSize;
  Graphics[ Flatten[{g, title}],
    PlotRange -> {{0, cols}, {0, rows}} DiagramSize,
    AspectRatio -> rows/cols]
]


DiagramBox[ Null, _ ] = {}

DiagramBox[ g_, {yoff_, xoff_} ] := Inset[ MmaRender[g],
  {xoff - 1, rows - yoff} DiagramSize, {0, 0}, {1, 1} DiagramSize ]

If[ $VersionNumber < 6,
  Inset[obj_, pos_, _, size_] := Rectangle[pos, pos + size, obj]
]


Title[ ] = {}

Title[ t_ ] := (
  rows += .3;
  LabelText[t, {.5 cols, rows - .12} DiagramSize, {0, 0}, 1.2, 0] )


PSRender[ DiagramGraphics[h___][pv__] ] :=
  "\n(" <> h <> ") Diagram\n" <> Transpose[PSRender[{pv}]]

PSRender[ Null ] = "\n() Diagram\n"

TeXRender[ DiagramGraphics[h___][pv__] ] :=
  "\n\\FADiagram{" <> h <> "}\n" <> Transpose[TeXRender[{pv}]]

TeXRender[ Null ] = "\n\\FADiagram{}\n"

MmaRender[ DiagramGraphics[h___][pv__] ] :=
Block[ {g = Transpose[MmaRender[{pv}]]},
  Graphics[
    Flatten[{
      PropagatorThickness, g[[1]],
      scope[ VertexThickness, g[[2]] ],
      DiagLabel[h]
    }] /. scope[a__] :> scope@@ Flatten[{a}] /. scope -> List,
    PlotRange -> {{0, DiagramSize}, {0, DiagramSize}} - DiagramBorder,
    AspectRatio -> 1 ]
]

DiagLabel[ ] = {}

DiagLabel[ t_ ] :=
  MmaRender[ LabelText[t, {.5 DiagramCanvas, -.5}, {0, -1}, .8] ]


PSRender[ VertexGraphics[cto_][xy_] ] :=
  {{}, PSString[xy, cto, "Vert\n"]}

TeXRender[ VertexGraphics[cto_][xy_] ] :=
  {{}, TeXString["\\FAVert", xy, "{", cto, "}\n"]}

MmaRender[ VertexGraphics[0][xy_] ] := {{}, Point[xy]}

MmaRender[ VertexGraphics[c_?Negative][xy_] ] := {{},
  { scope[ GrayLevel[Max[0, 1.1 + .3 c]], Disk[xy, CrossWidth] ],
    Circle[xy, CrossWidth] }}

MmaRender[ VertexGraphics[c_][xy_] ] := {{},
  { scope[ CounterThickness,
      Line[{xy - .5 CrossWidth, xy + .5 CrossWidth}],
      Line[{xy - {.5, -.5} CrossWidth, xy + {.5, -.5} CrossWidth}] ],
    Array[Circle[xy, (.25 # + .8) CrossWidth]&, c - 1] }}


PSRender[ PropagatorGraphics[type_, label_:0, arrow_:0][
  from_, to_, height_, labelpos_:0 ] ] :=
Block[ {dir, ommc, cs, ctr, rad, mid, dphi, line},
  line = PSString["{ ", type, "} ", arrow, height, from, to, "Prop\n"];

  If[ label =!= 0,
    CalcPropData[from, to, height];
    line = line <> PSRender[PropLabel[label, labelpos, arrow, type]] ];

  {line, {}}
]

TeXRender[ PropagatorGraphics[type_, label_:0, arrow_:0][
  from_, to_, height_, labelpos_:0 ] ] :=
Block[ {dir, ommc, cs, ctr, rad, mid, dphi, line},
  line = TeXString["\\FAProp", from, to,
    If[ NumberQ[height], {height, ""}, height],
    "{" <> StringDrop[PSString[type], -1] <> "}{", arrow, "}\n"];

  If[ label =!= 0,
    CalcPropData[from, to, height];
    line = line <> TeXRender[PropLabel[label, labelpos, arrow, type]] ];

  {line, {}}
]

MmaRender[ PropagatorGraphics[type_, label_:0, arrow_:0][
  from_, to_, height_, labelpos_:0 ] ] :=
Block[ {dir, ommc, cs, ctr, rad, mid, dphi, line, phi, damping, t, h, v},
  CalcPropData[from, to, height];
  If[ arrow =!= 0,
    damping[phi_] := (#/(# + DampingConst))&[(rad Abs[phi - ommc])^4],
  (* else *)
    _damping = 1 ];
  phi = Mod[ommc - dphi, 2 NPi];
  t = Flatten[{type}];
  h = Position[t, Straight | ScalarDash | GhostDash | Sine | Cycles, 1];
  dphi *= 2./Length[h];
  line = scope[MapAt[HalfLine[#, phi += dphi, dphi]&, t, h]];

  If[ arrow =!= 0,
    h = .5 arrow ArrowLength {Cos[dir], Sin[dir]};
    v = ArrowHeight cs;
    line = {line, Polygon[{mid + h, mid - h + v, mid - h - v}]} ];

  If[ label =!= 0,
    line = {line, MmaRender[PropLabel[label, labelpos, arrow, type]]} ];

  {line, {}}
]


(* CalcPropData computes the data necessary for actual drawing:
   dir  -- the direction along the propagator
   ommc -- the direction of the perpendicular bisector
   ctr  -- the center of the circle on which the prop lies,
   rad  -- the radius of the circle
   mid  -- the position on the middle of the prop
           (the blue square in the topology editor)
   dphi -- half the opening angle *)

CalcPropData[ from_, _, xy_List ] := (
  mid = xy;
  ctr = .5 (from + mid);
  rad = Distance[from, ctr];
  ommc = Orientation[from, ctr];
  dir = ommc + .5 NPi;
  cs = {Cos[ommc], Sin[ommc]};
  dphi = NPi
)

CalcPropData[ from_, to_, height_ ] :=
Block[ {lab, h},
  lab = .5 Distance[from, to];
  dir = Orientation[from, to];
  ommc = dir + If[height < 0, .5, -.5] NPi;
  cs = {Cos[ommc], Sin[ommc]};
  ctr = mid = .5 (from + to);
  If[ height != 0,
    dphi = 2. ArcTan[h = Abs[height]];
    rad = lab/Sin[dphi];
    dphi *= Sign[height];
    mid += h lab cs;
    h = If[h > 1., -1, 1],
  (* else *)
    rad = 20000.;
    dphi = ArcSin[lab/rad];
    h = 1
  ];
  ctr -= h Sqrt[rad^2 - lab^2] cs
]


HalfLine[ Sine, phi_, dphi_ ] :=
Block[ {arc, w, n},
  arc = rad Abs[dphi];
  w = 2. NPi (.5 + Max[1, Round[NCrestsSine arc]]);
  Line[
    Table[ arc = phi - n dphi;
      ctr + (rad - damping[arc] SineAmp Sin[n w]) {Cos[arc], Sin[arc]},
      {n, 0, 1, 1./Floor[NPoints arc]} ] ]
]

rshift = CyclesAmp - SineAmp

phadj = ArcCos[rshift/CyclesAmp]

sphadj = Sin[phadj]

HalfLine[ Cycles, phi_, dphi_ ] :=
Block[ {arc, w, n, phamp},
  arc = rad Abs[dphi];
  w = 2. (phadj + NPi Max[1, Round[NCrestsCycles arc]]);
  phamp = CyclesBreadth NPi Sign[dphi]/rad;
  Line[
    Table[ arc = n w - phadj;
      ctr + (rad + CyclesAmp Cos[arc] - rshift) *
        Through[{Cos, Sin}[
          phi - n dphi - phamp (Sin[arc] - (2 n - 1) sphadj)]],
      {n, 0, 1, 1./Floor[2 NPoints arc]} ] ]
]

HalfLine[ Straight, phi_, dphi_ ] :=
  If[ rad < 20000,
    Circle[ctr, rad, Sort[{phi - dphi, phi}]],
    Line[{ ctr + rad {Cos[phi - dphi], Sin[phi - dphi]},
           ctr + rad {Cos[phi], Sin[phi]} }] ]

HalfLine[ ScalarDash, phi_, dphi_ ] :=
  scope[ ScalarDashing, HalfLine[Straight, phi, dphi] ]

HalfLine[ GhostDash, phi_, dphi_ ] :=
  scope[ GhostDashing, HalfLine[Straight, phi, dphi] ]


PropLabel[ label_, labelpos_, arrow_, type_ ] :=
Block[ {rad, phi, s},
  {rad, phi} =
    If[ NumberQ[labelpos], {labelpos DefaultRadius, 0.}, labelpos ];
  s = Sign[rad Cos[phi]];
  cs *= s;
  Which[
    arrow =!= 0 || !FreeQ[type, Sine],
      mid += SineAmp cs,
    !FreeQ[type, Cycles],
      mid += If[s < 0, CyclesAmp + rshift, SineAmp] cs ];
  phi += ommc;
  cs = -Round[1.3 cs];
  LabelText[ label,
    mid + rad {Cos[phi], Sin[phi]} + .24 LabelFontSize cs,
    cs, 1, dir ]
]


LabelText[ in_ -> out_, pos_, align_, size_, dir_ ] :=
Block[ {t, l, cs = {Cos[dir], Sin[dir]}},
  t = Flatten[{in, "\\to", out}];
  l = .5 (Length[t] + 1);
  MapIndexed[LabelText[#1, pos + 4.5 (#2[[1]] - l) cs, align, size]&, t]
]

LabelText[ t_List, pos_, align_, size_, dir_ ] :=
Block[ {l = .5 (Length[t] + 1), cs = {Cos[dir], Sin[dir]}},
  MapIndexed[LabelText[#1, pos + 2 (#2[[1]] - l) cs, align, size]&, t]
]

PSRender[ LabelText[t_, pos_, align_, size_, ___] ] :=
  PSString["{ ", MapIndexed[PSChar, Flatten[{ToPS[t]}]], "} ",
    pos, align, size, "Label\n"]

PSChar[ _[], _ ] = {}

PSChar[ _[c_], {n_} ] := {"$(", c, psops[[n]]}

PSChar[ c_, {n_} ] := {"(", c, psops[[n]]}

psops = {")# ", ")_ ", ")^ ", ")~ "}


TeXRender[ LabelText[t_, pos_, align_, size_, ___] ] :=
Block[ {ComposedChar = TeXComposedChar},
  TeXString[
    "\\FALabel", pos, "[" <> Extract[texalign, align + {2, 2}] <> "]{" <>
    Which[size > 1, "\\large ", size < 1, "\\small ", True, ""] <>
    "$", t, "$}\n" ]
]

texalign = {{"bl", "l", "tl"},
            {"b",  "",  "t"},
            {"br", "r", "tr"}}

TeXComposedChar[ t_, sub_:Null, super_:Null, over_:Null ] :=
Block[ {tex = t},
  If[ sub =!= Null, tex = tex <> "_" <> ToString[sub] ];
  If[ super =!= Null, tex = tex <> "^" <> ToString[super] ];
  If[ over =!= Null, tex = ToString[over] <> " " <> tex ];
  tex
]

If[ $Notebooks,

MmaRender[ LabelText[t_, r__] ] :=
  NotebookChar[Flatten[{ToUnicode[t]}], r];

NotebookChar[ {t_, sub_:" ", super_:" ", over_:" "},
  pos_, align_, size_, ___ ] :=
Block[ {label = t},
  If[ over =!= " ", label = OverscriptBox[label, over] ];
  Which[
    sub =!= " " && super =!= " ",
      label = SubsuperscriptBox[label, sub, super],
    sub =!= " ",
      label = SubscriptBox[label, sub],
    super =!= " ",
      label = SuperscriptBox[label, super]
  ];
  Text[
    StyleForm[DisplayForm[label], FontFamily -> LabelFont,
      FontSize -> size fsize],
    pos, align ]
],

(* else $Notebooks *)

MmaRender[ LabelText[t_, r__] ] :=
  MapIndexed[ KernelChar[##, r]&, Flatten[{ToPS[t]}] ];

KernelChar[ _[], __ ] = {};

KernelChar[ t_, {n_}, pos_, align_, size_, ___ ] :=
Block[ {newpos, newalign, fscale = sizes[[n]] size},
  {newpos, newalign} = If[ size < 1, {pos, align},
    {pos + .24 size LabelFontSize (palign[[n]] - align), talign[[n]]} ];
  Text[ MmaChar[t], newpos, newalign ]
];

MmaChar[ _[c_] ] := FontForm[c, {"Symbol", fscale fsize}];

MmaChar[ c_ ] := FontForm[c, {LabelFont, fscale fsize}];

sizes = {1, .667, .667, 1};

palign = {{0, 0}, {1.2, -.9}, {1.2, .7}, {0, 1.2}};

talign = {{0, 0}, {-1, 0}, {-1, 0}, {0, -1}};

] (* endif $Notebooks *)


FindFlip[ h_, top_, rulz_ ] :=
Block[ {ntop = top /. rulz, ord, shapedata},
  ord = TopologyOrdering[ntop];
  FlipShape[h, ntop, ord,
    GetShape[TopologyCode[ ord[[1]] ]] /. Reverse/@ rulz]
]


FlipShape[ h_, top_, _[ntop_, map_], {vert_, prop_, labels_} ] :=
Block[ {ord = Ordering[Abs[map]], vord, vmap},
  vord = Flatten[(4 Abs[#] - 1 + Sign[#] {-1, 1})/2&/@ map];
  vmap = Thread[ Level[ntop, {2}] -> Level[top, {2}][[vord]] ];
  Throw[
    FAPrint[2, "  found shape through ", h, " flip"];
    { Sort[h[1]/@ vert /. vmap],
      Sequence@@ Transpose[MapThread[h[23],
        {prop[[ord]], labels[[ord]], map}]] } ]
]


TopBottom[1][ v_ -> {x_, y_} ] := v -> {x, DiagramCanvas - y}

TopBottom[23][ {x_, y_}, {r_, fi_}, _ ] := {{x, DiagramCanvas - y}, {r, -fi}}

TopBottom[23][ {x_, y_}, l_, _ ] := {{x, DiagramCanvas - y}, l}

TopBottom[23][ h_, {r_, fi_}, s_ ] := {-Sign[s] h, {r, -Sign[s] fi}} /; h != 0

TopBottom[23][ h_, l_, s_ ] := {-Sign[s] h, l} /; h != 0

TopBottom[23][ _, {r_, fi_}, s_ ] := {0, {r, If[s > 0, Pi - fi, -fi]}}

TopBottom[23][ _, l_, s_ ] := {0, -Sign[s] l}


LeftRight[1][ v_ -> {x_, y_} ] := v -> {DiagramCanvas - x, y}

LeftRight[23][ {x_, y_}, {r_, fi_}, _ ] := {{DiagramCanvas - x, y}, {r, -fi}}

LeftRight[23][ {x_, y_}, l_, _ ] := {{DiagramCanvas - x, y}, l}

LeftRight[23][ h_, {r_, fi_}, s_ ] := {-Sign[s] h, {r, -fi}} /; h != 0

LeftRight[23][ h_, l_, s_ ] := {-Sign[s] h, l} /; h != 0

LeftRight[23][ _, {r_, fi_}, s_ ] := {0, {r, If[s > 0, Pi - fi, -fi]}}

LeftRight[23][ _, l_, s_ ] := {0, -Sign[s] l}


AutoShape[ top_ ] :=
Block[ {in, out, vert, shapedata, props, ext, l, tree, mesh, mesh2,
vars, tadbr, tad, min, ok, c, ct, pt, shrink = {}, rev = {}, loops = {}},

	(* before embarking on any serious autoshaping, check if we know
	   how to paint the vertical or horizontal mirror image of top *)
  in = Sort[Cases[top, _[Incoming][v_, _] :> v]];
  out = Sort[Cases[top, _[Outgoing][v_, _] :> v]];
  vert = Join[out, in];

  FAPrint[1, "shaping topology ", TopologyCode[top]];

  FindFlip[ TopBottom, top, Thread[vert -> Reverse[Join[in, out]]] ];

  FindFlip[ LeftRight, top, Flatten[{Thread[vert -> Sort[vert]],
    Incoming -> Outgoing, Outgoing -> Incoming}] ];

  FAPrint[2, "  autoshaping"];

  props[_] = {};
  top /. Propagator -> pr;
  loops = Union[loops];

  Off[FindMinimum::fmmp, FindMinimum::fmcv, FindMinimum::precw,
    FindMinimum::fmgz, FindMinimum::sdprec, FindMinimum::lstol];

	(* a) fix the incoming and outgoing propagators on the left and
	      right side, respectively *)
  shapedata = Join[
    l = Length[props[Incoming]]/20.;
    MapIndexed[
      #1[[1]] -> {0, 20 - Round[(#2[[1]] - .5)/l]}&,
      props[Incoming] ],
    l = Length[props[Outgoing]]/20.;
    MapIndexed[
      #1[[1]] -> {20, 20 - Round[(#2[[1]] - .5)/l]}&,
      props[Outgoing] ] ];

	(* b) shrink loops to 1 point which is the center of an imaginary
	      circle on which the external points of the loop lie *)
  vert = Flatten[props[Tree]];
  tree = Fold[
    ( l = Flatten[props[#2]];
      ext[#2] = l = Select[vert, MemberQ[l, #]&];
      AppendTo[ rev, c = center[Length[l]][#2] -> Union[l] ];
      shrink = {shrink, c = Thread[Reverse[c]]};
      #1 /. c )&,
    props[Tree], loops ];
  shrink = Flatten[shrink];

	(* c) cut tadpole-like parts and minimize the length of the
	      remaining mesh of propagators *)
  mesh2 = Leaves@@ twig@@@ (tree /. shapedata) /. twig -> List;
  _tadbr = {};
  Cases[mesh2, branch[ctr:center[_][_], v_, ___] :>
    (tadbr[ctr] = Flatten[{tadbr[ctr], v /. rev}]), Infinity];
  mesh = mesh2 /. branch[__] :> Seq[];
  vert = Cases[mesh, leaf[a_] -> a];
  mesh = mesh /. _leaf :> Seq[];
  mesh = List@@ Fold[
    Replace[#1, x:{__, #2} :> Reverse[x], {1}] /.
      Leaves[x:{#2, __}, {#2, b__}] :> Leaves[{Sequence@@ Reverse[x], b}]&,
    mesh,
    vert = Union[vert, Cases[mesh, _[2, ___][_], {2}]]
  ] /. twig -> List;
  vert = (# -> CartesianVar[#])&/@
    Complement[Cases[mesh, _[__][_], {2}], vert];
  If[ Length[vert] =!= 0,
    dist = Distance2@@@ (mesh /. vert);
    dist = Expand[Apply[Plus, Outer[(#1 - #2)^2 &, dist, dist], {0, 1}]];
    vars = Flatten[Last/@ vert];
    min := FindMinimum@@ Prepend[{#, 10 + Random[]}&/@ vars, dist];
    Do[
      If[ ok = (Head[c = min] =!= FindMinimum &&
                Min[ l = Last/@ c[[2]] ] >= 0 &&
                Max[l] <= 20), Break[] ],
    {5} ];
    vert = vert /.
      If[ ok, MapAt[Round, #, 2]&/@ c[[2]], (# -> RandInt)&/@ vars ]
  ];
  shapedata = Join[shapedata, vert];
  shapedata = Flatten[ {shapedata,
    (l = (#[[-1]] - (c = #[[1]]))/(Length[#] - 1);
     MapIndexed[#1 -> c + l #2[[1]] &, Take[#, {2, -2}]])&/@
      (Select[mesh, Length[#] > 2 &] /. shapedata)} ];

	(* d) minimizing the straight distance of a tadpole to its
	      nearest vertex v would make the tadpole stick to v.
	      Therefore, the tadpole is given polar coordinates with
	      fixed radius and we try to construct the ideal angle with
	      respect to the lines joining at v by maximizing. *)
  mesh2 = (List@@ mesh2) /.
    leaf[br__] :> Sequence@@ Cases[{br}, branch[__]] /. shapedata;
  While[ Length[ tad = Cases[mesh2, branch[{_, _}, __]] ] =!= 0,
    mesh2 = Fold[FixTad, mesh2, tad] ];

	(* e) for each loop, distribute the external points of the loop
	      at the middle of the line from the center to the external
	      vertex and distribute the remaining points of the loop
	      on the imaginary circle around the center. *)
  ok = shapedata;
  Scan[
    Function[rul,
      If[ Length[ vert = Union[ext[ rul[[1, 1]] ]] ] === 1,
        c = SetTadpole[vert[[1]], rul],
        SetMiddle[#, rul]&/@ vert; c = rul[[2]] ];
      SetLoop[props[ rul[[1, 1]] ] /. shapedata, c] ],
    Select[shapedata, !FreeQ[#, center]&] ];

	(* f) last resort: randomize any remaining vertex *)
  shapedata = Join[ MapAt[Inside, #, 2]&/@ shapedata,
    (# -> {RandInt, RandInt})&/@
      Union[Cases[top /. shapedata, Vertex[__][_], Infinity]] ];

	(* g) give tadpoles and identical propagators curvature so that
	      they do not fall on top of each other *)
  pt[ _[Loop[n_]][v_, v_], _ ] :=
    center[ Length[ext[n]] ][n] /. shapedata /. center[_][_] :>
      Inside[ (v /. shapedata) + 4 Through[{Cos, Sin}[2. NPi Random[]]] ];
  pt[ _, 0 ] = 0;
  ct[ p_ ] := If[ (c = (Count[top, p] - 1)/2) === 0, 0,
    pt[ p, n_ ] = .8 n/c;
    ct[ p ] = c
  ];

  On[FindMinimum::fmmp, FindMinimum::fmcv, FindMinimum::precw,
    FindMinimum::fmgz, FindMinimum::sdprec, FindMinimum::lstol];

  { Select[shapedata, FreeQ[#, center]&],
    pt[#, ct[#]--]&/@ List@@ top,
    Table[1, {Length[top]}] }
]


pr[ Loop[l_] ][ from_, to_, ___ ] := (
  AppendTo[loops, l];
  AppendTo[props[l], {from, to}] )

pr[ type_ ][ from_, to_, ___ ] := (
  AppendTo[props[Tree], {from, to}];
  AppendTo[props[type], {from, to}] )


Attributes[ twig ] = {Orderless}

twig[ a:_[1, ___][_], b_ ] := branch[b, a]

Attributes[ Leaves ] = {Orderless, Flat}

Leaves[ branch[a_, b__], twig[a:_[2, ___][_], c_] ] :=
  Leaves[branch[c, a, b]]

Leaves[ br:branch[a_, __].., tw:twig[a_, _].. ] :=
  Switch[ Length[{tw}],
    1, Leaves[branch[ Sequence@@ DeleteCases[tw, a], a, br ]],
    2, Leaves[leaf[br, a], tw],
    _, Leaves[leaf[br], tw]
  ]

cutbranch[ vert__, br___branch ] :=
  Sequence[ br, Drop[{vert}, {2, -2}] ]


Inside[ xy_ ] := Max[Min[#, 20], 0]&/@ N[xy]

RandInt := Plus@@ Table[Random[Integer, 9], {2}] + 1


Distance[ p1_, p2_ ] := Block[ {d = p2 - p1}, Sqrt[d . d] ]

Distance2[ p1_, ___, p2_ ] := Block[ {d = p2 - p1}, d . d ]

Orientation[ p1_, p2_ ] := N[ArcTan@@ (p2 - p1)] /; p1 != p2

Orientation[ __ ] = 0


CartesianVar[ n_ ] := CartesianVar[n] = {Unique["X"], Unique["Y"]}


FixTad[ mesh_, br_ ] :=
Block[ {stem, root, c, phi, dphi, vert, min, dist},
  stem = DeleteCases[List@@ br, branch[__]];
  c = (root = stem[[1]]) + (Length[stem] + 1) {Cos[phi], Sin[phi]};
  vert = Cases[mesh, {a___, root, b___} :> Seq[a, b]];
  Switch[ Length[vert],
    0, phi = 0,
    1, phi = Orientation[vert[[1]], root],
    _, dist = -Plus@@ (Distance[c, #]&)/@ vert;
       Do[
         If[ Head[min = FindMinimum[dist, {phi, 1.57 + dphi}]] =!=
               FindMinimum &&
             Min[ min = c /. min[[2]] ] >= 0 &&
             Max[min] <= 20, c = min; Break[] ],
       {dphi, 0, 2. NPi, .5 NPi} ];
       phi = 2. NPi Random[]
  ];
  If[ Min[c] < 0 || Max[c] > 20, c = {RandInt, RandInt} ];
  c = (c - root)/(Length[stem] - 1);
  vert = MapIndexed[#1 -> root + #2[[1]] c &, Rest[stem]];
  shapedata = Join[shapedata, vert];
  mesh /. br -> cutbranch@@ br /. vert
]


SetLoop[ loop_, ctr_ ] :=
Block[ {vert, vars, rad, angle, off, min, c = 0},
  vert = Union[ Cases[loop, _[__][_], {2}] ];
  If[ Length[vert] =!= 0,
    rad = Union[ Cases[loop, _List, {2}] ];
    rad = If[ Length[rad] === 0, 5,
      Plus@@ (Distance[#, ctr]&/@ rad)/Length[rad] ];
    angle = 2. NPi/Length[vert];
    vars = (# ->
      ctr + rad Through[{Cos, Sin}[++c angle + off]])&/@ vert;
    min = FindMinimum@@ {
      -Plus@@ Distance@@@ (loop /. vars),
      {off, 2. NPi Random[]} };
    shapedata = Join[ shapedata,
      If[ Head[min] === FindMinimum, vars /. off -> 2. NPi Random[],
        c = vars /. min[[2]];
        If[ Length[Intersection[
          Round[Last/@ c], Round[Last/@ shapedata] ]] === 0,
          c,
          # -> ctr + Random[] rad Through[{Sin, Cos}[2. NPi Random[]]]&/@
            vert ] ]
    ]
  ];
]


SetMiddle[ vert_, ctr_ -> xy_ ] :=
Block[ {ex, mid},
  ex = DeleteCases[ Flatten[Select[props[Tree], !FreeQ[#, vert]&]], vert ];
  If[ Length[ex] =!= 0,
    If[ Length[ mid = Complement[ex, tadbr[ctr]] ] =!= 0, ex = mid ];
    ex = ex /. shrink /. shapedata;
    While[ Length[ex] =!= 0 && (mid = Plus@@ ex/Length[ex]) == xy,
      ex = Rest[ex] ];
    If[ Length[ex] === 0, mid = {RandInt, RandInt} ];
    AppendTo[shapedata, vert -> .6 xy + .4 mid]
  ]
]


SetTadpole[ vert_, ctr_ -> xy_ ] :=
Block[ {adj, max, new, a1, a2},
  adj = Select[tree, !FreeQ[#, ctr]&] /. ok;
  shapedata = shapedata /. ctr -> vert;
  If[ Length[adj] === 1,
    new = 2.6 xy - .8 Plus@@ adj[[1]],
  (* else *)
    adj = Orientation@@@ (adj /. {a_, xy} -> {xy, a});
    max = -1;
    Outer[
      If[ (new = Abs[#2 - #1]) > max, max = new; a1 = #1; a2 = #2 ]&,
      adj, adj ];
    new = xy + 4 Through[{Cos, Sin}[.5 (a1 + a2)]] ];
  If[ Distance2[xy, max = Inside[new]] < 2.8,
    shapedata = shapedata /. xy -> xy - .7 (new - max) ];
  AppendTo[shapedata, ctr -> max];
  max
]


vcode = Characters["\
abcdefghijklmnopqrstuvwxyz\
ABCDEFGHIJKLMNOPQRSTUVWXYZ\
0987654321"]

pcode[ _[Incoming | External][c__] ] := {{c}, {}, {}}

pcode[ _[Outgoing][c__] ] := {{}, {c}, {}}

pcode[ _[c__] ] := {{}, {}, {c}}

TopologyCode[ top:P$Topology ] :=
  StringJoin/@ Transpose[pcode/@ Apply[vcode[[#]]&, List@@ top, {2}]] /.
    "" -> "0"


MkDir[ dir_String ] := dir /; FileType[dir] === Directory

MkDir[ dir_String ] := Check[CreateDirectory[dir], Abort[]]

MkDir[ dirs__String ] := Fold[MkDir[ToFileName[##]]&, {}, {dirs}]


Attributes[ GetShape ] = {HoldRest}

GetShape[ {topcode__}, ___ ] :=
Block[ {shapedata = ShapeData[topcode]},
  shapedata /; Head[shapedata] === List
]

GetShape[ {ext__, int_}, ___ ] :=
Block[ {file = ToFileName[{$ShapeDataDir, ext}, int <> ".m"]},
  (ShapeData[ext, int] = Get[file]) /; FileType[file] === File
]

GetShape[ {topcode__}, alt_ ] := ShapeData[topcode] = alt


PutShape[ shapedata_, {ext__, int_} ] :=
  Put[ShapeData[ext, int] = shapedata,
    ToFileName[MkDir[$ShapeDataDir, ext], int <> ".m"]]


ToJava[ p__, n_?NumberQ ] := ToJava[p, {n DefaultRadius, 0.}]

ToJava[ _[from_, from_], {xc_, yc_}, {xl_, yl_} ] :=
  {-from, xc, yc, xl, yl}

ToJava[ _[from_, to_], height_, {xl_, yl_} ] :=
  {from, to, height, xl, yl}


(* call the topology editor *)

Shape::wait =
"Starting Java and the topology editor.  This may take a moment."

Shape::notopedit =
"Could not load the topology editor.  Make sure you have J/Link and Java
installed."

Shape::javaerror =
"Could not open a topology-editor window."

Shape[ tops:TopologyList[___][___] | TopologyList[___] ] :=
Block[ {remaining = Length[tops]},
  MapIndexed[
    (--remaining; FAPrint[ 2, "> Top. ", #2[[1]] ]; Shape[#1])&,
    List@@ tops ]
]

Shape[ top:P$Topology -> _ ] := Shape[top]

Shape[ top:P$Topology, auto_:0 ] :=
Block[ {edittop, topcode, shapedata, arg1, arg2, exitcode},
  edittop = Take[#, 2]&/@ Topology@@ top /.
    External -> Incoming /. Vertex[e_, _] -> Vertex[e];
  topcode = TopologyCode[edittop];
  res = auto;
  shapedata = GetShape[topcode, --res; Catch[AutoShape[edittop]]];
  If[ res > 0, Return[shapedata] ];

  If[ editorclass === False,
    Message[Shape::wait];
    Needs["JLink`"];
    JLink`InstallJava[];
    AppendTo[ JLink`$ExtraClassPath,
      ToFileName[$FeynArtsProgramDir, "TopologyEditor.jar"] ];
    editorclass = JLink`LoadClass["de.FeynArts.TopologyEditor"];
    If[ Head[editorclass] =!= JLink`JavaClass, Message[Shape::notopedit] ]
  ];
  If[ Head[editorclass] =!= JLink`JavaClass, Return[shapedata] ];

  If[ !JLink`JavaObjectQ[editor],
    editor = JLink`JavaNew[editorclass];
    If[ !JLink`JavaObjectQ[editor],
      Message[Shape::javaerror];
      Return[shapedata] ];
  ];

  arg1 = Last/@ shapedata[[1]];
  arg2 = MapThread[ ToJava,
    { List@@ edittop /. MapIndexed[ First/@ Rule[##]&, shapedata[[1]] ],
      shapedata[[2]],
      shapedata[[3]] } ];
  editor@putShapeData[ N[Flatten[arg1]], N[Flatten[arg2]] ];

  exitcode = JLink`DoModal[];
  If[ exitcode === 0,
    shapedata = MapAt[
      MapThread[Rule, {First/@ shapedata[[1]], #}]&,
      editor@getShapeData[], 1 ];
    PutShape[shapedata, topcode]
  ];

  If[ remaining === 0 || exitcode === 2,
    editor@closeWindow[];
    JLink`ReleaseObject[editor];
    editor = False;
    If[ exitcode === 2, Abort[] ];
  ];

  shapedata
]


editorclass = editor = False

remaining = 0


ToPS[ "\\alpha" ] = SymbolChar["a"];
ToPS[ "\\beta" ] = SymbolChar["b"];
ToPS[ "\\gamma" ] = SymbolChar["g"];
ToPS[ "\\delta" ] = SymbolChar["d"];
ToPS[ "\\epsilon" ] = ToPS[ "\\varepsilon" ] = SymbolChar["e"];
ToPS[ "\\zeta" ] = SymbolChar["z"];
ToPS[ "\\eta" ] = SymbolChar["h"];
ToPS[ "\\theta" ] = SymbolChar["q"];
ToPS[ "\\vartheta" ] = SymbolChar["J"];
ToPS[ "\\iota" ] = SymbolChar["i"];
ToPS[ "\\kappa" ] = SymbolChar["k"];
ToPS[ "\\lambda" ] = SymbolChar["l"];
ToPS[ "\\mu" ] = SymbolChar["m"];
ToPS[ "\\nu" ] = SymbolChar["n"];
ToPS[ "\\xi" ] = SymbolChar["x"];
ToPS[ "\\pi" ] = SymbolChar["p"];
ToPS[ "\\varpi" ] = SymbolChar["v"];
ToPS[ "\\rho" ] = ToPS[ "\\varrho" ] = SymbolChar["r"];
ToPS[ "\\sigma" ] = SymbolChar["s"];
ToPS[ "\\varsigma" ] = SymbolChar["V"];
ToPS[ "\\tau" ] = SymbolChar["t"];
ToPS[ "\\upsilon" ] = SymbolChar["u"];
ToPS[ "\\phi" ] = SymbolChar["f"];
ToPS[ "\\varphi" ] = SymbolChar["j"];
ToPS[ "\\chi" ] = SymbolChar["c"];
ToPS[ "\\psi" ] = SymbolChar["y"];
ToPS[ "\\omega" ] = SymbolChar["w"];
ToPS[ "\\Gamma" ] = SymbolChar["G"];
ToPS[ "\\Delta" ] = SymbolChar["D"];
ToPS[ "\\Theta" ] = SymbolChar["Q"];
ToPS[ "\\Lambda" ] = SymbolChar["L"];
ToPS[ "\\Xi" ] = SymbolChar["X"];
ToPS[ "\\Pi" ] = SymbolChar["P"];
ToPS[ "\\Sigma" ] = SymbolChar["S"];
ToPS[ "\\Upsilon" ] = SymbolChar["\241"];
ToPS[ "\\Phi" ] = SymbolChar["F"];
ToPS[ "\\Psi" ] = SymbolChar["X"];
ToPS[ "\\Omega" ] = SymbolChar["W"];
ToPS[ "\\infty" ] = SymbolChar["\245"];
ToPS[ "\\pm" ] = SymbolChar["\261"];
ToPS[ "\\partial" ] = SymbolChar["\266"];
ToPS[ "\\leq" ] = SymbolChar["\243"];
ToPS[ "\\geq" ] = SymbolChar["\263"];
ToPS[ "\\times" ] = SymbolChar["\264"];
ToPS[ "\\otimes" ] = SymbolChar["\304"];
ToPS[ "\\oplus" ] = SymbolChar["\305"];
ToPS[ "\\nabla" ] = SymbolChar["\321"];
ToPS[ "\\neq" ] = SymbolChar["\271"];
ToPS[ "\\equiv" ] = SymbolChar["\272"];
ToPS[ "\\approx" ] = SymbolChar["\273"];
ToPS[ "\\ldots" ] = SymbolChar["\274"];
ToPS[ "\\in" ] = SymbolChar["\316"];
ToPS[ "\\notin" ] = SymbolChar["\317"];
ToPS[ "\\sim" ] = SymbolChar["\176"];
ToPS[ "\\sqrt" ] = SymbolChar["\326"];
ToPS[ "\\propto" ] = SymbolChar["\265"];
ToPS[ "\\subset" ] = SymbolChar["\314"];
ToPS[ "\\supset" ] = SymbolChar["\311"];
ToPS[ "\\subseteq" ] = SymbolChar["\315"];
ToPS[ "\\supseteq" ] = SymbolChar["\312"];
ToPS[ "\\bullet" ] = SymbolChar["\267"];
ToPS[ "\\perp" ] = SymbolChar["^"];
ToPS[ "\\simeq" ] = SymbolChar["@"];
ToPS[ "\\vee" ] = SymbolChar["\332"];
ToPS[ "\\wedge" ] = SymbolChar["\331"];
ToPS[ "\\leftrightarrow" ] = SymbolChar["\253"];
ToPS[ "\\leftarrow" ] = SymbolChar["\254"];
ToPS[ "\\rightarrow" ] = ToPS[ "\\to" ] = SymbolChar["\256"];
ToPS[ "\\uparrow" ] = SymbolChar["\255"];
ToPS[ "\\downarrow" ] = SymbolChar["\257"];
ToPS[ "\\Leftarrow" ] = SymbolChar["\334"];
ToPS[ "\\Rightarrow" ] = SymbolChar["\336"];
ToPS[ "\\Uparrow" ] = SymbolChar["\335"];
ToPS[ "\\Downarrow" ] = SymbolChar["\337"];
ToPS[ "\\bar" ] = "\305";
ToPS[ "\\hat" ] = "\303";
ToPS[ "\\tilde" ] = "\304";
ToPS[ "\\dot" ] = "\307";
ToPS[ "\\ddot" ] = "\310";
ToPS[ "\\vec" ] = SymbolChar["\256"];
ToPS[ "\\prime" ] = "'";
ToPS[ "\\#" ] = "#";
ToPS[ "\\&" ] = "&";
ToPS[ "\\$" ] = "$";
ToPS[ "\\%" ] = "%";
ToPS[ "\\_" ] = "_";
ToPS[ "-" ] = SymbolChar["-"];
ToPS[ Null ] = SymbolChar[];
ToPS[ ComposedChar[t__] ] := ToPS/@ {t};
ToPS[ c_ ] := ToString[c]

ToUnicode[ "\\alpha" ] = "\[Alpha]";
ToUnicode[ "\\beta" ] = "\[Beta]";
ToUnicode[ "\\gamma" ] = "\[Gamma]";
ToUnicode[ "\\delta" ] = "\[Delta]";
ToUnicode[ "\\epsilon" ] = "\[Epsilon]";
ToUnicode[ "\\varepsilon" ] = "\[CurlyEpsilon]";
ToUnicode[ "\\zeta" ] = "\[Zeta]";
ToUnicode[ "\\eta" ] = "\[Eta]";
ToUnicode[ "\\theta" ] = "\[Theta]";
ToUnicode[ "\\vartheta" ] = "\[CurlyTheta]";
ToUnicode[ "\\iota" ] = "\[Iota]";
ToUnicode[ "\\kappa" ] = "\[Kappa]";
ToUnicode[ "\\lambda" ] = "\[Lambda]";
ToUnicode[ "\\mu" ] = "\[Mu]";
ToUnicode[ "\\nu" ] = "\[Nu]";
ToUnicode[ "\\xi" ] = "\[Xi]";
ToUnicode[ "\\pi" ] = "\[Pi]";
ToUnicode[ "\\varpi" ] = "\[CurlyPi]";
ToUnicode[ "\\rho" ] = "\[Rho]";
ToUnicode[ "\\varrho" ] = "\[CurlyRho]";
ToUnicode[ "\\sigma" ] = "\[Sigma]";
ToUnicode[ "\\varsigma" ] = "\[FinalSigma]";
ToUnicode[ "\\tau" ] = "\[Tau]";
ToUnicode[ "\\upsilon" ] = "\[Upsilon]";
ToUnicode[ "\\phi" ] = "\[Phi]";
ToUnicode[ "\\varphi" ] = "\[CurlyPhi]";
ToUnicode[ "\\chi" ] = "\[Chi]";
ToUnicode[ "\\psi" ] = "\[Psi]";
ToUnicode[ "\\omega" ] = "\[Omega]";
ToUnicode[ "\\Gamma" ] = "\[CapitalGamma]";
ToUnicode[ "\\Delta" ] = "\[CapitalDelta]";
ToUnicode[ "\\Theta" ] = "\[CapitalTheta]";
ToUnicode[ "\\Lambda" ] = "\[CapitalLambda]";
ToUnicode[ "\\Xi" ] = "\[CapitalXi]";
ToUnicode[ "\\Pi" ] = "\[CapitalPi]";
ToUnicode[ "\\Sigma" ] = "\[CapitalSigma]";
ToUnicode[ "\\Upsilon" ] = "\[CapitalUpsilon]";
ToUnicode[ "\\Phi" ] = "\[CapitalPhi]";
ToUnicode[ "\\Psi" ] = "\[CapitalPsi]";
ToUnicode[ "\\Omega" ] = "\[CapitalOmega]";
ToUnicode[ "\\infty" ] = "\[Infinity]";
ToUnicode[ "\\pm" ] = "\[PlusMinus]";
ToUnicode[ "\\partial" ] = "\[PartialD]";
ToUnicode[ "\\leq" ] = "\[LessEqual]";
ToUnicode[ "\\geq" ] = "\[GreaterEqual]";
ToUnicode[ "\\times" ] = "\[Times]";
ToUnicode[ "\\otimes" ] = "\[CircleTimes]";
ToUnicode[ "\\oplus" ] = "\[CirclePlus]";
ToUnicode[ "\\nabla" ] = "\[Del]";
ToUnicode[ "\\neq" ] = "\[NotEqual]";
ToUnicode[ "\\equiv" ] = "\[Congruent]";
ToUnicode[ "\\approx" ] = "\[TildeTilde]";
ToUnicode[ "\\ldots" ] = "\[Ellipsis]";
ToUnicode[ "\\in" ] = "\[Element]";
ToUnicode[ "\\notin" ] = "\[NotElement]";
ToUnicode[ "\\sim" ] = "\[Tilde]";
ToUnicode[ "\\sqrt" ] = "\[Sqrt]";
ToUnicode[ "\\propto" ] = "\[Proportional]";
ToUnicode[ "\\subset" ] = "\[Subset]";
ToUnicode[ "\\supset" ] = "\[Superset]";
ToUnicode[ "\\subseteq" ] = "\[SubsetEqual]";
ToUnicode[ "\\supseteq" ] = "\[SupersetEqual]";
ToUnicode[ "\\bullet" ] = "\[FilledSmallCircle]";
ToUnicode[ "\\perp" ] = "\[RightAngle]";
ToUnicode[ "\\simeq" ] = "\[TildeEqual]";
ToUnicode[ "\\vee" ] = "\[Vee]";
ToUnicode[ "\\wedge" ] = "\[Wedge]";
ToUnicode[ "\\leftrightarrow" ] = "\[LeftRightArrow]";
ToUnicode[ "\\leftarrow" ] = "\[LeftArrow]";
ToUnicode[ "\\rightarrow" ] = ToUnicode[ "\\to" ] = "\[RightArrow]";
ToUnicode[ "\\uparrow" ] = "\[UpArrow]";
ToUnicode[ "\\downarrow" ] = "\[DownArrow]";
ToUnicode[ "\\Leftarrow" ] = "\[DoubleLeftArrow]";
ToUnicode[ "\\Rightarrow" ] = "\[DoubleRightArrow]";
ToUnicode[ "\\Uparrow" ] = "\[DoubleUpArrow]";
ToUnicode[ "\\Downarrow" ] = "\[DoubleDownArrow]";
ToUnicode[ "\\bar" ] = "-";
ToUnicode[ "\\hat" ] = "^";
ToUnicode[ "\\tilde" ] = "\[Tilde]";
ToUnicode[ "\\dot" ] = "\[CenterDot]";
ToUnicode[ "\\ddot" ] = "\[CenterDot]\[CenterDot]";
ToUnicode[ "\\vec" ] = "\[RightVector]";
ToUnicode[ "\\prime" ] = "'";
ToUnicode[ "\\#" ] = "#";
ToUnicode[ "\\&" ] = "&";
ToUnicode[ "\\$" ] = "$";
ToUnicode[ "\\%" ] = "%";
ToUnicode[ "\\_" ] = "_";
ToUnicode[ Null ] = " ";
ToUnicode[ ComposedChar[t__] ] := ToUnicode/@ {t};
ToUnicode[ c_ ] := ToString[c]

End[]

