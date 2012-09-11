(*
	Convert3to31.m
		rewrites the GraphInfo files of FeynArts 3
		as individual files for each topology in the
		ShapeData directory
		last modified 10 Mar 03 th
*)


Needs["FeynArts`"]

If[ FileType[$ShapeDataDir] === Directory,
  Print["Conversion already performed (", $ShapeDataDir, " exists)."];
  Abort[]
]


If[ !ValueQ[$TopologyDataDir],
  $TopologyDataDir = ToFileName[{$FeynArtsDir, "GraphInfo"}] ]


sep = StringTake[ToFileName["d", "f"], {2, -2}]

BaseName[file_] :=
Block[ {p1, p2},
  p1 = Flatten[{0, StringPosition[file, sep]}][[-1]] + 1;
  p2 = Flatten[{0, StringPosition[file, "."]}][[-1]] - 1;
  If[ p2 < p1, p2 = -1 ];
  StringTake[file, {p1, p2}]
]


LoadOld[file_] := (
  Get[file];
  ToExpression[BaseName[file]]
)

cats = LoadOld/@ FileNames["*", $TopologyDataDir]


SaveNew[ _[cat_[top_]], shapedata_ ] :=
Block[ {topcode = TopologyCode[top], shapefile},
  If[ Head[tag[topcode]] =!= tag,
    Print[topcode, " in more than one category: ", tag[topcode], " and ", cat];
    Abort[]
  ];
  tag[topcode] = cat;
  shapefile = ToFileName[$ShapeDataDir, topcode <> ".m"];
  Put[shapedata, shapefile];
]

CreateDirectory[$ShapeDataDir]

Apply[SaveNew, DownValues/@ cats, {2}];

