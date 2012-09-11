(*
	Convert31to32.m
		renames the ShapeData files of FeynArts 3.1
		for the new naming scheme of FeynArts 3.2
		last modified 1 Jul 03 th
*)


Needs["FeynArts`"]

shapeDir = ToFileName[Directory[], "ShapeData"]

If[ FileType[shapeDir] =!= Directory,
  Print["No directory ", shapeDir, " found"];
  Abort[]
]

backupDir = ToFileName[Directory[], "ShapeData-3.1"]

If[ FileType[backupDir] === Directory,
  Print["Conversion already performed (", backupDir, " exists)."];
  Abort[]
]


RenameFile[shapeDir, backupDir]


MkDir[dir_String] := dir /; FileType[dir] === Directory

MkDir[dir_String] := Check[CreateDirectory[dir], Abort[]]

MkDir[dirs__String] := Fold[MkDir[ToFileName[##]]&, {}, {dirs}]


newName[s_String] :=
Block[ {c, p, int = "0", out = "0", in = "0"},
  c = Characters[StringReplace[s, ".m" -> ""]];
  p = Position[c, "+"];
  If[ p =!= {},
    p = p[[1, 1]];
    in = str[Drop[c, p]];
    c = Take[c, p - 1];
  ];
  p = Position[c, "-"];
  If[ p =!= {},
    p = p[[1, 1]];
    out = str[Drop[c, p]];
    c = Take[c, p - 1];
  ];
  int = str[c] <> ".m";
  CopyFile[s, ToFileName[MkDir[shapeDir, in, out], int]]
]

str[c_] := StringJoin[Partition[c, 2]//Reverse] /. "" -> "0"



SetDirectory[backupDir]

newName/@ FileNames[]

