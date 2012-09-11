(*
	SMNLG.mod
		Add-on model file for the SM in Non-Linear Gauge
		needs Lorentzbgf.gen
		based on arXiv:0710.1999
		by Thomas Gajdosik and Jurgis Pasukonis
		last modified 1 Apr 12 by Thomas Hahn
*)


LoadModel["SM"]


ReplaceCouplings[

	(* structural changes re Lorentzbgf.mod *)

  C[ s1 V[j1], s2 V[j2], s3 V[j3] ] == Join[
    C[ s1 V[j1], s2 V[j2], s3 V[j3] ],
    {{0, 0}, {0, 0}, {0, 0}} ],

  C[ s1 S[j1], s2 S[j2], s3 V[j3] ] == Join[
    C[ s1 S[j1], s2 S[j2], s3 V[j3] ],
   -C[ s1 S[j1], s2 S[j2], s3 V[j3] ] ]

]


ReplaceCouplings[

	(* V-V-V-V *)

  C[ -V[3], V[3], V[2], V[2] ] += -I EL^2 CW^2/SW^2 *
    Gbeta^2/GaugeXi[W] {{0, 0}, {1, 0}, {1, 0}},

  C[ -V[3], V[3], V[1], V[2] ] += I EL^2 CW/SW *
    Galpha Gbeta/GaugeXi[W] {{0, 0}, {1, 0}, {1, 0}},

  C[ -V[3], V[3], V[1], V[1] ] += -I EL^2 *
    Galpha^2/GaugeXi[W] {{0, 0}, {1, 0}, {1, 0}},

	(* V-V-V *)

  C[ V[1], -V[3], V[3] ] += -I EL *
    Galpha/GaugeXi[W] {{0, 0}, {1, 0}, {-1, 0}, {0, 0}},

  C[ V[2], -V[3], V[3] ] += I EL CW/SW *
    Gbeta/GaugeXi[W] {{0, 0}, {1, 0}, {-1, 0}, {0, 0}},

	(* S-S-S-S *)

  C[ S[1], S[1], S[2], S[2] ] += -I EL^2 /(4 SW^2 MW^2) *
    2 Gepsilon^2 MZ^2 GaugeXi[Z] {{1, 0}},

  C[ S[1], S[1], S[3], -S[3] ] += -I EL^2 /(4 SW^2 MW^2) *
    2 Gdelta^2 MW^2 GaugeXi[W] {{1, 0}},

  C[ S[2], S[2], S[3], -S[3] ] += -I EL^2 /(4 SW^2 MW^2) *
    2 Gkappa^2 MW^2 GaugeXi[W] {{1, 0}},

	(* S-S-S *)

  C[ S[1], S[2], S[2] ] += -I EL /(2 SW MW) *
    2 Gepsilon MZ^2 GaugeXi[Z] {{1, 0}},

  C[ S[3], S[1], -S[3] ] += -I EL /(2 SW MW) *
    2 Gdelta MW^2 GaugeXi[W] {{1, 0}},

	(* S-S-V-V *)

  C[ S[1], -S[3], V[3], V[2] ] += -I EL^2/(2 CW) *
    CW^2/SW^2 Gbeta Gdelta {{1, 0}},

  C[ S[1], S[3], -V[3], V[2] ] += -I EL^2/(2 CW) *
    CW^2/SW^2 Gbeta Gdelta {{1, 0}},

  C[ S[1], S[3], -V[3], V[1] ] += -I EL^2/(2 SW) *
    (-Galpha Gdelta) {{1, 0}},

  C[ S[1], -S[3], V[3], V[1] ] += -I EL^2/(2 SW) *
    (-Galpha Gdelta) {{1, 0}},

  C[ S[3], S[2], V[2], -V[3] ] += EL^2/(2 CW) *
    CW^2/SW^2 Gbeta Gkappa {{1, 0}},

  C[ -S[3], S[2], V[2], V[3] ] += -EL^2/(2 CW) *
    CW^2/SW^2 Gbeta Gkappa {{1, 0}},

  C[ S[3], S[2], V[1], -V[3] ] += EL^2/(2 SW) *
    (-Galpha Gkappa) {{1, 0}},

  C[ -S[3], S[2], V[1], V[3] ] += -EL^2/(2 SW) *
    -(Galpha Gkappa) {{1, 0}},

	(* S-S-V *)

  C[ S[2], S[1], V[2] ] += EL/(2 CW SW) *
    Gepsilon {{-1, 0}, {-1, 0}},

  C[ S[3], S[1], -V[3] ] += -I EL/(2 SW) *
    Gdelta {{-1, 0}, {-1, 0}},

  C[ -S[3], S[1], V[3] ] += I EL/(2 SW) *
    Gdelta {{-1, 0}, {-1, 0}},

  C[ S[3], S[2], -V[3] ] += EL/(2 SW) *
    Gkappa {{-1, 0}, {-1, 0}},

  C[ -S[3], S[2], V[3] ] += EL/(2 SW) *
    Gkappa {{-1, 0}, {-1, 0}},

	(* S-V-V *)

  C[ -S[3], V[3], V[2] ] += -I EL MW SW/CW *
    CW^2/SW^2 Gbeta {{1, 0}},

  C[ S[3], -V[3], V[2] ] += -I EL MW SW/CW *
    CW^2/SW^2 Gbeta {{1, 0}},

  C[ -S[3], V[3], V[1] ] += -I EL MW * 
    (-Galpha) {{1, 0}},

  C[ S[3], -V[3], V[1] ] += -I EL MW * 
    (-Galpha) {{1, 0}},

	(* U-U-V:  G(+) . { p1_mu3, p2_mu3 } *)

  C[ -U[3], U[3], V[1] ] += -I EL *
    (-Galpha) {{0, 0}, {1, 0}},

  C[ -U[4], U[4], V[1] ] += I EL *
    (-Galpha) {{0, 0}, {1, 0}},

  C[ -U[3], U[3], V[2] ] += I EL CW/SW *
    (-Gbeta) {{0, 0}, {1, 0}},

  C[ -U[4], U[4], V[2] ] += -I EL CW/SW *
    (-Gbeta) {{0, 0}, {1, 0}},

  C[ -U[3], U[2], V[3] ] += -I EL CW/SW *
    Gbeta {{0, 0}, {1, 0}},

  C[ -U[4], U[2], -V[3] ] += I EL CW/SW *
    Gbeta {{0, 0}, {1, 0}},

  C[ -U[3], U[1], V[3] ] += I EL *
    Galpha {{0, 0}, {1, 0}},

  C[ -U[4], U[1], -V[3] ] += -I EL *
    Galpha {{0, 0}, {1, 0}},

	(* S-U-U *)

  C[ S[1], -U[2], U[2] ] += -I EL MZ GaugeXi[Z]/(2 SW CW) *
    Gepsilon {{1, 0}},

  C[ S[1], -U[3], U[3] ] += -I EL MW GaugeXi[W]/(2 SW) *
    Gdelta {{1, 0}},

  C[ S[1], -U[4], U[4] ] += -I EL MW GaugeXi[W]/(2 SW) *
    Gdelta {{1, 0}},

  C[ S[2], -U[4], U[4] ] += EL MW GaugeXi[W]/(2 SW) *
    (-Gkappa) {{1, 0}},

  C[ S[2], -U[3], U[3] ] += -EL MW GaugeXi[W]/(2 SW) *
    (-Gkappa) {{1, 0}},

  C[ -S[3], -U[4], U[2] ] += I EL MW GaugeXi[W]/(2 CW SW) *
    (-Gkappa) {{1, 0}},

  C[ S[3], -U[3], U[2] ] += I EL MW GaugeXi[W]/(2 CW SW) *
    (-Gkappa) {{1, 0}}

]


M$CouplingMatrices = Join[ M$CouplingMatrices, {

	(* U-U-V-V *)

  C[ -U[4], U[1], V[1], -V[3] ] == -I EL^2 * 
    {{Galpha}}, 

  C[ -U[3], U[1], V[1], V[3] ] == -I EL^2 * 
    {{Galpha}}, 

  C[ -U[4], U[1], V[2], -V[3] ] == I EL^2 (CW/SW) * 
    {{Gbeta}}, 

  C[ -U[3], U[1], V[2], V[3] ] == I EL^2 (CW/SW) * 
    {{Gbeta}}, 

  C[ -U[4], U[2], V[1], -V[3] ] == I EL^2 (CW/SW) * 
    {{Galpha}}, 

  C[ -U[3], U[2], V[1], V[3] ] == I EL^2 (CW/SW) * 
    {{Galpha}}, 

  C[ -U[4], U[2], V[2], -V[3] ] == -I EL^2 (CW^2/SW^2) * 
    {{Gbeta}}, 

  C[ -U[3], U[2], V[2], V[3] ] == -I EL^2 (CW^2/SW^2) * 
    {{Gbeta}}, 

  C[ -U[3], U[3], V[3], -V[3] ] == -I EL^2 * 
    {{Galpha + CW^2/SW^2 Gbeta}},

  C[ -U[4], U[4], V[3], -V[3] ] == -I EL^2 * 
    {{Galpha + CW^2/SW^2*Gbeta}},

  C[ -U[3], U[4], V[3], V[3] ] == 2 I EL^2 * 
    {{Galpha + CW^2/SW^2 Gbeta}},

  C[ -U[4], U[3], -V[3], -V[3] ] == 2 I EL^2 * 
    {{Galpha + CW^2/SW^2 Gbeta}},

  C[ -U[3], U[3], V[1], V[1] ] == 2 I EL^2 * 
    {{Galpha}},

  C[ -U[4], U[4], V[1], V[1] ] == 2 I EL^2 * 
    {{Galpha}},

  C[ -U[3], U[3], V[1], V[2] ] == -I EL^2 (CW/SW) * 
    {{Galpha + Gbeta}},

  C[ -U[4], U[4], V[1], V[2] ] == -I EL^2 (CW/SW) * 
    {{Galpha + Gbeta}},

  C[ -U[3], U[3], V[2], V[2] ] == 2 I EL^2 (CW^2/SW^2) * 
    {{Gbeta}},

  C[ -U[4], U[4], V[2], V[2] ] == 2 I EL^2 (CW^2/SW^2) * 
    {{Gbeta}},

	(* S-S-U-U *)

  C[ S[1], S[1], -U[2], U[2] ] == -I EL^2 GaugeXi[Z]/(2 SW^2 CW^2) *
    {{Gepsilon}},

  C[ S[2], S[2], -U[2], U[2] ] == I EL^2 GaugeXi[Z]/(2 SW^2 CW^2) *
    {{Gepsilon}},

  C[ S[3], S[1], -U[2], U[4] ] == I EL^2 GaugeXi[Z]/(4 SW^2 CW) *
    {{Gepsilon}},

  C[ -S[3], S[1], -U[2], U[3] ] == I EL^2 GaugeXi[Z]/(4 SW^2 CW) *
    {{Gepsilon}},

  C[ S[3], S[2], -U[2], U[4] ] == EL^2 GaugeXi[Z]/(4 SW^2 CW) *
    {{Gepsilon}},

  C[ -S[3], S[2], -U[2], U[3] ] == - EL^2 GaugeXi[Z]/(4 SW^2 CW) *
    {{Gepsilon}},

  C[ -S[3], S[1], -U[4], U[1] ] == I EL^2 GaugeXi[W]/(2 SW) *
    {{Gdelta}},

  C[ S[3], S[1], -U[3], U[1] ] == I EL^2 GaugeXi[W]/(2 SW) *
    {{Gdelta}},

  C[ -S[3], S[2], -U[4], U[1] ] == EL^2 GaugeXi[W]/(2 SW) *
    {{Gkappa}},

  C[ S[3], S[2], -U[3], U[1] ] == -EL^2 GaugeXi[W]/(2 SW) *
    {{Gkappa}},

  C[ -S[3], S[1], -U[4], U[2] ] == -I EL^2 GaugeXi[W]/(4 SW^2 CW) *
    {{ Gkappa + Gdelta (CW^2 - SW^2) }},

  C[ S[3], S[1], -U[3], U[2] ] == -I EL^2 GaugeXi[W]/(4 SW^2 CW) *
    {{ Gkappa + Gdelta (CW^2 - SW^2) }},

  C[ -S[3], S[2], -U[4], U[2] ] == -EL^2 GaugeXi[W]/(4 SW^2 CW) *
    {{ Gdelta + Gkappa (CW^2 - SW^2) }},

  C[ S[3], S[2], -U[3], U[2] ] == EL^2 GaugeXi[W]/(4 SW^2 CW) *
    {{ Gdelta + Gkappa (CW^2 - SW^2) }},

  C[ S[1], S[1], -U[3], U[3] ] == -I EL^2 GaugeXi[W]/(2 SW^2) *
    {{ Gdelta }},

  C[ S[1], S[1], -U[4], U[4] ] == -I EL^2 GaugeXi[W]/(2 SW^2) *
    {{ Gdelta }},

  C[ S[2], S[2], -U[3], U[3] ] == -I EL^2 GaugeXi[W]/(2 SW^2) *
    {{ Gkappa }},

  C[ S[2], S[2], -U[4], U[4] ] == -I EL^2 GaugeXi[W]/(2 SW^2) *
    {{ Gkappa }},

  C[ S[2], S[1], -U[3], U[3] ] == EL^2 GaugeXi[W]/(4 SW^2) *
    {{ Gkappa - Gdelta }},

  C[ S[2], S[1], -U[4], U[4] ] == -EL^2 GaugeXi[W]/(4 SW^2) *
    {{ Gkappa - Gdelta }},

  C[ S[3], -S[3], -U[3], U[3] ] == I EL^2 GaugeXi[W]/(4 SW^2) *
    {{ Gkappa + Gdelta }},

  C[ S[3], -S[3], -U[4], U[4] ] == I EL^2 GaugeXi[W]/(4 SW^2) *
    {{ Gkappa + Gdelta }},

  C[ S[3], S[3], -U[3], U[4] ] == I EL^2 GaugeXi[W]/(2 SW^2) *
    {{ Gkappa - Gdelta }},

  C[ -S[3], -S[3], -U[4], U[3] ] == I EL^2 GaugeXi[W]/(2 SW^2) *
    {{ Gkappa - Gdelta }}

} ]

