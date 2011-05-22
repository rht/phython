% ears.mp -- do not edit, generated automatically by feynSample.tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 50*1.0pt#, 30*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vtop (__t0, __t1, __t2, __t3);
vbottom (__b0, __b1, __b2, __b3);
vconnect ("phantom", __t1, __v1, __b1);
vconnect ("phantom", __t2, __v2, __b2);
vfreeze;
vconnect ("plain,right", __v1, __v2, __v1);
vconnect ("plain,tension=0.8,right=270", __v1, __v1);
vconnect ("plain,tension=0.8,right", __v2, __v2);
vdot (__v1, __v2);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
