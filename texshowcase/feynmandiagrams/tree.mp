% tree.mp -- do not edit, generated automatically by feynSample.tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 40*1.0pt#, 30*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vbottom (__i1, __d1, __o1);
vright (__o0, __o2, __o3);
vconnect ("fermion", __i1, __v1, __o1);
vfreeze;
vconnect ("fermion", __o2, __v2, __o3);
vconnect ("photon,tension=1.5", __v1, __v2);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
