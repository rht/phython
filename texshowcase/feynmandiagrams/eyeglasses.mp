% eyeglasses.mp -- do not edit, generated automatically by feynSample.tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 50*1.0pt#, 30*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vleft (__i);
vright (__o);
vconnect ("phantom,tension=10", __i, __i1);
vconnect ("phantom,tension=10", __o, __o1);
vconnect ("fermion,left,tension=0.4", __i1, __v1, __i1);
vconnect ("fermion,right,tension=0.4", __o1, __v2, __o1);
vconnect ("dashes", __v1, __v2);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
