% simple_cross.mp -- do not edit, generated automatically by feynSample.tex
input feynmp
require_RCS_revision "1.30";
vleft (__i1, __i2);
vright (__o1, __o2);
vconnect ("fermion", __i1, __v1, __o2);
vconnect ("fermion", __i2, __v2, __o1);
vconnect ("photon", __v1, __v2);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
