* communication with external channel
* should be in CalcFeynAmp.frm
* presently in own file until Windows version of FORM gets fixed
* (#toexternal not ignored inside #if 0 ... #endif)


b Den, intM;
.sort

keep brackets;

#setexternal `PIPE1_'

#call MapAll(ToMma)
#toexternal "#"

drop;
.sort

#fromexternal
.global

#call MapAll(FromMma)
