* run.F
* defines parameter settings for one run, then includes main.F
* this file is part of FormCalc
* last modified 7 Mar 11 th


* Whether to run in debugging mode

#define DEBUG


* Uncomment the following to check UV- and IR-finiteness

c#define DELTA 1D7
c#define MUDIM 1D100
c#define LAMBDA 1D10

* The following sets all ren. constants to zero for debugging
* or if calculating a tree-level reaction.

c#define NO_RENCONST


* The model initialization code is included here.  Options for the
* model_XXX.F file (e.g. NO_SQUARK_MIXING for model_mssm.F) must be
* defined before the #include line.

* options for model_sm.F:
c#define NO_CKM

* options for model_mssm.F:
c#define NO_SQUARK_MIXING
c#define HIGGS_MASSES TREE
c#define HIGGS_MASSES SIMPLE
c#define HIGGS_MASSES FEYNHIGGS
c#define SM_HIGGS
c#define NMFV
c#define NO_EXCLUSION_LIMITS
c#define NO_GUT_RELATION

* options for model_thdm.F:
c#define MODEL_TYPE_I
c#define MODEL_TYPE_II
c#define NO_EXCLUSION_LIMITS

#define MODEL "model_mssm.F"


* The LOOPn set up the model parameters for the calculation, either
* - as fixed declarations, e.g.
*       MH = 115
* - or as do-loops terminating on label 1, e.g.
*       do 1 MH = 100, 500, 10

#define LOOP1 TB = 1.5D0
#define LOOP2 MA0 = 700
#define LOOP3 MSusy = 1000
#define LOOP4 MUE = -1000
#define LOOP5 Mino2 = 100
#define LOOP6 At = MSusy
#define LOOP7 Ab = At
#define LOOP8 Atau = At
#define LOOP9
#define LOOP10
#define LOOP11
#define LOOP12
#define LOOP13
#define LOOP14
#define LOOP15
#define LOOP16
#define LOOP17
#define LOOP18
#define LOOP19
#define LOOP20


* The PRINTn definitions contain output statements with which changes
* in parameters (typically the ones included in the scan) can be
* written to the output file and screen.  For example,
*       SHOW "MH", MH
* writes out the values of MH and MT in each pass of the loop.

#define PRINT1 SHOW "TB", TB
#define PRINT2 SHOW "MA0", MA0
#define PRINT3 SHOW "M_2", ReM_2
#define PRINT4 SHOW "MUE", ReMUE
#define PRINT5 SHOW "MSusy", MSusy
#define PRINT6 SHOW "At", ReAt
#define PRINT7
#define PRINT8
#define PRINT9
#define PRINT10
#define PRINT11
#define PRINT12
#define PRINT13
#define PRINT14
#define PRINT15
#define PRINT16
#define PRINT17
#define PRINT18
#define PRINT19
#define PRINT20


* Cuts and increments.  May be defined in terms of sqrtS, the CMS energy.
* Note: cuts are numbered by the integration they belong to, which is
* from the last particle downward, e.g. COSTH1CUT refers to cos(th) of
* the *last* particle.

#define SQRTSSTEP 10

#define COSTH1CUT .01D0
#ifdef DEBUG
#define XCOSTH1STEP 1/3D0
#else
#define XCOSTH1STEP .02D0
#endif

#define XMREM1STEP .1D0
#define E1MIN ESOFTMAX


* The integration parameters:
* - the desired relative and absolute accuracy

#define RELACCURACY 1D-3
#define ABSACCURACY 1D-10

* The multidimensional integration routines also need

#define METHOD DIVONNE
#define USERDATA bogus
#define VERBOSE 1
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 50000

* for Vegas:
#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE ""

* for Suave:
#define NNEW 1000
#define FLATNESS 50

* for Divonne:
#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 1D-6
#define MAXCHISQ 10
#define MINDEVIATION .25D0

* for Cuhre:
#define KEY 0

* Now go for it:

#include "main.F"

