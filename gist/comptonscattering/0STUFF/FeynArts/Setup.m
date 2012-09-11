(*
	Setup.m
		FeynArts startup file
		last modified 17 Mar 08 th

Here you can set up your own changes and enhancements to FeynArts,
e.g. some particular options you always want set, or $SVMixing = True.
It is a good idea to do this here since changing the FeynArts code 
directly is inherently unportable.

*)


$Verbose = 2

$ModelPath = { Directory[],
  ToFileName[{Directory[], "Models"}],
  ToFileName[{$FeynArtsDir, "Models"}] }

$ModelDebug = False

$ModelDebugForm = Short[#, 5]&

$ShapeDataDir = ToFileName[{$FeynArtsDir, "ShapeData"}]

$SVMixing = False

$CounterTerms = True

$FermionLines = True

	(* eliminate those `>' in front of continuation lines so
	   one can cut and paste more easily *)
Format[ Continuation[_] ] = "    "

