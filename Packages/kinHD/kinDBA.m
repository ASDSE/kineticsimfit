(* Wolfram Language Package *)

BeginPackage["kinHD`"]
(* Exported symbols added here with SymbolName::usage *)



kinHDParaCacheSignalall::usage = "";


Begin["`Private`"] (* Begin Private Context *)


(*kinetics IDA HD,D,HG,G,H*)

eqkin = {
   hd'[t] == -khd2*hd[t] + khd1*h[t]*d[t] ,
   h'[t]  == khd2*hd[t] - khd1*h[t]*d[t],
   d'[t]  == khd2*hd[t] - khd1*h[t]*d[t]};
(*signal*)
kinHDParaCacheSignalall[fh0_, fd0_, fsig0_, fsighd_, fsigd_, ftmin_, ftmax_] :=
  Block[{subHD, solv, h0c, d0c, g0c, hg0c, hd0c, initials, stepsize, startsignal},
    startsignal=fsig0+fsigd*fd0;
   initials = {h[0] == fh0, d[0] == fd0, hd[0] == 0, sig[0] == startsignal};
   solv = ParametricNDSolveValue[{eqkinsignal, initials}, sig, {t, ftmin, ftmax}, {khd1, khd2,sig0,sighd,sigd}];
   Return[solv];];


End[] (* End Private Context *)

EndPackage[]
