(* Wolfram Language Package *)

BeginPackage["kinHD`"]
(* Exported symbols added here with SymbolName::usage *)



kinHDParaCacheHD::usage = "";
kinHDParaCacheD::usage = "";
kinHDParaCacheH::usage = "";


Begin["`Private`"] (* Begin Private Context *)


(*kinetics IDA HD,D,HG,G,H*)

eqkin = {
   hd'[t] == -khd2*hd[t] + khd1*h[t]*d[t] ,
   h'[t]  == khd2*hd[t] - khd1*h[t]*d[t],
   d'[t]  == khd2*hd[t] - khd1*h[t]*d[t]};
(*signal*)
kinHDParaCacheHD[fh0_, fd0_, ftmin_, ftmax_] := kinHDParaCacheHD[fh0, fd0, ftmin, ftmax] =
  Block[{subHD, solv, h0c, d0c, g0c, hg0c, hd0c, initials, stepsize, startsignal},
   initials = {h[0] == fh0, d[0] == fd0, hd[0] == 0};
   solv = ParametricNDSolveValue[{eqkin, initials}, hd, {t, ftmin, ftmax}, {khd1, khd2}];
   Return[solv];];

kinHDParaCacheD[fh0_, fd0_, ftmin_, ftmax_] := kinHDParaCacheD[fh0, fd0, ftmin, ftmax] =
  Block[{subHD, solv, h0c, d0c, g0c, hg0c, hd0c, initials, stepsize, startsignal},
   initials = {h[0] == fh0, d[0] == fd0, hd[0] == 0};
   solv = ParametricNDSolveValue[{eqkin, initials}, d, {t, ftmin, ftmax}, {khd1, khd2}];
   Return[solv];];

kinHDParaCacheH[fh0_, fd0_, ftmin_, ftmax_] := kinHDParaCacheH[fh0, fd0, ftmin, ftmax] =
  Block[{subHD, solv, h0c, d0c, g0c, hg0c, hd0c, initials, stepsize, startsignal},
   initials = {h[0] == fh0, d[0] == fd0, hd[0] == 0};
   solv = ParametricNDSolveValue[{eqkin, initials}, h, {t, ftmin, ftmax}, {khd1, khd2}];
   Return[solv];];


End[] (* End Private Context *)

EndPackage[]
