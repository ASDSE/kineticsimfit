(* Wolfram Language Package *)

BeginPackage["kinHD`"]
(* Exported symbols added here with SymbolName::usage *)



kinHDParaCacheSignalall::usage = "";


Begin["`Private`"] (* Begin Private Context *)
eqthermo = {h0 == h + hd, d0 == d + hd, kd == hd/(h*d)};

sthdCacheHDKdD0[fkd_?NumericQ, fd0_?NumericQ] := sthdCacheHDKdD0[fkd, fd0] =
   Block[{eli, solv},
    eli = Eliminate[eqthermo /. {d0 -> fd0, kd -> fkd}, {h, d}];
    solv = hd /. NSolve[eli, hd];
    Return[solv];];

sthdCacheResultHD[fkd_?NumericQ, fd0_?NumericQ, fh0_?NumericQ] := sthdCacheResultHD[fkd, fd0, fh0] =
  Block[{sol, e, hs0, ds0, brac},
   If[fh0 == 0, hs0 = 1.*10^-15, hs0 = fh0];
   If[fd0 == 0, ds0 = 1.*10^-15, ds0 = fd0];
   sol = sthdCacheHDKdD0[fkd, ds0] /. {h0 -> hs0};
   e = Select[sol, 0 <= # <= hs0 && 0 <= # <= ds0 &];
   brac = First[e];
   Return[brac];];

IDAStart[fkd_?NumericQ, fd0_?NumericQ, fh0_?NumericQ, fkga0_?NumericQ, fga0_?NumericQ] := IDAStart[fkd, fd0, fh0, fkga0, fga0] =
  Block[{h0c, d0c, g0c, hg0c, hd0c, initials, stepsize, sol},
   hd0c = sthdCacheResultHD[fkd, fd0, fh0];
   h0c = fh0 - hd0c;
   d0c = fd0 - hd0c;
   hg0c = 0;
   g0c = fga0;
   sol = {hd0c, h0c, d0c, hg0c, g0c};
   Return[N[sol]]
   ];

GDAStart[fkd_?NumericQ, fd0_?NumericQ, fh0_?NumericQ, fkga0_?NumericQ, fga0_?NumericQ] := GDAStart[fkd, fd0, fh0, fkga0, fga0] =
 Block[{h0c, d0c, g0c, hg0c, hd0c, initials, stepsize, sol},
  hd0c = 0;
  hg0c = sthdCacheResultHD[fkga0, fga0, fh0];
  g0c = fga0 - hg0c;
	h0c = fh0 - hg0c;
	d0c = fd0;
  sol = {hd0c, h0c, d0c, hg0c, g0c};
  Return[N[sol]]
  ];

SBAStart[fkd_?NumericQ, fd0_?NumericQ, fh0_?NumericQ, fkga0_?NumericQ, fga0_?NumericQ] := SBAStart[fkd, fd0, fh0, fkga0, fga0] =
 Block[{h0c, d0c, g0c, hg0c, hd0c, initials, stepsize, sol},
  hd0c = 0;
  hg0c = 0;
  g0c = fga0;
	h0c = fh0;
	d0c = fd0;
  sol = {hd0c, h0c, d0c, hg0c, g0c};
  Return[N[sol]]
  ];

(*kinetics IDA HD,D,HG,G,H*)

eqkinsignal = {
   hd'[t] == (-khd1/kequhd)*hd[t] + khd1*h[t]*d[t] ,
   h'[t] == (khd1/kequhd)*hd[t] - khd1*h[t]*d[t],
   d'[t] == khd1/kequhd*hd[t] - khd1*h[t]*d[t],
  sig[t] == sig0 + sighd*hd[t] + sigd*d[t]};
(*signal*)
kinHDParaCacheSignalall[fkd_, fk1hd_, fh0_, fd0_, fsig0_, fsighd_, fsigd_, ftmin_, ftmax_] := kinHDParaCacheSignalall[fkd, fk1hd, fh0, fd0, fsig0, fsighd, fsigd, ftmin, ftmax] =
  Block[{subHD, solv, h0c, d0c, g0c, hg0c, hd0c, initials, stepsize},
   initials = {h[0] == fh0, d[0] == fd0, hd[0] == 0};
   solv = ParametricNDSolveValue[{eqkinsignal, initials}, sig, {t, ftmin, ftmax}, {khd1, kequhd,sig0,sighd,sigd}];
   Return[solv];];

kinIDAParaCacheDall[fkd_, fkga_, fk1hd_, fk1hga_, fh0_, fd0_, fga0_, ftmin_, ftmax_] := kinIDAParaCacheDall[fkd, fkga, fk1hd, fk1hga, fh0, fd0, fga0, ftmin, ftmax] =
  Block[{subHD, solv, h0c, d0c, g0c, hg0c, hd0c, initials, stepsize},
   {hd0c, h0c, d0c, hg0c, g0c} = IDAStart[fkd, fd0, fh0, fkga, fga0];
   initials = {h[0] == h0c, g[0] == g0c, d[0] == d0c, hg[0] == hg0c, hd[0] == hd0c};
   solv = ParametricNDSolveValue[{eqkin, initials}, d, {t, ftmin, ftmax}, {khd1, khg1, kequhd, kequhg}];
   Return[solv];];

kinIDAParaCacheHall[fkd_, fkga_, fk1hd_, fk1hga_, fh0_, fd0_, fga0_, ftmin_, ftmax_] := kinIDAParaCacheHall[fkd, fkga, fk1hd, fk1hga, fh0, fd0, fga0, ftmin, ftmax] =
  Block[{subHD, solv, h0c, d0c, g0c, hg0c, hd0c, initials, stepsize},
   {hd0c, h0c, d0c, hg0c, g0c} = IDAStart[fkd, fd0, fh0, fkga, fga0];
   initials = {h[0] == h0c, g[0] == g0c, d[0] == d0c, hg[0] == hg0c, hd[0] == hd0c};
   solv = ParametricNDSolveValue[{eqkin, initials}, h, {t, ftmin, ftmax}, {khd1, khg1, kequhd, kequhg}];
   Return[solv];];

kinIDAParaCacheHGall[fkd_, fkga_, fk1hd_, fk1hga_, fh0_, fd0_, fga0_, ftmin_, ftmax_] := kinIDAParaCacheHGall[fkd, fkga, fk1hd, fk1hga, fh0, fd0, fga0, ftmin, ftmax] =
  Block[{subHD, solv, h0c, d0c, g0c, hg0c, hd0c, initials, stepsize},
   {hd0c, h0c, d0c, hg0c, g0c} = IDAStart[fkd, fd0, fh0, fkga, fga0];
   initials = {h[0] == h0c, g[0] == g0c, d[0] == d0c, hg[0] == hg0c, hd[0] == hd0c};
   solv = ParametricNDSolveValue[{eqkin, initials}, hg, {t, ftmin, ftmax}, {khd1, khg1, kequhd, kequhg}];
   Return[solv];];

kinIDAParaCacheGall[fkd_, fkga_, fk1hd_, fk1hga_, fh0_, fd0_, fga0_, ftmin_, ftmax_] := kinIDAParaCacheGall[fkd, fkga, fk1hd, fk1hga, fh0, fd0, fga0, ftmin, ftmax] =
  Block[{subHD, solv, h0c, d0c, g0c, hg0c, hd0c, initials, stepsize},
   {hd0c, h0c, d0c, hg0c, g0c} = IDAStart[fkd, fd0, fh0, fkga, fga0];
   initials = {h[0] == h0c, g[0] == g0c, d[0] == d0c, hg[0] == hg0c, hd[0] == hd0c};
   solv = ParametricNDSolveValue[{eqkin, initials}, g, {t, ftmin, ftmax}, {khd1, khg1, kequhd, kequhg}];
   Return[solv];];

(*GDA*)
kinGDAParaCacheHDall[fkd_, fkga_, fk1hd_, fk1hga_, fh0_, fd0_, fga0_, ftmin_, ftmax_] := kinGDAParaCacheHDall[fkd, fkga, fk1hd, fk1hga, fh0, fd0, fga0, ftmin, ftmax] =
 Block[{subHD, solv, h0c, d0c, g0c, hg0c, hd0c, initials, stepsize},
  {hd0c, h0c, d0c, hg0c, g0c} = GDAStart[fkd, fd0, fh0, fkga, fga0];
  initials = {h[0] == h0c, g[0] == g0c, d[0] == d0c, hg[0] == hg0c, hd[0] == hd0c};
  solv = ParametricNDSolveValue[{eqkin, initials}, hd, {t, ftmin, ftmax}, {khd1, khg1, kequhd, kequhg}];
  Return[solv];];


kinGDAParaCacheDall[fkd_, fkga_, fk1hd_, fk1hga_, fh0_, fd0_, fga0_, ftmin_, ftmax_] := kinGDAParaCacheDall[fkd, fkga, fk1hd, fk1hga, fh0, fd0, fga0, ftmin, ftmax] =
Block[{subHD, solv, h0c, d0c, g0c, hg0c, hd0c, initials, stepsize},
{hd0c, h0c, d0c, hg0c, g0c} = GDAStart[fkd, fd0, fh0, fkga, fga0];
initials = {h[0] == h0c, g[0] == g0c, d[0] == d0c, hg[0] == hg0c, hd[0] == hd0c};
solv = ParametricNDSolveValue[{eqkin, initials}, d, {t, ftmin, ftmax}, {khd1, khg1, kequhd, kequhg}];
Return[solv];];

kinGDAParaCacheHall[fkd_, fkga_, fk1hd_, fk1hga_, fh0_, fd0_, fga0_, ftmin_, ftmax_] := kinGDAParaCacheHall[fkd, fkga, fk1hd, fk1hga, fh0, fd0, fga0, ftmin, ftmax] =
Block[{subHD, solv, h0c, d0c, g0c, hg0c, hd0c, initials, stepsize},
{hd0c, h0c, d0c, hg0c, g0c} = GDAStart[fkd, fd0, fh0, fkga, fga0];
initials = {h[0] == h0c, g[0] == g0c, d[0] == d0c, hg[0] == hg0c, hd[0] == hd0c};
solv = ParametricNDSolveValue[{eqkin, initials}, h, {t, ftmin, ftmax}, {khd1, khg1, kequhd, kequhg}];
Return[solv];];

kinGDAParaCacheHGall[fkd_, fkga_, fk1hd_, fk1hga_, fh0_, fd0_, fga0_, ftmin_, ftmax_] := kinGDAParaCacheHGall[fkd, fkga, fk1hd, fk1hga, fh0, fd0, fga0, ftmin, ftmax] =
Block[{subHD, solv, h0c, d0c, g0c, hg0c, hd0c, initials, stepsize},
{hd0c, h0c, d0c, hg0c, g0c} = GDAStart[fkd, fd0, fh0, fkga, fga0];
initials = {h[0] == h0c, g[0] == g0c, d[0] == d0c, hg[0] == hg0c, hd[0] == hd0c};
solv = ParametricNDSolveValue[{eqkin, initials}, hg, {t, ftmin, ftmax}, {khd1, khg1, kequhd, kequhg}];
Return[solv];];

kinGDAParaCacheGall[fkd_, fkga_, fk1hd_, fk1hga_, fh0_, fd0_, fga0_, ftmin_, ftmax_] := kinGDAParaCacheGall[fkd, fkga, fk1hd, fk1hga, fh0, fd0, fga0, ftmin, ftmax] =
Block[{subHD, solv, h0c, d0c, g0c, hg0c, hd0c, initials, stepsize},
{hd0c, h0c, d0c, hg0c, g0c} = GDAStart[fkd, fd0, fh0, fkga, fga0];
initials = {h[0] == h0c, g[0] == g0c, d[0] == d0c, hg[0] == hg0c, hd[0] == hd0c};
solv = ParametricNDSolveValue[{eqkin, initials}, g, {t, ftmin, ftmax}, {khd1, khg1, kequhd, kequhg}];
Return[solv];];

(*SBA*)
kinSBAParaCacheHDall[fkd_, fkga_, fk1hd_, fk1hga_, fh0_, fd0_, fga0_, ftmin_, ftmax_] := kinSBAParaCacheHDall[fkd, fkga, fk1hd, fk1hga, fh0, fd0, fga0, ftmin, ftmax] =
 Block[{subHD, solv, h0c, d0c, g0c, hg0c, hd0c, initials, stepsize},
  {hd0c, h0c, d0c, hg0c, g0c} = SBAStart[fkd, fd0, fh0, fkga, fga0];
  initials = {h[0] == h0c, g[0] == g0c, d[0] == d0c, hg[0] == hg0c, hd[0] == hd0c};
  solv = ParametricNDSolveValue[{eqkin, initials}, hd, {t, ftmin, ftmax}, {khd1, khg1, kequhd, kequhg}];
  Return[solv];];

kinSBAParaCacheDall[fkd_, fkga_, fk1hd_, fk1hga_, fh0_, fd0_, fga0_, ftmin_, ftmax_] := kinSBAParaCacheDall[fkd, fkga, fk1hd, fk1hga, fh0, fd0, fga0, ftmin, ftmax] =
Block[{subHD, solv, h0c, d0c, g0c, hg0c, hd0c, initials, stepsize},
{hd0c, h0c, d0c, hg0c, g0c} = SBAStart[fkd, fd0, fh0, fkga, fga0];
initials = {h[0] == h0c, g[0] == g0c, d[0] == d0c, hg[0] == hg0c, hd[0] == hd0c};
solv = ParametricNDSolveValue[{eqkin, initials}, d, {t, ftmin, ftmax}, {khd1, khg1, kequhd, kequhg}];
Return[solv];];

kinSBAParaCacheHall[fkd_, fkga_, fk1hd_, fk1hga_, fh0_, fd0_, fga0_, ftmin_, ftmax_] := kinSBAParaCacheHall[fkd, fkga, fk1hd, fk1hga, fh0, fd0, fga0, ftmin, ftmax] =
Block[{subHD, solv, h0c, d0c, g0c, hg0c, hd0c, initials, stepsize},
{hd0c, h0c, d0c, hg0c, g0c} = SBAStart[fkd, fd0, fh0, fkga, fga0];
initials = {h[0] == h0c, g[0] == g0c, d[0] == d0c, hg[0] == hg0c, hd[0] == hd0c};
solv = ParametricNDSolveValue[{eqkin, initials}, h, {t, ftmin, ftmax}, {khd1, khg1, kequhd, kequhg}];
Return[solv];];

kinSBAParaCacheHGall[fkd_, fkga_, fk1hd_, fk1hga_, fh0_, fd0_, fga0_, ftmin_, ftmax_] := kinSBAParaCacheHGall[fkd, fkga, fk1hd, fk1hga, fh0, fd0, fga0, ftmin, ftmax] =
Block[{subHD, solv, h0c, d0c, g0c, hg0c, hd0c, initials, stepsize},
{hd0c, h0c, d0c, hg0c, g0c} = SBAStart[fkd, fd0, fh0, fkga, fga0];
initials = {h[0] == h0c, g[0] == g0c, d[0] == d0c, hg[0] == hg0c, hd[0] == hd0c};
solv = ParametricNDSolveValue[{eqkin, initials}, hg, {t, ftmin, ftmax}, {khd1, khg1, kequhd, kequhg}];
Return[solv];];

kinSBAParaCacheGall[fkd_, fkga_, fk1hd_, fk1hga_, fh0_, fd0_, fga0_, ftmin_, ftmax_] := kinSBAParaCacheGall[fkd, fkga, fk1hd, fk1hga, fh0, fd0, fga0, ftmin, ftmax] =
Block[{subHD, solv, h0c, d0c, g0c, hg0c, hd0c, initials, stepsize},
{hd0c, h0c, d0c, hg0c, g0c} = SBAStart[fkd, fd0, fh0, fkga, fga0];
initials = {h[0] == h0c, g[0] == g0c, d[0] == d0c, hg[0] == hg0c, hd[0] == hd0c};
solv = ParametricNDSolveValue[{eqkin, initials}, g, {t, ftmin, ftmax}, {khd1, khg1, kequhd, kequhg}];
Return[solv];];

End[] (* End Private Context *)

EndPackage[]
