(* Wolfram Language Package *)

BeginPackage["kinHDCompStart`"]
(* Exported symbols added here with SymbolName::usage *)


IDAStart::usage = "";
GDAStart::usage = "";
SBAStart::usage = "";
sthdCacheResultHD::usage = "";


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


IDAParameters[fkd_?NumericQ, fd0_?NumericQ, fh0_?NumericQ, fkga0_?NumericQ, fga0_?NumericQ] := IDAStart[fkd, fd0, fh0, fkga0, fga0] =
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


End[] (* End Private Context *)

EndPackage[]
