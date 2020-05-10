(* Wolfram Language Package *)

BeginPackage["kintothermo`"]
(* Exported symbols added here with SymbolName::usage *)
kinIDAHDdata::usage = "";
Begin["`Private`"] (* Begin Private Context *)

eqthermo = {h0 == h + hd, d0 == d + hd, kd == hd/(h*d)};

eqkin = {hd'[t] == (-khd1/kequhd)*hd[t] + khd1*h[t]*d[t],
   hg'[t] == (khg1*h[t]*g[t]) - (khg1/kequhg*hg[t]),
   h'[t] == (khg1/kequhg*hg[t]) - (khg1*h[t]*g[t]) + (khd1/kequhd)*
      hd[t] - khd1*h[t]*d[t],
   g'[t] == ((khg1/kequhg)*hg[t]) - (khg1*h[t]*g[t]),
   d'[t] == khd1/kequhd*hd[t] - khd1*h[t]*d[t]};

sthdCacheHDKdD0[fkd_?NumericQ, fd0_?NumericQ] :=
  sthdCacheHDKdD0[fkd, fd0] =
   Block[{eli, solv},
    eli = Eliminate[eqthermo /. {d0 -> fd0, kd -> fkd}, {h, d}];
    solv = hd /. NSolve[eli, hd];
    Return[solv];];

sthdCacheResultHD[fkd_?NumericQ, fd0_?NumericQ, fh0_?NumericQ] :=
  sthdCacheResultHD[fkd, fd0, fh0] =
   Block[{sol, e, hs0, ds0, brac},
    If[fh0 == 0, hs0 = 1.*10^-15, hs0 = fh0];
    If[fd0 == 0, ds0 = 1.*10^-15, ds0 = fd0];
    sol = sthdCacheHDKdD0[fkd, ds0] /. {h0 -> hs0};
    e = Select[sol, 0 <= # <= hs0 && 0 <= # <= ds0 &];
    brac = First[e];
    Return[brac];];

IDAStart[fkd_?NumericQ, fd0_?NumericQ, fh0_?NumericQ, fkga0_?NumericQ,
    fga0_?NumericQ] :=
  IDAStart[fkd, fd0, fh0, fkga0, fga0] =
   Block[{h0c, d0c, g0c, hg0c, hd0c, initials, stepsize, sol},
    hd0c = sthdCacheResultHD[fkd, fd0, fh0];
    h0c = fh0 - hd0c;
    d0c = fd0 - hd0c;
    hg0c = 0;
    g0c = fga0;
    sol = {hd0c, h0c, d0c, hg0c, g0c};
    Return[N[sol]]];

parametricHDfun[finitials_, ftmin_, ftmax_] :=
  parametricHDfun[finitials, ftmin, ftmax] =
   Block[{solv},
    solv =
     ParametricNDSolveValue[{eqkin, finitials},
      hd, {t, ftmin, ftmax}, {khd1, khg1, kequhd, kequhg}];
    Return[solv];
    ];

parametricHGfun[finitials_, ftmin_, ftmax_] :=
  parametricHGfun[finitials, ftmin, ftmax] =
   Block[{solv},
    solv =
     ParametricNDSolveValue[{eqkin, finitials},
      hg, {t, ftmin, ftmax}, {khd1, khg1, kequhd, kequhg}];
    Return[solv];
    ];

kinIDAHDdata[fkd_, fkga_, fk1hd_, fk1hga_, fh0_, fd0_, fga0_, ftmin_, ftmax_, ftimesteps_, fconcsteps_] :=
  kinIDAHDdata[fkd, fkga, fk1hd, fk1hga, fh0, fd0, fga0, ftmin, ftmax, ftimesteps, fconcsteps] =
    Block[{subHD, solv, h0c, d0c, g0c, hg0c, gc, g0added, hd0c, initials, concsteps, dataHDkin, ga0s, c, timedata, hgdata, summary, list, hd0cache, dataHDtherm, concdata},
      {hd0c, h0c, d0c, hg0c, g0c} = IDAStart[fkd, fd0, fh0, fkga, fga0];
      initials = {h[0] == h0c, g[0] == g0c, d[0] == d0c, hg[0] == hg0c, hd[0] == hd0c};
      summary = {hd0c, h0c, d0c, hg0c, g0c};
      gc = 0;
      g0c = 0;
      dataHDkin = {};
      dataHDtherm = {{g0c, hd0c}};

      For[c = 1., c <= fconcsteps, c++,
          g0added = (fga0/fconcsteps)*c;
          g0c = (fga0/fconcsteps) + gc;
          timedata =
          Table[
          {x + ftmax*(c - 1),
          parametricHDfun[{h[0] == h0c, g[0] == g0c, d[0] == d0c, hg[0] == hg0c, hd[0] == hd0c}, ftmin, ftmax][fk1hd, fk1hga, fkd, fkga][x]},
          {x, 0., ftmax, ftmax/ftimesteps}
          ];

          dataHDkin = Join[dataHDkin, timedata];
          hd0cache = hd0c;
          hd0c = parametricHDfun[{h[0] == h0c, g[0] == g0c, d[0] == d0c, hg[0] == hg0c, hd[0] == hd0c}, ftmin, ftmax][fk1hd, fk1hga, fkd, fkga][ftmax];

          concdata = {{g0added, hd0c}};
          dataHDtherm = Join[dataHDtherm, concdata];
          hg0c = parametricHGfun[{h[0] == h0c, g[0] == g0c, d[0] == d0c, hg[0] == hg0c, hd[0] == hd0cache}, ftmin, ftmax][fk1hd, fk1hga, fkd, fkga][ftmax];
          gc = (fga0/fconcsteps)*(c + 1) - hg0c;
          d0c = fd0 - hd0c;
          h0c = fh0 - (hg0c + hd0c);
       ];


list = {dataHDtherm, dataHDkin};

Return[list];

];

End[] (* End Private Context *)

EndPackage[]
