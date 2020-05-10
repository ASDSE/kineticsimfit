(* Wolfram Language Package *)

BeginPackage["thermoHD`"]
(* Exported symbols added here with SymbolName::usage *)
sthdCacheSignal::usage = "
	st-hg-Signal[kd,d0,h0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";

sthdCacheSignallog::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";

sthdCacheResultHD::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";

sthdCacheHDData::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";

sthdCacheResultH::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";

sthdCacheResultD::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";


sthdCacheSignalData::usage = "
	st-hg-Signal[[kd,kga,h0,d0,ga0,sig0,sigHD,sigD,Ga0steps]
	Creates thermodynamic Signal data for the HD2 + Ga -> HGa + 2D binding scenario";

sthdCacheHDKdD0::usage="";

sthdCacheHKdD0::usage="";

sthdCacheDKdD0::usage="";

Begin["`Private`"] (* Begin Private Context *)

eqthermo = {
   xh0 == xh + xhd,
   xd0 == xd + xhd,
   kd == xhd/(xh*xd)
   };



sthdCacheHDKdD0[fkd_, fd0_] := sthdCacheHDKdD0[fkd, fd0]=
	Module[{eliHD, solvHD},
  	eliHD = Eliminate[eqthermo /. {d0 -> fd0, kd -> fkd}, {h, d}];
  	solvHD = hd /. NSolve[eliHD, hd];
  	Return[solvHD];
  	];

sthdCacheResultHD[fkd_?NumericQ, fd0_?NumericQ, fh0_?NumericQ] :=
	Module[{sol, e, hs0, ds0, brac},
  		If[fh0 == 0, hs0 = 1.*10^-15, hs0 = fh0];
  		If[fd0 == 0, ds0 = 1.*10^-15, ds0 = fd0];
  		sol = sthdCacheHDKdD0[fkd, ds0] /. {h0 -> hs0};
  		e = Select[sol, 0 <= # <= hs0 && 0 <= # <= ds0 &];
  		brac = First[e];
  		Return[brac];
	];


sthdCacheHKdD0[fkd_, fd0_] := sthdCacheHKdD0[fkd, fd0]=
	Module[{eliHD, solvHD},
  	eliHD = Eliminate[eqthermo /. {d0 -> fd0, kd -> fkd}, {hd, d}];
  	solvHD = h /. NSolve[eliHD, h];
  	Return[solvHD];
  	];

sthdCacheResultH[fkd_?NumericQ, fd0_?NumericQ, fh0_?NumericQ] :=
	Module[{sol, e, hs0, ds0, brac},
  		If[fh0 == 0, hs0 = 1.*10^-15, hs0 = fh0];
  		If[fd0 == 0, ds0 = 1.*10^-15, ds0 = fd0];
  		sol = sthdCacheHKdD0[fkd, ds0] /. {h0 -> hs0};
  		e = Select[sol, 0 <= # <= hs0 &];
  		brac = First[e];
  		Return[brac];
	];


sthdCacheDKdD0[fkd_, fd0_] := sthdCacheDKdD0[fkd, fd0]=
	Module[{eliHD, solvHD},
  	eliHD = Eliminate[eqthermo /. {d0 -> fd0, kd -> fkd}, {hd, h}];
  	solvHD = d /. NSolve[eliHD, d];
  	Return[solvHD];
  	];

sthdCacheResultD[fkd_?NumericQ, fd0_?NumericQ, fh0_?NumericQ] :=
	Module[{sol, e, hs0, ds0, brac},
  		If[fh0 == 0, hs0 = 1.*10^-15, hs0 = fh0];
  		If[fd0 == 0, ds0 = 1.*10^-15, ds0 = fd0];
  		sol = sthdCacheDKdD0[fkd, ds0] /. {h0 -> hs0};
  		e = Select[sol, 0 <= # <= ds0 &];
  		brac = First[e];
  		Return[brac];
	];




sthdCacheSignal[fkd_?NumericQ, fd0_?NumericQ, fh0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ] :=
	Module[{sol},
   		sol = fsig0 + fsigHD*sthdCacheResultHD[fkd, fd0, fh0] + fsigD*sthdCacheResultD[fkd, fd0, fh0]
   ];

sthdCacheSignallog[fpkd_?NumericQ, fd0_?NumericQ, fh0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ] :=
	Module[{sol},
   		sol = fsig0 + fsigHD*sthdCacheResultHD[(10^fpkd), fd0, fh0] + fsigD*sthdCacheResultD[(10^fpkd), fd0, fh0]
   ];

sthdCacheSignalData[fkd_?NumericQ,  fd0_?NumericQ, fh0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ, fstep_?NumericQ] :=
  	Module[{sol, i},
   		sol = Table[{N[i], sthdCacheSignal[fkd,  fd0, i, fsig0, fsigHD, fsigD]}, {i, 0, fh0, fstep}];
   		Return[sol];
   	];

sthdCacheHDData[fkd_?NumericQ,  fd0_?NumericQ, fh0_?NumericQ, fstep_?NumericQ] :=
  	Module[{sol, i},
   		sol = Table[{N[i], sthdCacheResultHD[fkd, fd0, i]}, {i, 0, fh0, fstep}];
   		Return[sol];
   	];

sthdCacheDData[fkd_?NumericQ,  fd0_?NumericQ, fh0_?NumericQ, fstep_?NumericQ] :=
  	Module[{sol, i},
   		sol = Table[{N[i],  sthdCacheResultD[fkd, fd0, i]}, {i, 0, fh0, fstep}];
   		Return[sol];
   	];


sthdCacheHData[fkd_?NumericQ,  fd0_?NumericQ, fh0_?NumericQ, fstep_?NumericQ] :=
  	Module[{sol, i},
   		sol = Table[{N[i],  sthdCacheResultH[fkd, fd0, i]}, {i, 0, fh0, fstep}];
   		Return[sol];
   	];

End[] (* End Private Context *)

EndPackage[]
