(* Wolfram Language Package *)

BeginPackage["thermoHD`"]
(* Exported symbols added here with SymbolName::usage *)  
sthdSignal::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]	
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";
	
sthdSignallog::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]	
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";
	
sthdHD::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]	
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";
	
sthdHDData::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]	
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";
	
sthdDData::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]	
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";

sthdHData::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]	
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";
	
		
sthdSignalData::usage = "
	st-hg-Signal[[kd,kga,h0,d0,ga0,sig0,sigHD,sigD,Ga0steps]	
	Creates thermodynamic Signal data for the HD2 + Ga -> HGa + 2D binding scenario";


Begin["`Private`"] (* Begin Private Context *) 

eqthermo = {
   h0 == h + hd,
   d0 == d + hd,
   kd == hd/(h*d)
   };
   
furearrangedHD[fkd_, fh0_, fd0_] := 
 	 Eliminate[eqthermo /. {h0 -> fh0, d0 -> fd0, kd -> fkd}, {h, d}
    	];
    	
sthdHD[fkd_?NumericQ, fh0_?NumericQ, fd0_?NumericQ] := 
	Module[{sol, e, hs0, brac},
   		If[fh0 == 0, hs0 = 1*10^-15, hs0 = fh0];
   		sol = hd /. NSolve[furearrangedHD[fkd, hs0, fd0], hd, Reals];
   		e = Select[sol, 0 <= # <= hs0 && 0 <= # <= fd0 &];
   		brac=First[e];
   		Return[brac]
   	];
   	
furearrangedH[fkd_, fh0_, fd0_] :=
		Eliminate[eqthermo /. {h0 -> fh0, d0 -> fd0, kd -> fkd}, {hd, d}
    	];
    	
sthdH[fkd_?NumericQ, fh0_?NumericQ, fd0_?NumericQ] := 
	Module[{sol, e, hs0, brac},
   		If[fh0 == 0, hs0 = 1*10^-15, hs0 = fh0];
   		sol = h /. NSolve[furearrangedH[fkd, hs0, fd0], h, Reals];
   		e = Select[sol, 0 <= # <= hs0 &];
   		brac=First[e];
   		Return[e]
   	];

furearrangedD[fkd_, fh0_, fd0_] :=
		Eliminate[eqthermo /. {h0 -> fh0, d0 -> fd0, kd -> fkd}, {h, hd}
     ];
     
sthdD[fkd_?NumericQ, fh0_?NumericQ, fd0_?NumericQ] := 
	Module[{sol, e, hs0, brac},
   		If[fh0 == 0, hs0 = 1*10^-15, hs0 = fh0];
   		sol = d /. NSolve[furearrangedD[fkd, hs0, fd0], d, Reals];
   		e = Select[sol, 0 <= # <= fd0 &];
   		brac=First[e];
   		Return[brac]
  	];



sthdSignal[fkd_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ] := 
	Module[{sol, hs0},
   		If[fh0 == 0, hs0 = 1*10^-15, hs0 = fh0];
   		sol = fsig0 + fsigHD*sthdHD[fkd, hs0, fd0] + fsigD*sthdD[fkd, hs0, fd0]
   ];
   
sthdSignallog[fpkd_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ] := 
	Module[{sol, hs0},
   		If[fh0 == 0, hs0 = 1*10^-15, hs0 = fh0];
   		sol = fsig0 + fsigHD*sthdHD[(10^fpkd), hs0, fd0] + fsigD*sthdD[(10^fpkd), hs0, fd0]
   ];

sthdSignalData[fkd_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ, fstep_?NumericQ] := 
  	Module[{sol, i}, 
   		sol = Table[{N[i], sthdSignal[fkd, i, fd0, fsig0, fsigHD, fsigD]}, {i, 0, fh0, fstep}];
   		Return[sol];
   	];

sthdHDData[fkd_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fstep_?NumericQ] := 
  	Module[{sol, i}, 
   		sol = Table[{N[i], sthdHD[fkd, i, fd0]}, {i, 0, fh0, fstep}];
   		Return[sol];
   	];
   	
sthdDData[fkd_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fstep_?NumericQ] := 
  	Module[{sol, i}, 
   		sol = Table[{N[i], sthdD[fkd, i, fd0]}, {i, 0, fh0, fstep}];
   		Return[sol];
   	];   	
   	
   	
sthdHData[fkd_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fstep_?NumericQ] := 
  	Module[{sol, i}, 
   		sol = Table[{N[i], sthdH[fkd, i, fd0]}, {i, 0, fh0, fstep}];
   		Return[sol];
   	];    	

End[] (* End Private Context *)

EndPackage[]