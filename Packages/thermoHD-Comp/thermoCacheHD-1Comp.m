(* Wolfram Language Package *)

BeginPackage["thermoHDIDA1`"]
(* Exported symbols added here with SymbolName::usage *)
sthdIDAcacheHD::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";

sthdIDAcacheHDresult::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";

sthdIDAcacheHDData::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";

sthdIDAcacheD::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";

sthdIDAcacheDresult::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";

sthdIDAcacheDData::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";

sthdIDAcacheSignal::usage = "st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";

sthdIDAcacheSignallog::usage = "st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";

sthdIDAcacheSignalData::usage = "st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";


sthdiIDAcacheHD::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";

sthdiIDAcacheHDresult::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";

sthdiIDAcacheHDData::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";



sthdiIDAcacheD::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";

sthdiIDAcacheDresult::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";

sthdiIDAcacheDData::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";

sthdiIDAcacheSignal::usage = "";

sthdiIDAcacheSignallog::usage = "";

sthdiIDAcacheSignalData::usage = "";

sthdIDAcacheHGresult::usage = "";
sthdIDAcacheHGData::usage = "";


sthdGDAcacheHDresult::usage = "";
sthdGDAcacheHGresult::usage = "";
sthdGDAcacheDresult::usage = "";
sthdGDAcacheHDData::usage = "";
sthdGDAcacheHGData::usage = "";
sthdGDAcacheDData::usage = "";
sthdGDAcacheSignal::usage = "";
sthdGDAcacheSignalData::usage = "";


Begin["`Private`"] (* Begin Private Context *)

eqthermo = {
   h0 == h + hd + hga,
   d0 == d + hd,
   ga0 == ga + hga,
   kga == hga/(h*ga),
   kd == hd/(h*d)
   };

(* IDA HD*)
		sthdIDAcacheHD[fkd_, fkga_, fh0_, fd0_] :=  sthdIDAcacheHD[fkd, fkga, fh0, fd0] =
		   	Module[{eliHD, solvHD},
		    	eliHD =
		     Eliminate[
		      eqthermo /. {d0 -> fd0, h0 -> fh0, kd -> fkd, kga -> fkga}, {h, d, hga, ga}];
		    	solvHD = hd /. NSolve[eliHD, hd];
		    	Return[solvHD];
		    	];

		sthdIDAcacheHDresult[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ] :=
		  	Module[{sol, e, hs0, gas0, ds0, brac},
		  		Off[LessEqual::nord];
		   		If[fh0 == 0, hs0 = 1.*10^-15, hs0 = fh0];
		   		If[fd0 == 0, ds0 = 1.*10^-15, ds0 = fd0];
		   		If[fga0 == 0, gas0 = 1.*10^-15, gas0 = fga0];
		   		sol = sthdIDAcacheHD[fkd, fkga, hs0, ds0] /. {ga0 -> gas0};
		   		e = Select[Re[sol], 0 <= # <= hs0 && 0 <= # <= ds0 &];
		   		brac = First[e];
		   	Return[brac];];

		sthdIDAcacheHDData[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fstep_?NumericQ] :=

		  	Module[{sol, i},

		   		sol = Table[{N[i], sthdIDAcacheHDresult[fkd, fkga, fh0, fd0, i]}, {i, 0, fga0, fstep}];
		  		Return[sol];
		   ];

(* IDA HG*)
		sthdIDAcacheHG[fkd_, fkga_, fh0_, fd0_] :=  sthdIDAcacheHG[fkd, fkga, fh0, fd0] =
		   	Module[{eliHG, solvHG},
		    	eliHG =
		     Eliminate[
		      eqthermo /. {d0 -> fd0, h0 -> fh0, kd -> fkd, kga -> fkga}, {h, d, hd, ga}];
		    	solvHG = hga /. NSolve[eliHG, hga];
		    	Return[solvHG];
		    	];

		sthdIDAcacheHGresult[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ] :=
		  	Module[{sol, e, hs0, gas0, ds0, brac},
		  		Off[LessEqual::nord];
		   		If[fh0 == 0, hs0 = 1.*10^-15, hs0 = fh0];
		   		If[fd0 == 0, ds0 = 1.*10^-15, ds0 = fd0];
		   		If[fga0 == 0, gas0 = 1.*10^-15, gas0 = fga0];
		   		sol = sthdIDAcacheHG[fkd, fkga, hs0, ds0] /. {ga0 -> gas0};
		   		e = Select[Re[sol], 0 <= # <= hs0 && 0 <= # <= gas0 &];
		   		brac = First[e];
		   	Return[brac];];

		sthdIDAcacheHGData[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fstep_?NumericQ] :=

		  	Module[{sol, i},

		   		sol = Table[{N[i], sthdIDAcacheHGresult[fkd, fkga, fh0, fd0, i]}, {i, 0, fga0, fstep}];
		  		Return[sol];
		   ];


(* IDA Dye*)

		sthdIDAcacheD[fkd_, fkga_, fh0_, fd0_] :=  sthdIDAcacheD[fkd, fkga, fh0, fd0] =
		   	Module[{eliHD, solvHD},
		    	eliHD = Eliminate[eqthermo /. {d0 -> fd0, h0 -> fh0, kd -> fkd, kga -> fkga}, {h, hd, hga, ga}];
		    	solvHD = d /. NSolve[eliHD, d];
		    	Return[solvHD];
		    	];

		sthdIDAcacheDresult[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ] :=
		  	Module[{sol, e, hs0, gas0, ds0, brac},
		  		Off[LessEqual::nord];
		   		If[fh0 == 0, hs0 = 1.*10^-15, hs0 = fh0];
		   		If[fd0 == 0, ds0 = 1.*10^-15, ds0 = fd0];
		   		If[fga0 == 0, gas0 = 1.*10^-15, gas0 = fga0];
		   		sol = sthdIDAcacheD[fkd, fkga, hs0, fd0] /. {ga0 -> gas0};
		   		e = Select[Re[sol], 0 <= # <= ds0 &];
		   		brac = First[e];
		   		Return[brac];
		   	];


		sthdIDAcacheDData[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fstep_?NumericQ] :=

		  	Module[{sol, i},

		   		sol = Table[{N[i], sthdIDAcacheDresult[fkd, fkga, fh0, fd0, i]}, {i, 0, fga0, fstep}];
		  		Return[sol];
		   ];

(* IDA Signal*)

		sthdIDAcacheSignal[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ] :=
			Module[{sol, solsigHD,solsigD},
		   		If[fsigHD != 0, solsigHD = sthdIDAcacheHDresult[fkd, fkga, fh0, fd0, fga0]];
				If[fsigD != 0, solsigD = sthdIDAcacheDresult[fkd, fkga, fh0, fd0, fga0]];
				sol = fsig0 + fsigHD*solsigHD + fsigD*solsigD;
				Return[sol];
		   ];

		sthdIDAcacheSignallog[fpkd_?NumericQ, fpkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ] :=
			Module[{sol, solsigHD,solsigD},
		   		If[fsigHD != 0, solsigHD = sthdIDAcacheHDresult[(10^fpkd), (10^fpkga), fh0, fd0, fga0]];
				If[fsigD != 0, solsigD = sthdIDAcacheDresult[(10^fpkd), (10^fpkga), fh0, fd0, fga0]];
				sol = fsig0 + fsigHD*solsigHD + fsigD*solsigD;
				Return[sol];
		   	];

		sthdIDAcacheSignalData[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ, fstep_?NumericQ] :=
		  	Module[{sol, i},
		   		sol = Table[{N[i], sthdIDAcacheSignal[fkd, fkga, fh0, fd0, i, fsig0, fsigHD, fsigD]}, {i, 0, fga0, fstep}];
		   		Return[sol];
		   	];

(* iIDA HD*)

		sthdiIDAcacheHD[fkd_, fkga_, fga0_] := sthdiIDAcacheHD[fkd, fkga, fga0] =
		   	Module[{eliHD, solvHD},
		   		eliHD = Eliminate[eqthermo /. {ga0 -> fga0, kd -> fkd, kga -> fkga}, {h, d, hga, ga}];
		    	solvHD = hd /. NSolve[eliHD, hd];
		    	Return[solvHD];
		    ];

		sthdiIDAcacheHDresult[fkd_?NumericQ, fkga_?NumericQ, fhd0_?NumericQ, fga0_?NumericQ] :=
		  	Module[{sol, e, hds0, gas0, brac},
		   		Off[LessEqual::nord];
		   		If[fhd0 == 0, hds0 = 1.*10^-15, hds0 = fhd0];
		   		If[fga0 == 0, gas0 = 1.*10^-15, gas0 = fga0];
		   		sol = sthdiIDAcacheHD[fkd, fkga, gas0] /. {h0 -> hds0, d0 -> hds0};
		   		e = Select[Re[sol], 0 <= # <= hds0 &];
		   		brac = First[e];
		   		Return[brac];
		   	];

		sthdiIDAcacheHDData[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fstep_?NumericQ] :=
		  	Module[{sol, i},
		   		sol = Table[{N[i], sthdiIDAcacheHDresult[fkd, fkga, i, i, fga0]}, {i, 0, fh0, fstep}];
		   		Return[sol];
		   	];

(* iIDA Dye*)

		sthdiIDAcacheD[fkd_, fkga_, fga0_] := sthdiIDAcacheD[fkd, fkga, fga0] =
		   	Module[{eliHD, solvHD},
		   		eliHD = Eliminate[eqthermo /. {ga0 -> fga0, kd -> fkd, kga -> fkga}, {h, hd, hga, ga}];
		    	solvHD = d /. NSolve[eliHD, d];
		    	Return[solvHD];
		    ];

		sthdiIDAcacheDresult[fkd_?NumericQ, fkga_?NumericQ, fhd0_?NumericQ, fga0_?NumericQ] :=
		  	Module[{sol, e, hds0, gas0, brac},
		   		Off[LessEqual::nord];
		   		If[fhd0 == 0, hds0 = 1.*10^-15, hds0 = fhd0];
		   		If[fga0 == 0, gas0 = 1.*10^-15, gas0 = fga0];
		   		sol = sthdiIDAcacheD[fkd, fkga, gas0] /. {h0 -> hds0, d0 -> hds0};
		   		e = Select[Re[sol], 0 <= # <= hds0 &];
		   		brac = First[e];
		   		Return[brac];
		   	];

		sthdiIDAcacheDData[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fstep_?NumericQ] :=
		  	Module[{sol, i},
		   		sol = Table[{N[i], sthdiIDAcacheDresult[fkd, fkga, i, i, fga0]}, {i, 0, fh0, fstep}];
		   		Return[sol];
		   	];

(* iIDA Signal*)

		sthdiIDAcacheSignal[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ] :=
			Module[{sol, solsigHD,solsigD},
		   		If[fsigHD != 0, solsigHD = sthdiIDAcacheHDresult[fkd, fkga, fh0, fd0, fga0]];
				If[fsigD != 0, solsigD = sthdiIDAcacheDresult[fkd, fkga, fh0, fd0, fga0]];
				sol = fsig0 + fsigHD*solsigHD + fsigD*solsigD;
				Return[sol];
		   ];

		sthdiIDAcacheSignallog[fpkd_?NumericQ, fpkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ] :=
			Module[{sol, solsigHD,solsigD},
		   		If[fsigHD != 0, solsigHD = sthdiIDAcacheHDresult[(10^fpkd), (10^fpkga), fh0, fd0, fga0]];
				If[fsigD != 0, solsigD = sthdiIDAcacheDresult[(10^fpkd), (10^fpkga), fh0, fd0, fga0]];
				sol = fsig0 + fsigHD*solsigHD + fsigD*solsigD;
				Return[sol];
		   	];

		sthdiIDAcacheSignalData[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ, fstep_?NumericQ] :=
		  	Module[{sol, i},
		   		sol = Table[{N[i], sthdiIDAcacheSignal[fkd, fkga, i, i, fga0, fsig0, fsigHD, fsigD]}, {i, 0, fd0, fstep}];
		   		Return[sol];
		   	];

(* GDA*)
(* Elimination*)
(* GDA HD*)
sthdGDAcacheHD[fkd_, fkga_, fh0_, fga0_] :=  sthdGDAcacheHD[fkd, fkga, fh0, fga0] =
   	Module[{eliHD, solvHD},
    	eliHD =
     Eliminate[
      eqthermo /. {ga0 -> fga0, h0 -> fh0, kd -> fkd, kga -> fkga}, {h, d, hga, ga}];
    	solvHD = hd /. NSolve[eliHD, hd];
    	Return[solvHD];
    	];

(* GDA HG*)
sthdGDAcacheHG[fkd_, fkga_, fh0_, fga0_] :=  sthdGDAcacheHG[fkd, fkga, fh0, fga0] =
   	Module[{eliHG, solvHG},
    	eliHG =
     Eliminate[
      eqthermo /. {ga0 -> fga0, h0 -> fh0, kd -> fkd, kga -> fkga}, {h, d, hd, ga}];
    	solvHG = hga /. NSolve[eliHG, hga];
    	Return[solvHG];
    	];

(* GDA Dye*)

sthdGDAcacheD[fkd_, fkga_, fh0_, fga0_] :=  sthdGDAcacheD[fkd, fkga, fh0, fga0] =
   	Module[{eliHD, solvHD},
    	eliHD = Eliminate[eqthermo /. {ga0 -> fga0, h0 -> fh0, kd -> fkd, kga -> fkga}, {h, hd, hga, ga}];
    	solvHD = d /. NSolve[eliHD, d];
    	Return[solvHD];
    	];

(* Solution*)
(* GDA HD*)
sthdGDAcacheHDresult[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ] := sthdGDAcacheHDresult[fkd, fkga, fh0, fd0, fga0]=
  	Module[{sol, e, hs0, gas0, ds0, brac},
  		Off[LessEqual::nord];
   		If[fh0 == 0, hs0 = 1.*10^-15, hs0 = fh0];
   		If[fd0 == 0, ds0 = 1.*10^-15, ds0 = fd0];
   		If[fga0 == 0, gas0 = 1.*10^-15, gas0 = fga0];
   		sol = sthdGDAcacheHD[fkd, fkga, hs0, gas0] /. {d0 -> ds0};
   		e = Select[Re[sol], 0 <= # <= hs0 && 0 <= # <= ds0 &];
   		brac = First[e];
   	Return[brac];];

(* GDA HG*)
sthdGDAcacheHGresult[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ] :=sthdGDAcacheHGresult[fkd, fkga, fh0, fd0, fga0]=
  	Module[{sol, e, hs0, gas0, ds0, brac},
  		Off[LessEqual::nord];
   		If[fh0 == 0, hs0 = 1.*10^-15, hs0 = fh0];
   		If[fd0 == 0, ds0 = 1.*10^-15, ds0 = fd0];
   		If[fga0 == 0, gas0 = 1.*10^-15, gas0 = fga0];
   		sol = sthdGDAcacheHG[fkd, fkga, hs0, gas0] /. {d0 -> ds0};
   		e = Select[Re[sol], 0 <= # <= hs0 && 0 <= # <= gas0 &];
   		brac = First[e];
   	Return[brac];];

(* GDA D*)
sthdGDAcacheDresult[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ] :=sthdGDAcacheDresult[fkd, fkga, fh0, fd0, fga0]=
  	Module[{sol, e, hs0, gas0, ds0, brac},
  		Off[LessEqual::nord];
   		If[fh0 == 0, hs0 = 1.*10^-15, hs0 = fh0];
   		If[fd0 == 0, ds0 = 1.*10^-15, ds0 = fd0];
   		If[fga0 == 0, gas0 = 1.*10^-15, gas0 = fga0];
   		sol = sthdGDAcacheD[fkd, fkga, hs0, gas0] /. {d0 -> ds0};
   		e = Select[Re[sol], 0 <= # <= ds0 &];
   		brac = First[e];
   		Return[brac];
   	];

(* Data *)
sthdGDAcacheHDData[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fstep_?NumericQ] :=

  	Module[{sol, i},

   		sol = Table[{N[i], sthdGDAcacheHDresult[fkd, fkga, fh0, i, fga0]}, {i, 0, fd0, fstep}];
  		Return[sol];
   ];


sthdGDAcacheHGData[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fstep_?NumericQ] :=

  	Module[{sol, i},

   		sol = Table[{N[i], sthdGDAcacheHGresult[fkd, fkga, fh0, i, fga0]}, {i, 0, fd0, fstep}];
  		Return[sol];
   ];

sthdGDAcacheDData[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fstep_?NumericQ] :=

  	Module[{sol, i},

   		sol = Table[{N[i], sthdGDAcacheDresult[fkd, fkga, fh0, i, fga0]}, {i, 0, fd0, fstep}];
  		Return[sol];
   ];

(* GDA Signal*)

sthdGDAcacheSignal[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ] :=
	Module[{sol, solsigHD,solsigD},
   		If[fsigHD != 0, solsigHD = sthdGDAcacheHDresult[fkd, fkga, fh0, fd0, fga0]];
		If[fsigD != 0, solsigD = sthdGDAcacheDresult[fkd, fkga, fh0, fd0, fga0]];
		sol = fsig0 + fsigHD*solsigHD + fsigD*solsigD;
		Return[sol];
   ];

sthdGDAcacheSignallog[fpkd_?NumericQ, fpkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ] :=
	Module[{sol, solsigHD,solsigD},
   		If[fsigHD != 0, solsigHD = sthdGDAcacheHDresult[(10^fpkd), (10^fpkga), fh0, fd0, fga0]];
		If[fsigD != 0, solsigD = sthdGDAcacheDresult[(10^fpkd), (10^fpkga), fh0, fd0, fga0]];
		sol = fsig0 + fsigHD*solsigHD + fsigD*solsigD;
		Return[sol];
   	];

sthdGDAcacheSignalData[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ, fstep_?NumericQ] :=
  	Module[{sol, i},
   		sol = Table[{N[i], sthdGDAcacheSignal[fkd, fkga, fh0, i, fga0, fsig0, fsigHD, fsigD]}, {i, 0, fd0, fstep}];
   		Return[sol];
   	];


End[] (* End Private Context *)

EndPackage[]
