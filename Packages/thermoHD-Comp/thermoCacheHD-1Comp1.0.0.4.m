(* Wolfram Language Package *)

BeginPackage["thermoHDIDA1`"]
(* Exported symbols added here with SymbolName::usage *)
sthdSBAcacheSignal::usage="";
sthdSBAcacheSignalData::usage="";
sthdSBAcacheHDData::usage="";
sthdSBAcacheDData::usage="";
sthdSBAcacheHDresult::usage="";
sthdSBAcacheDresult::usage="";

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


sthdiADAcacheSignal::usage = "";

sthdiADAcacheSignallog::usage = "";

sthdiADAcacheSignalData::usage = "";

sthdADAcacheSignal::usage = "";

sthdADAcacheSignallog::usage = "";

sthdADAcacheSignalData::usage = "";


sthdADAcacheHDData::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";
sthdADAcacheDData::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";


volumefunctionIDA::usage = "
	st-hg-Signal[kd,kga,h0,d0,ga0,sig0,sigHD,sigD]
	Gives numeric solution for Host in a HD2 + Ga -> HGa + 2D competition scenario";


Begin["`Private`"] (* Begin Private Context *)

eqthermo = {
   h0 == h + hd + hga,
   d0 == d + hd,
   ga0 == ga + hga,
   kga == hga/(h*ga),
   kd == hd/(h*d)
   };


   volumefunctionIDA[fdstart_, fhstart_,fgastart_, fVstart_, fVinj_] :=
  	Module[{sol, d0, h0, ga0},
   	d0 = fVstart*fdstart/(fVstart + fVinj); (*Cell*)
   	h0 = fVstart*fhstart/(fVstart + fVinj); (*Cell*)
   	ga0 = fVinj*fgastart/(fVstart + fVinj); (*Titrant*)
   	sol = {d0, h0, ga0};
   	Return[sol];
   	];


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
		   		sol = sthdIDAcacheHD[fkd, fkga, hs0, fd0] /. {ga0 -> gas0};
		   		e = Select[Re[sol], 0 <= # <= hs0 && 0 <= # <= ds0 &];
		   		brac = First[e];
		   	Return[brac];];

		sthdIDAcacheHDData[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fstep_?NumericQ] :=

		  	Module[{sol, i},

		   		sol = Table[{N[i], sthdIDAcacheHDresult[fkd, fkga, fh0, fd0, i]}, {i, 0, fga0, fstep}];
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

		sthdiIDAcacheHDData[fkd_?NumericQ, fkga_?NumericQ, fhd0_?NumericQ, fga0_?NumericQ, fstep_?NumericQ] :=
		  	Module[{sol, i},
		   		sol = Table[{N[i], sthdiIDAcacheHDresult[fkd, fkga, i, fga0]}, {i, 0, fhd0, fstep}];
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

		sthdiIDAcacheDData[fkd_?NumericQ, fkga_?NumericQ, fhd0_?NumericQ, fga0_?NumericQ, fstep_?NumericQ] :=
		  	Module[{sol, i},
		   		sol = Table[{N[i], sthdiIDAcacheDresult[fkd, fkga, i, fga0]}, {i, 0, fhd0, fstep}];
		   		Return[sol];
		   	];

(* iIDA Signal*)

		sthdiIDAcacheSignal[fkd_?NumericQ, fkga_?NumericQ, fhd0_?NumericQ, fga0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ] :=
			Module[{sol, solsigHD,solsigD},
				If[fhd0==0,sol=fsig0,
		   			If[fsigHD != 0, solsigHD = sthdiIDAcacheHDresult[fkd, fkga, fhd0, fga0]];
					If[fsigD != 0, solsigD = sthdiIDAcacheDresult[fkd, fkga, fhd0, fga0]];
					sol = fsig0 + fsigHD*solsigHD + fsigD*solsigD;
				];
				Return[sol];
		   ];

		sthdiIDAcacheSignallog[fpkd_?NumericQ, fpkga_?NumericQ, fhd0_?NumericQ, fga0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ] :=
			Module[{sol, solsigHD,solsigD},
		   		If[fsigHD != 0, solsigHD = sthdiIDAcacheHDresult[(10^fpkd), (10^fpkga), fhd0, fga0]];
				If[fsigD != 0, solsigD = sthdiIDAcacheDresult[(10^fpkd), (10^fpkga), fhd0, fga0]];
				sol = fsig0 + fsigHD*solsigHD + fsigD*solsigD;
				Return[sol];
		   	];

		sthdiIDAcacheSignalData[fkd_?NumericQ, fkga_?NumericQ, fhd0_?NumericQ, fga0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ, fstep_?NumericQ] :=
		  	Module[{sol, i},
		   		sol = Table[{N[i], sthdiIDAcacheSignal[fkd, fkga, i, fga0, fsig0, fsigHD, fsigD]}, {i, 0, fhd0, fstep}];
		   		Return[sol];
		   	];

(* ADA HD*)
		sthdADAcacheHD[fkd_, fkga_, fh0_, fga0_] :=  sthdADAcacheHD[fkd, fkga, fh0, fga0] =
		   	Module[{eliHD, solvHD},
		    	eliHD =
		     Eliminate[
		      eqthermo /. {ga0 -> fga0, h0 -> fh0, kd -> fkd, kga -> fkga}, {h, d, hga, ga}];
		    	solvHD = hd /. NSolve[eliHD, hd];
		    	Return[solvHD];
		    	];

		sthdADAcacheHDresult[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ] :=
		  	Module[{sol, e, hs0, gas0, ds0, brac},
		  		Off[LessEqual::nord];
		   		If[fh0 == 0, hs0 = 1.*10^-15, hs0 = fh0];
		   		If[fd0 == 0, ds0 = 1.*10^-15, ds0 = fd0];
		   		If[fga0 == 0, gas0 = 1.*10^-15, gas0 = fga0];
		   		sol = sthdADAcacheHD[fkd, fkga, hs0, gas0] /. {d0 -> ds0};
		   		e = Select[Re[sol], 0 <= # <= hs0 && 0 <= # <= ds0 &];
		   		brac = First[e];
		   	Return[brac];];

		sthdADAcacheHDData[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fstep_?NumericQ] :=

		  	Module[{sol, i},

		   		sol = Table[{N[i], sthdADAcacheHDresult[fkd, fkga, fh0, i, fga0]}, {i, 0, fd0, fstep}];
		  		Return[sol];
		   ];


(* ADA Dye*)

		sthdADAcacheD[fkd_, fkga_, fh0_, fga0_] :=  sthdADAcacheD[fkd, fkga, fh0, fga0] =
		   	Module[{eliHD, solvHD},
		    	eliHD = Eliminate[eqthermo /. {ga0 -> fga0, h0 -> fh0, kd -> fkd, kga -> fkga}, {h, hd, hga, ga}];
		    	solvHD = d /. NSolve[eliHD, d];
		    	Return[solvHD];
		    	];

		sthdADAcacheDresult[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ] :=
		  	Module[{sol, e, hs0, gas0, ds0, brac},
		  		Off[LessEqual::nord];
		   		If[fh0 == 0, hs0 = 1.*10^-15, hs0 = fh0];
		   		If[fd0 == 0, ds0 = 1.*10^-15, ds0 = fd0];
		   		If[fga0 == 0, gas0 = 1.*10^-15, gas0 = fga0];
		   		sol = sthdADAcacheD[fkd, fkga, hs0, gas0] /. {d0 -> ds0};
		   		e = Select[Re[sol], 0 <= # <= ds0 &];
		   		brac = First[e];
		   		Return[brac];
		   	];


		sthdADAcacheDData[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fstep_?NumericQ] :=

		  	Module[{sol, i},

		   		sol = Table[{N[i], sthdADAcacheDresult[fkd, fkga, fh0, i, fga0]}, {i, 0, fd0, fstep}];
		  		Return[sol];
		   ];

(* ADA Signal*)

		sthdADAcacheSignal[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ] :=
			Module[{sol, solsigHD,solsigD},
		   		If[fsigHD != 0, solsigHD = sthdADAcacheHDresult[fkd, fkga, fh0, fd0, fga0]];
				If[fsigD != 0, solsigD = sthdADAcacheDresult[fkd, fkga, fh0, fd0, fga0]];
				sol = fsig0 + fsigHD*solsigHD + fsigD*solsigD;
				Return[sol];
		   ];

		sthdADAcacheSignallog[fpkd_?NumericQ, fpkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ] :=
			Module[{sol, solsigHD,solsigD},
		   		If[fsigHD != 0, solsigHD = sthdADAcacheHDresult[(10^fpkd), (10^fpkga), fh0, fd0, fga0]];
				If[fsigD != 0, solsigD = sthdADAcacheDresult[(10^fpkd), (10^fpkga), fh0, fd0, fga0]];
				sol = fsig0 + fsigHD*solsigHD + fsigD*solsigD;
				Return[sol];
		   	];

		sthdADAcacheSignalData[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ, fstep_?NumericQ] :=
		  	Module[{sol, i},
		   		sol = Table[{N[i], sthdADAcacheSignal[fkd, fkga, fh0, i, fga0, fsig0, fsigHD, fsigD]}, {i, 0, fd0, fstep}];
		   		Return[sol];
		   	];

(* iADA HD*)

		sthdiADAcacheHD[fkd_, fkga_, fd0_] := sthdiADAcacheHD[fkd, fkga, fd0] =
		   	Module[{eliHD, solvHD},
		   		eliHD = Eliminate[eqthermo /. {d0 -> fd0, kd -> fkd, kga -> fkga}, {h, d, hga, ga}];
		    	solvHD = hd /. NSolve[eliHD, hd];
		    	Return[solvHD];
		    ];

		sthdiADAcacheHDresult[fkd_?NumericQ, fkga_?NumericQ, fhga0_?NumericQ, fd0_?NumericQ] :=
		  	Module[{sol, e, hgas0, ds0, brac},
		   		Off[LessEqual::nord];
		   		If[fhga0 == 0, hgas0 = 1.*10^-15, hgas0 = fhga0];
		   		If[fd0 == 0, ds0 = 1.*10^-15, ds0 = fd0];
		   		sol = sthdiADAcacheHD[fkd, fkga, ds0] /. {h0 -> hgas0, ga0 -> hgas0};
		   		e = Select[Re[sol], 0 <= # <= hgas0 && 0 <= # <= ds0 &];
		   		brac = First[e];
		   		Return[brac];
		   	];

		sthdiADAcacheHDData[fkd_?NumericQ, fkga_?NumericQ, fhga0_?NumericQ, fd0_?NumericQ, fstep_?NumericQ] :=
		  	Module[{sol, i},
		   		sol = Table[{N[i], sthdiADAcacheHDresult[fkd, fkga, i, fd0]}, {i, 0, fhga0, fstep}];
		   		Return[sol];
		   	];

(* iADA Dye*)

		sthdiADAcacheD[fkd_, fkga_, fd0_] := sthdiADAcacheD[fkd, fkga, fd0] =
		   	Module[{eliHD, solvHD},
		   		eliHD = Eliminate[eqthermo /. {d0 -> fd0, kd -> fkd, kga -> fkga}, {h, hd, hga, ga}];
		    	solvHD = d /. NSolve[eliHD, d];
		    	Return[solvHD];
		    ];

		sthdiADAcacheDresult[fkd_?NumericQ, fkga_?NumericQ, fhga0_?NumericQ, fd0_?NumericQ] :=
		  	Module[{sol, e, hgas0, ds0, brac},
		   		Off[LessEqual::nord];
		   		If[fhga0 == 0, hgas0 = 1.*10^-15, hgas0 = fhga0];
		   		If[fd0 == 0, ds0 = 1.*10^-15, ds0 = fd0];
		   		sol = sthdiADAcacheD[fkd, fkga, ds0] /. {h0 -> hgas0, ga0 -> hgas0};
		   		e = Select[Re[sol], 0 <= # <= hgas0 && 0 <= # <= ds0 &];
		   		brac = First[e];
		   		Return[brac];
		   	];

		sthdiADAcacheDData[fkd_?NumericQ, fkga_?NumericQ, fhga0_?NumericQ, fd0_?NumericQ, fstep_?NumericQ] :=
		  	Module[{sol, i},
		   		sol = Table[{N[i], sthdiADAcacheDresult[fkd, fkga, i, fd0]}, {i, 0, fhga0, fstep}];
		   		Return[sol];
		   	];

(* iADA Signal*)

		sthdiADAcacheSignal[fkd_?NumericQ, fkga_?NumericQ, fhga0_?NumericQ, fd0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ] :=
			Module[{sol, solsigHD,solsigD},
				If[fhga0==0,
					solsigHD = 0;
					solsigD = fd0;,
		   				If[fsigHD != 0, solsigHD = sthdiADAcacheHDresult[fkd, fkga, fhga0, fd0]];
						If[fsigD != 0, solsigD = sthdiADAcacheDresult[fkd, fkga, fhga0, fd0]];
				];
				sol = fsig0 + fsigHD*solsigHD + fsigD*solsigD;
				Return[sol];
		   ];

		sthdiADAcacheSignallog[fpkd_?NumericQ, fpkga_?NumericQ, fhga0_?NumericQ, fd0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ] :=
			Module[{sol, solsigHD,solsigD},
		   		If[fsigHD != 0, solsigHD = sthdiADAcacheHDresult[(10^fpkd), (10^fpkga), fhga0, fd0]];
				If[fsigD != 0, solsigD = sthdiADAcacheDresult[(10^fpkd), (10^fpkga), fhga0, fd0]];
				sol = fsig0 + fsigHD*solsigHD + fsigD*solsigD;
				Return[sol];
		   	];

		sthdiADAcacheSignalData[fkd_?NumericQ, fkga_?NumericQ, fhga0_?NumericQ, fd0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ, fstep_?NumericQ] :=
		  	Module[{sol, i},
		   		sol = Table[{N[i], sthdiADAcacheSignal[fkd, fkga, i, fd0, fsig0, fsigHD, fsigD]}, {i, 0, fhga0, fstep}];
		   		Return[sol];
		   	];



(* SBA HD*)
	sthdSBAcacheHD[fkd_, fkga_, fd0_, fga0_] :=  sthdSBAcacheHD[fkd, fkga, fd0, fga0] =
	   	Module[{eliHD, solvHD},
	    	eliHD =
	     Eliminate[
	      eqthermo /. {ga0 -> fga0, d0 -> fd0, kd -> fkd, kga -> fkga}, {h, d, hga, ga}];
	    	solvHD = hd /. NSolve[eliHD, hd];
	    	Return[solvHD];
	    	];

	sthdSBAcacheHDresult[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ] := sthdSBAcacheHDresult[fkd, fkga, fh0, fd0, fga0] =
	  	Module[{sol, e, hs0, gas0, ds0, brac},
	  		Off[LessEqual::nord];
	   		If[fh0 == 0, hs0 = 1.*10^-15, hs0 = fh0];
	   		If[fd0 == 0, ds0 = 1.*10^-15, ds0 = fd0];
	   		If[fga0 == 0, gas0 = 1.*10^-15, gas0 = fga0];
	   		sol = sthdSBAcacheHD[fkd, fkga, ds0, gas0] /. {h0 -> hs0};
	   		e = Select[Re[sol], 0 <= # <= hs0 && 0 <= # <= ds0 &];
	   		brac = First[e];
	   	Return[brac];];

	sthdSBAcacheHDData[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fstep_?NumericQ] := sthdSBAcacheHDData[fkd, fkga, fh0, fd0, fga0, fstep]=

	  	Module[{sol, i},

	   		sol = Table[{N[i], sthdSBAcacheHDresult[fkd, fkga, i, fd0, fga0]}, {i, 0, fh0, fstep}];
	  		Return[sol];
	   ];


(* SBA Dye*)

	sthdSBAcacheD[fkd_, fkga_, fd0_, fga0_] :=  sthdSBAcacheD[fkd, fkga, fd0, fga0] =
	   	Module[{eliHD, solvHD},
	    	eliHD = Eliminate[eqthermo /. {ga0 -> fga0, d0 -> fd0, kd -> fkd, kga -> fkga}, {h, hd, hga, ga}];
	    	solvHD = d /. NSolve[eliHD, d];
	    	Return[solvHD];
	    	];

	sthdSBAcacheDresult[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ] := sthdSBAcacheDresult[fkd, fkga, fh0, fd0, fga0]=
	  	Module[{sol, e, hs0, gas0, ds0, brac},
	  		Off[LessEqual::nord];
	   		If[fh0 == 0, hs0 = 1.*10^-15, hs0 = fh0];
	   		If[fd0 == 0, ds0 = 1.*10^-15, ds0 = fd0];
	   		If[fga0 == 0, gas0 = 1.*10^-15, gas0 = fga0];
	   		sol = sthdSBAcacheD[fkd, fkga, ds0, gas0] /. {h0 -> hs0};
	   		e = Select[Re[sol], 0 <= # <= ds0 &];
	   		brac = First[e];
	   		Return[brac];
	   	];


	sthdSBAcacheDData[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fstep_?NumericQ] :=sthdSBAcacheDData[fkd, fkga, fh0, fd0, fga0, fstep]=

	  	Module[{sol, i},

	   		sol = Table[{N[i], sthdSBAcacheDresult[fkd, fkga, i, fd0, fga0]}, {i, 0, fh0, fstep}];
	  		Return[sol];
	   ];

(* SBA Signal*)

	sthdSBAcacheSignal[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ] := sthdSBAcacheSignal[fkd, fkga, fh0, fd0, fga0, fsig0, fsigHD, fsigD]=
		Module[{sol, solsigHD,solsigD},
	   		If[fsigHD != 0, solsigHD = sthdSBAcacheHDresult[fkd, fkga, fh0, fd0, fga0]];
			If[fsigD != 0, solsigD = sthdSBAcacheDresult[fkd, fkga, fh0, fd0, fga0]];
			sol = fsig0 + fsigHD*solsigHD + fsigD*solsigD;
			Return[sol];
	   ];

	sthdSBAcacheSignallog[fpkd_?NumericQ, fpkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ] :=
		Module[{sol, solsigHD,solsigD},
	   		If[fsigHD != 0, solsigHD = sthdSBAcacheHDresult[(10^fpkd), (10^fpkga), fh0, fd0, fga0]];
			If[fsigD != 0, solsigD = sthdSBAcacheDresult[(10^fpkd), (10^fpkga), fh0, fd0, fga0]];
			sol = fsig0 + fsigHD*solsigHD + fsigD*solsigD;
			Return[sol];
	   	];

	sthdSBAcacheSignalData[fkd_?NumericQ, fkga_?NumericQ, fh0_?NumericQ, fd0_?NumericQ, fga0_?NumericQ, fsig0_?NumericQ, fsigHD_?NumericQ, fsigD_?NumericQ, fstep_?NumericQ] := sthdSBAcacheSignalData[fkd, fkga, fh0, fd0, fga0, fsig0, fsigHD, fsigD, fstep]=
	  	Module[{sol, i},
	   		sol = Table[{N[i], sthdSBAcacheSignal[fkd, fkga, i, fd0, fga0, fsig0, fsigHD, fsigD]}, {i, 0, fh0, fstep}];
	   		Return[sol];
	   	];


End[] (* End Private Context *)

EndPackage[]
