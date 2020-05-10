(* Wolfram Language Package *)

BeginPackage["datagen`"]
(* Exported symbols added here with SymbolName::usage *)
datagen::usage = "
	datagen[dataset-xy,y-minuent]
	";
normalize::usage = "
	normalize[value,max,min,[1-None,2-to N[0.1],3-Devide by Max]]
	";


Begin["`Private`"] (* Begin Private Context *)

datagen[dataset_, minuent_] := Module[{sol},
   sol = dataset /. {conc_, y_} -> {conc, minuent - y};
   Return[sol]
   ];

normalize[fvalue_, fmax_, fmin_, fnorm_] := Module[{a, b, newvalue},
   If[fnorm == 1,
    newvalue = fvalue;];
   If[fnorm == 2,
    a = 1/(fmax - fmin);
    b = 1 - a*fmax;
    newvalue = a*fvalue + b;];
   If[fnorm == 3,
    newvalue = fvalue/fmax;];
	   Return[newvalue]
   ];

End[] (* End Private Context *)

EndPackage[]
