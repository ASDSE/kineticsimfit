(* Wolfram Language Package *)

BeginPackage["thermoCacheComp`"]
(* Exported symbols added here with SymbolName::usage *)
funTrimList::usage = "funTrimList[k]";

compElimFunHD::usage = "funTrimList[parm]";

compResultFunHD::usage = "funTrimList[input]";

Begin["`Private`"] (* Begin Private Context *)



compConstrFun[k_] :=
  Module[{eqBase, h0eq, i, parBase, eliBase, g0eq, kgeq},
   eqBase = {h0 == h + hd, d0 == d + hd, kd == hd/(h*d)};
   h0eq = {h0, h + hd};
   parBase = {h0, d0, kd};
   eliBase = {h, d};
   For[i = 0, i <= k, i++,
    If[i > 0,
      g0eq = {ToExpression["g0" <> ToString[i]] ==
         ToExpression["g" <> ToString[i]] +
          ToExpression["hg" <> ToString[i]]};
      kgeq = {ToExpression["kg" <> ToString[i]] ==
         ToExpression[
           "hg" <> ToString[i]]/(ToExpression["g" <> ToString[i]]*
            ToExpression["h"])};
      AppendTo[eqBase, g0eq[[1]]];
      AppendTo[parBase, Symbol["g0" <> ToString@i]];
      AppendTo[eliBase, Symbol["g" <> ToString@i]];
      AppendTo[eqBase, kgeq[[1]]];
      AppendTo[parBase, Symbol["kg" <> ToString@i]];
      AppendTo[eliBase, Symbol["hg" <> ToString@i]];
      h0eq = {h0eq[[1]],
        h0eq[[2]] + ToExpression["hg" <> ToString[i]]};
      eqBase = ReplacePart[eqBase, 1 -> h0eq[[1]] == h0eq[[2]]];

      ];

    ];

   Return[{eqBase, parBase, eliBase}]
   ];

compElimFunHD[parm_] :=
  compElimFunHD[parm] = Module[{result, sol, brac, l, eliHD},
    l = Length[parm];
    result = compConstrFun[(l - 3)/2];
    eliHD =
     Eliminate[result[[1]] /. Thread[result[[2]] -> parm],
      result[[3]]];
    sol = hd /. NSolve[eliHD, hd];
    Return[sol];

    ];

compResultFunHD[input_] :=
  compResultFunHD[input] =
   Module[{result, sol, brac, l, par, posZero, trimPar},
    par = N[input];
    Off[Select::normal];
    l = Length[par];
    If[(AllTrue[par, NumericQ] && OddQ[l]),
     If[(par[[1]] == 0. || par[[2]] == 0.),
      sol = 0.;,
      	If[MemberQ[Drop[par, 3], 0.],
        	trimPar = funTrimList[par];
        	sol = compElimFunHD[trimPar];,
        	sol = compElimFunHD[par];
        	];
      ];

     	result =
      Select[Re[sol], 0. <= # <= par[[1]] && 0. <= # <= par[[2]] &];
     	brac = First[result];
     , brac = "False"];

    Return[brac];
    ];
End[] (* End Private Context *)

EndPackage[]
