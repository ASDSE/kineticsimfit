(* Wolfram Language Package *)

BeginPackage["inputlinsim`", { "inputlinsim`"}]
(* Exported symbols added here with SymbolName::usage *)


ilheading::usage;
ilinputfield::usage;
ilinfofield::usage;
placeholder::usage;
fitilheading::usage;
fitilinputfield::usage;
fitilinfofield::usage;


Begin["`Private`"] (* Begin Private Context *)

placeholder:= Block[{sol}, sol={"","",""}]

ilheading := Block[{heading},
   heading = {"Input parameters ", "Input Value", "", "", "Input Value", "", "", "Input Value"}; Return[heading]];

ilinputfield[ftype_, fname_, fvar_] :=
  Block[{after},
   after = Switch[ftype,
     "concentration",
      { Style["[" <> fname <> "\!\(\*SubscriptBox[\(]\), \(0\)]\) =", Bold, FontSize -> 16 , TextAlignment -> Right],
        InputField[Dynamic[fvar], FieldSize -> 8],
        "M"
      },
     "signal",
      {
       Style["Signal-" <> fname <> " =", Bold, FontSize -> 16 , TextAlignment -> Right],
       InputField[Dynamic[fvar], FieldSize -> 8],
       ""
     },
     "binding",
     {Style["\!\(\*SubscriptBox[\(K\), \(equ\)]\)(" <> fname <> ") =", Bold, FontSize -> 16 , TextAlignment -> Right],
      InputField[Dynamic[fvar], FieldSize -> 8],
      "\!\(\*SuperscriptBox[\(M\), \(-1\)]\)"
      },
     "kin",
     {Style[
        "\!\(\*SubscriptBox[\(k\), \(1\)]\)(" <> fname <> ") =", Bold, FontSize -> 16 , TextAlignment -> Right],
         InputField[Dynamic[fvar], FieldSize -> 8],
         "\!\(\*SuperscriptBox[\(M\), \(-1\)]\) \ \!\(\*SuperscriptBox[\(s\), \(-1\)]\)"
      },
     "kout",
      {Style["\!\(\*SubscriptBox[\(k\), \(-1\)]\)(" <> fname <> ") =", Bold, FontSize -> 16 , TextAlignment -> Right],
      InputField[Dynamic[fvar, None], FieldSize -> 8, Background -> LightGray],
      "\!\(\*SuperscriptBox[\(s\), \(-1\)]\)"
      },
     "time",
      {Style["\!\(\*SubscriptBox[\(t\), \(" <> fname <>"\)]\) =", Bold, FontSize -> 16 , TextAlignment -> Right],
      InputField[Dynamic[fvar], FieldSize -> 8],
      "sec"
      },
     "step",
      {Style["\!\(\*SubscriptBox[\("<>fname<>"\), \(steps\)]\) =", Bold, FontSize -> 16 , TextAlignment -> Right],
      InputField[Dynamic[fvar], FieldSize -> 8],
      "number"
      },
     "end",
      {Style["\!\(\*SubscriptBox[\("<>fname<>"\), \(endpoint\)]\) =", Bold, FontSize -> 16 , TextAlignment -> Right],
      InputField[Dynamic[fvar, None], FieldSize -> 8, Background -> LightGray],
      "number"
      }
     ];


     Return[after]];

ilinfofield[ftype_, fname_, fvar_] :=
  Block[{after}, after = Switch[ftype,
     "concentration", {Style["initial " <> fname <> " concentration",
       Thin, FontSize -> 14], Dynamic[fvar*1000000], "\[Mu]M"},
     "signal", {Style["Signal of " <> fname <> " complex", Thin,
       FontSize -> 14 ], Slider[Dynamic[fvar], {0, 10000000.}, ImageSize -> Small], ""},
     "binding", {Style["affinity of " <> fname <> " for host", Thin ,
       FontSize -> 14], Dynamic[Log10[fvar]], ""},
     "kin", {Style[fname <> " ingression rate", Bold, FontSize -> 14 ,
        TextAlignment -> Right],
      "\!\(\*SubscriptBox[\(k\), \(in\)]\)", ""},
     "kout", {Style[fname <> " egression rate", Bold, FontSize -> 14 ,
        TextAlignment -> Right],
      "\!\(\*SubscriptBox[\(k\), \(out\)]\)", ""},
     "time",{Style[fname <>" time", FontSize -> 14 ], Dynamic[fvar*1000], "msec"},
     "step", {Style["stepsize =", FontSize -> 14 ],Dynamic[fvar], fname},
     "end", {Style["endpoint =", FontSize -> 14 ], Dynamic[fvar], fname}
     ]; Return[after]];

fitilheading := Block[{heading},
  heading = {"Quantitiy", "Value", "Error", "Unit", "Fixation"}; Return[heading]];

fitilinputfield[ftype_, fname_, fvar_, fvarerror_, fvarfixation_, ffixation_] :=
 Block[{after,fixation},

  after = Switch[ftype,
    "concentration",
     {Style["[" <> fname <> "\!\(\*SubscriptBox[\(]\), \(0\)]\) =", Bold, FontSize -> 16 , TextAlignment -> Right],
       InputField[Dynamic[fvar], FieldSize -> 8],
       InputField[Dynamic[fvarerror, None], FieldSize -> 8, Background -> LightGray],
       "M",
        If[ffixation,RadioButtonBar[Dynamic[fvarfixation], {1 -> "Fixed", 2 -> "not fixed"}],""]
     },
    "signal",
     {
      Style["Signal-" <> fname <> " =", Bold, FontSize -> 16 , TextAlignment -> Right],
      InputField[Dynamic[fvar], FieldSize -> 8],
      InputField[Dynamic[fvarerror, None], FieldSize -> 8, Background -> LightGray],
      "a.u.",
       If[ffixation,RadioButtonBar[Dynamic[fvarfixation], {1 -> "Fixed", 2 -> "not fixed"}],""]
    },
   "autosignal",
    {
     Style["Signal-" <> fname <> " =", Bold, FontSize -> 16 , TextAlignment -> Right],
     InputField[Dynamic[fvar, None], FieldSize -> 8, FieldSize -> 8, Background -> LightGray],
     InputField[Dynamic[fvarerror, None], FieldSize -> 8, Background -> LightGray],
     "a.u.",
      If[ffixation,RadioButtonBar[Dynamic[fvarfixation], {1 -> "Fixed", 2 -> "not fixed"}],""]
   },
    "binding",
    {Style["\!\(\*SubscriptBox[\(K\), \(equ\)]\)(" <> fname <> ") =", Bold, FontSize -> 16 , TextAlignment -> Right],
     InputField[Dynamic[fvar], FieldSize -> 8],
     InputField[Dynamic[fvarerror, None], FieldSize -> 8, Background -> LightGray],
     "\!\(\*SuperscriptBox[\(M\), \(-1\)]\)",
      If[ffixation,RadioButtonBar[Dynamic[fvarfixation], {1 -> "Fixed", 2 -> "not fixed"}],""]

     },
    "kin",
    {Style[
       "\!\(\*SubscriptBox[\(k\), \(1\)]\)(" <> fname <> ") =", Bold, FontSize -> 16 , TextAlignment -> Right],
        InputField[Dynamic[fvar], FieldSize -> 8],
        InputField[Dynamic[fvarerror, None], FieldSize -> 8, Background -> LightGray],
        "\!\(\*SuperscriptBox[\(M\), \(-1\)]\) \ \!\(\*SuperscriptBox[\(s\), \(-1\)]\)",
         If[ffixation,RadioButtonBar[Dynamic[fvarfixation], {1 -> "Fixed", 2 -> "not fixed"}],""]

     },
    "kout",
     {Style["\!\(\*SubscriptBox[\(k\), \(-1\)]\)(" <> fname <> ") =", Bold, FontSize -> 16 , TextAlignment -> Right],
     InputField[Dynamic[fvar, None], FieldSize -> 8, Background -> LightGray],
     InputField[Dynamic[fvarerror, None], FieldSize -> 8, Background -> LightGray],
     "\!\(\*SuperscriptBox[\(s\), \(-1\)]\)",
      If[ffixation,RadioButtonBar[Dynamic[fvarfixation], {1 -> "Fixed", 2 -> "not fixed"}],""]

     }
    ]; Return[after]];


fitilinfofield[ftype_, fname_, fvar_, fvarerror_, fvarfixation_, ffixation_] :=
  Block[{after}, after = Switch[ftype,
     "concentration", {Style["initial " <> fname <> " concentration",
       Thin, FontSize -> 14], Dynamic[fvar*1000000], Dynamic[fvarerror*1000000], "\[Mu]M"},
     "signal", {Style["Signal of " <> fname <> " complex", Thin,
       FontSize -> 14 ], Slider[Dynamic[fvar], {0, 10000000.}, ImageSize -> Small], ""},
     "autosignal", {Style["Signal of " <> fname <> " complex", Thin,
       FontSize -> 14 ], "", ""},
     "binding", {Style["affinity of " <> fname <> " for host", Thin ,
       FontSize -> 14], Dynamic[Log10[fvar]], Dynamic[Log10[fvarerror]]},
     "kin", {Style[fname <> " ingression rate constant", Bold, FontSize -> 14 ,
        TextAlignment -> Right],
      "\!\(\*SubscriptBox[\(k\), \(in\)]\)", ""},
     "kout", {Style[fname <> " egression rate constant", Bold, FontSize -> 14 ,
        TextAlignment -> Right],
      "\!\(\*SubscriptBox[\(k\), \(out\)]\)", ""},
     "time",{Style[fname <>" time", FontSize -> 14 ], Dynamic[fvar*1000], "msec"},
     "step", {Style["stepsize =", FontSize -> 14 ],Dynamic[fvar], fname},
     "end", {Style["endpoint =", FontSize -> 14 ], Dynamic[fvar], fname}
     ]; Return[after]];
End[] (* End Private Context *)

EndPackage[]
