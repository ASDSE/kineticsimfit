(* Wolfram Language Package *)

BeginPackage["styles`", { "styles`"}]
(* Exported symbols added here with SymbolName::usage *)


stylesButtonCreate::usage;

stylesButtonGenericStyle::usage;

headerStyling::usage;

exampleStyling::usage;

panelstyling::usage;

Begin["`Private`"] (* Begin Private Context *)

headerStyling[fname_] :=
 Module[{},
  Style[fname, Bold, Underlined, FontSize -> 16 ,
   TextAlignment -> Right, FontColor -> Gray]]

exampleStyling[fname_] :=
 Module[{},
  Style[fname, Bold, Italic, FontSize -> 16 , TextAlignment -> Right,
   FontColor -> LightGray]]

panelstyling[fname_] :=
 Module[{},
  Style[fname, Bold, FontSize -> 16 , TextAlignment -> Left]]



stylesButtonName[]:=stylesButtonName[]=Module[{buttonDefaultImage},

buttonDefaultImage=Import[FileNameJoin[{ParentDirectory[NotebookDirectory[]], "Images",
  "buttonDefault.png"}]];
  Return[buttonDefaultImage]
   ];

stylesButtonCreate[fName_, fTag_, fToolTip_, fEnabled_]:=stylesButtonCreate[fName, fTag, fToolTip, fEnabled]=Module[{buttonAppearance,button,nb,img},
	img=stylesButtonName[];
	buttonAppearance =
 	Framed[Style[fName, "Section", FontSize -> 18],
 		RoundingRadius -> 15, Background -> LightGray,
 		FrameStyle -> Directive[LightGray, 12]];

		If[fEnabled == 1,
    		button =
      Tooltip[MouseAppearance[
        Button[buttonAppearance,
         { nb = EvaluationNotebook[],
 				NotebookFind[nb, fTag, All, CellTags],
 			SelectionEvaluate[nb]}, Appearance -> None], "LinkHand"],
       fToolTip];,
    		button =
      Tooltip[MouseAppearance[
        Button[buttonAppearance,
         { nb = EvaluationNotebook[],
 		NotebookFind[nb, fTag, All, CellTags],
 		SelectionEvaluate[nb]}, Appearance -> None,
         Background -> LightGray, Enabled -> False], "LinkHand"],
       fToolTip];
    	];

   		Return[button]

];


stylesButtonGenericStyle[fName_, fToolTip_]:=stylesButtonGenericStyle[fName,fToolTip]=Module[{buttonAppearance},

	buttonAppearance = Tooltip[Framed[Style[fName, "Section", FontSize -> 18], RoundingRadius -> 15, Background -> LightGray, FrameStyle -> Directive[LightGray, 12]],fToolTip];
	Return[buttonAppearance]
	];

End[] (* End Private Context *)

EndPackage[]
