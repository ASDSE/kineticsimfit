(* Wolfram Language Package *)

(* Created by the Wolfram Workbench 27.04.2020 *)

BeginPackage["kineticSimFit`"]
(* Exported symbols added here with SymbolName::usage *)

createButton::usage = "
	Name of the Button, Path to File, ToolTip, Enabled=1";

Begin["`Private`"]
(* Implementation of the package *)

createButton[fName_, fPath_, fToolTip_, fEnabled_] :=
  	Module[{button, nb},
   	If[fEnabled == 1,
    		button =
      Tooltip[MouseAppearance[
        Button[Style[fName, "Section", FontSize -> 18],
         nb = NotebookOpen[fPath], Appearance -> None], "LinkHand"],
       fToolTip];,
    		button =
      Tooltip[MouseAppearance[
        Button[Style[fName, "Section", FontSize -> 18],
         nb = NotebookOpen[fPath], Appearance -> None,
         Background -> LightGray, Enabled -> False], "LinkHand"],
       fToolTip];
    	];

   		Return[button]

   		];

End[]

EndPackage[]
