(* Wolfram Language Package *)

BeginPackage["buttons`", { "buttons`"}]
(* Exported symbols added here with SymbolName::usage *)


simulationbutton::usage;
importbutton::usage;
fitbutton::usage;
browsebutton::usage;

savepath::usage;
transferpath::usage;

Needs["styles`"];

Begin["`Private`"] (* Begin Private Context *)
simulationbutton :=
 MouseAppearance[
  Button[stylesButtonGenericStyle["Simulate",
    "Simulate the isotherm"], { nb = EvaluationNotebook[],
    NotebookFind[nb, "simulation", All, CellTags],
    SelectionEvaluate[nb]}, Appearance -> None], "LinkHand"];

importbutton :=
 MouseAppearance[
  Button[stylesButtonGenericStyle["Import/Update Data",
    "Import the dataset in ASCII format"], { nb = EvaluationNotebook[],
    NotebookFind[nb, "import", All, CellTags],
    SelectionEvaluate[nb]}, Appearance -> None], "LinkHand"];

fitbutton :=
 MouseAppearance[
  Button[stylesButtonGenericStyle["Fit the Data",
    "Fit a model to the data"], { nb = EvaluationNotebook[],
    NotebookFind[nb, "fitting", All, CellTags],
    SelectionEvaluate[nb]}, Appearance -> None], "LinkHand"];




browsebutton[fvar_] :=
    MouseAppearance[
    Tooltip[
      Framed[
        Style[
          FileNameSetter[
            Dynamic[fvar],
            "Open", {"Readable files" -> {"*.TXT"}, "All files" -> {"*"}}
          ],
          "Section", FontSize -> 18
        ],
      RoundingRadius -> 15, Background -> LightGray,
      FrameStyle -> Directive[LightGray, 12]], "Browse the file"],
    "LinkHand"];




tempsavebutton[fname_, fparameterslist_] :=
  Block[{path,button},
  path = FileNameJoin[{ParentDirectory[NotebookDirectory[]],"Parameter", fname <> ".mx"}];
  button = MouseAppearance[
            Button[stylesButtonGenericStyle["Save Parameters", "Saves the set values of the parameters"], {
              Export[path,fparameterslist]}, Appearance -> None], "LinkHand"];
  Return[button]
  ];

savepath[fname_]:=FileNameJoin[{ParentDirectory[NotebookDirectory[]],"Parameter", fname <> ".mx"}];
transferpath[fname_]:=FileNameJoin[{ParentDirectory[NotebookDirectory[]],"Parameter", "transfer" <> StringDrop[fname,3] <> ".mx"}];

temploadbutton[fname_, fparameterslist_] :=
  Block[{path,button,parameter},
  path = FileNameJoin[{ParentDirectory[NotebookDirectory[]],"Parameter", fname <> ".mx"}];
  button = MouseAppearance[
            Button[stylesButtonGenericStyle["Load Saved Parameters", "Loads the previously saved values"], {
              Import[path]}, Appearance -> None], "LinkHand"];
  Return[path]
  ];


End[] (* End Private Context *)

EndPackage[]
