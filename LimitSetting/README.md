### Running Resonant Limits
```
bbgg2DFit jest.json
```

Prepare your limit trees on different folders for each mass selection. Name each folder with the mass selection mass hypothesis at the end:   
```
ResonantAnalysis_250   
ResonantAnalysis_260   
ResonantAnalysis_320
...
```

##### jest.json:   
* Edit "dir" in data and signal with the path to your directories containing limit trees, and replace the hypothesis mass by "MASS"; the code will replace "MASS" by the masses defined in the "mass" section of the json file.   
* Replace the location of your signalModelCard (LimitSetting/Models/models_2D_sig.rs)   


* If you want to change the value of "minMggMassFit" "maxMggMassFit" ... etc for one Mass in particular just add the line in "signal":	
"param_Mass" :[minMggMassFit,maxMggMassFit,minMjjMassFit,maxMjjMassFit,minSigFitMgg,maxSigFitMgg,minSigFitMjj,maxSigFitMjj,minHigMggFit,maxHigMggFit,minHigMjjFit,maxHigMjj],
* If you want to run runCombine and BrazilianFlag at the same time than bbgg2DFit, just put runCombine and doBrazilianFlag accordingly.

* MyFolder is the name of the folder in witch bbgg2DFit runCombine and BrazilianFlag will put all their outputs. 
The name of the directory will be MyFolder_v{version} with {version} the number provided in MyJsonFile
If you don't provide MyFolder argument bbgg2DFit will create bbggToolsResults_v{version} by default.

* If you want to run RunCombine and/or BrazilianFlag alone run :
	runCombine MyJsonFile MyFolder
	BrazilianFlag MyJsonFile MyFolder
    


#### Make Low mass+High mass plot:   
```
python scripts/LowHighResPlotter.py \
-L <Folder with low mass results> \
-H <Folder with high mass results> \
-l <Int. Lumi.> \
-G/-R (Graviton or Radion) \
-o (if plot observed)
```   

#### Make signal region fits with errors:   
```
python scripts/MakeBkgPlot.py \
-w DataWorkspace.root \
-c 0 [cat0 or cat1] \
-o mgg,mjj [observables] \
-l 2.70 [lumi] \
-a "pp#rightarrowX#rightarrowHH#rightarrowb#bar{b}#gamma#gamma|M_{X} = $i GeV Selection" [legend] \
-b 80,40 [bins for respective observables]   
```     

### Running Non-Resonant Limits
```
bbggNonResFit NonRes.json
```   

Follow same procedure of Resonant limits, but with NonRes.json.   
Instead of defining masses, define the NonRes nodes you want to run:
* 0: Box only
* 1: SM
* 2, 3, 4, 5...

Instead of separating the directories in terms of masses, separate in terms of HighMass or LowMass selection:
```
NonResonantAnalysis_HighMass
NonResonantAnalysis_LowMass
```
Edit "dir" in the json file to contain the full path to these directories, replacing HighMass/LowMass by "TYPE".

## Datacards for Combine
The datacards to be used are created based on template datacards stored in LimitSetting/Models directory:
```
HighMassResDatacardModel.txt
LowMassResDatacardModel.txt
NonResDatacardModel.txt
```
The script that reads the fitting results (expected number of signal events, which order of pol was used, etc) is in LimitSetting/scripts/DataCardMaker.py.   
The code is still creating the old version of datacards (Run1), this should be removed.


