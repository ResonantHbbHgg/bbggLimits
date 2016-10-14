### Making Limit Trees

The code to produce the Limit Trees is located at
*src/bbggLTMaker.cc*. In order to run it, we use the python script
under *scripts/pyLimitTreeMaker.py*, which exists in the *$PATH* after
scram build. To run it just do: 
```
pyLimitTreeMaker.py -i fList.txt -o outDir
```
where `fList.txt` is a list of root files to be run over, and
`outDir` is where the output trees will be created. The ```fList.txt``` can be obtained with the following simple command: `ls /path/to/flatTrees/output* > fList.txt`.


Other options can be specified:
* `-f <input File>` or `-i <Text file with a List of Root files full paths>`
* `-o  <output location>` - directory will be created.
* `--min  <min mtot>`, `--max <max mtot>`
* `--scale  <Lumi*CrossSection*SF/NEvts>` - scale; should be 1 for data.
* `--photonCR`  - do photon control region.
* `--KF`  - use Mtot_KF to cut on mass window.
* `--MX` -  use MX to cut on mass window; choose either `--MX` or `--KF`!.
* `--tilt`  - select tilted mass window.
* `--doNoCat`  - no categorization, all is *cat0*.
* `--btagWP <WP>` - set btagging working point for categories.
* `--doCatMixed` -  do categories with mixed btagging;  Cat0: 2>low, Cat1: 1<low+1>high
* `--singleCat`  - only one category, High Mass analysis.
* `--doBVariation <VAR>`  - apply b-tagging SF factors: 0, 1 or -1.
* `--doPhoVariation <VAR>`  - Apply photon SF factors: 0, 1 or -1.
* `--cosThetaStar <VAR>`  - cut on CosTheta Star variable

In order to make limit trees from all sample use these scripts:
```
makeAllTrees.py -x nonres
```

### Using C++ script to make Limit Trees (will be depricated soon):
The script is located in *bin/LimitTreeMaker.C*
```
LimitTreeMaker OPTIONS
```   

##### LimitTreeMaker Options:   
* -i <input list of files, text file with root files full paths> ( or -inputFile <single root file> )   
* -o <output location>   
* -min <min mtot> -max <max mtot>   
* -scale <Lumi*CrossSection*SF/NEvts, 1 for data>   
* -photonCR (do photon control region)   
* -KF (use Mtot_KF to cut on mass window)   
* -MX (use MX to cut on mass window) (choose either -MX or -KF!)   
* -tilt (select tilted mass window)   
* -doNoCat (no categorization, all is cat0)   
* -btagWP <WP> (set btagging working point for categories   
* -doCatMixed (do categories with mixed btagging - cat0: 2>low, cat1: 1<low+1>high)   
* -singleCat (only one category, High Mass analysis)   
* -doBVariation <VAR> (Apply b-tagging SF factors: 0, 1 or -1)
* -doPhoVariation <VAR> (Apply photon SF factors: 0, 1 or -1)
* -cosThetaStar <VAR> (Cut on CosTheta Star)
                                
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
bbggNonRes NonRes.json
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

