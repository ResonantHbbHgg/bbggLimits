# bbggLimits
Package for computing limits for the Run II analyses

### Instalation
First, setup the environment with the Higgs Combine tools: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHiggsAnalysisCombinedLimit#For_end_users_that_don_t_need_to   
Currently working with 74X (check latest on HiggsCombine twiki).   

Get bbggLimits:   
```
cd ${CMSSW_BASE}/src/HiggsAnalysis/
git clone git@github.com:ResonantHbbHgg/bbggLimits.git
cd bbggLimits
scramv1 b -j 10
```

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
                                



### Set the Limits

For setting the Resonant limits, follow instructions in *LimitSetting* sub-directory.

For Non-Resonant limit, stay here and run:
```
pyNonResLimits.py -f NonRes.json

```
The results of the combined limit should appear in *outDir/CombinedCard_Node_X/* subdirectories.

Good luck!
