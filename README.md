# Package for computing limits for the Run II analyses

## Instalation
First, setup the environment with the Higgs Combine tools: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHiggsAnalysisCombinedLimit#For_end_users_that_don_t_need_to   
Currently working with 74X (check latest on HiggsCombine twiki).   

Get bbggLimits:   
```
cd ${CMSSW_BASE}/src/HiggsAnalysis/
git clone git@github.com:ResonantHbbHgg/bbggLimits.git
cd bbggLimits
scramv1 b -j 10
```

## Making Limit Trees

In order to make limit trees from all samples use these script:
```
makeAllTrees.py -x nonres [--NRW]
```

### Working examples

Make non-res shape benchmark points trees (MVA based with 400 Mhh threshold):
```
makeAllTrees.py -x nonres \   
-d /afs/cern.ch/work/r/rateixei/work/DiHiggs/flashgg_Moriond17/CMSSW_8_0_25/src/flashgg/bbggTools/test/RunJobs/Regression_Data/Hadd \   
-s /afs/cern.ch/work/r/rateixei/work/DiHiggs/flashgg_Moriond17/CMSSW_8_0_25/src/flashgg/bbggTools/test/RunJobs/Regression_Signal/Hadd/ \   
-f LT_NonRes_MVABased400Reg_ \   
--doCatMVA --MVAHMC0 0.960 --MVAHMC1 0.6 --MVALMC0 0.96 --MVALMC1 0.750 --massNR 400   
```   
   
Make non-res shape benchmark points trees (cut based with 400 Mhh threshold and cut on cos theta star):   
```
makeAllTrees.py -x nonres \   
-d /afs/cern.ch/work/r/rateixei/work/DiHiggs/flashgg_Moriond17/CMSSW_8_0_25/src/flashgg/bbggTools/test/RunJobs/newData_HHTagger/Hadd \   
-s /afs/cern.ch/work/r/rateixei/work/DiHiggs/flashgg_Moriond17/CMSSW_8_0_25/src/flashgg/bbggTools/test/RunJobs/Signal_HHTagger400/Hadd/ \   
-f LT_NonRes_CatBased400CTS_ --massNR 400 --ctsCut 0.8   
```   
   
Make resonant limit trees with low mass categorization:   
```
makeAllTrees.py -x res \   
-d /afs/cern.ch/work/r/rateixei/work/DiHiggs/flashgg_Moriond17/CMSSW_8_0_25/src/flashgg/bbggTools/test/RunJobs/newData_HHTagger/Hadd \   
-s /afs/cern.ch/work/r/rateixei/work/DiHiggs/flashgg_Moriond17/CMSSW_8_0_25/src/flashgg/bbggTools/test/RunJobs/Signal_HHTagger400/Hadd/ \   
-f LT_ResLMnW   
```   
   
Make resonant limit trees with high mass categorization:   
```
makeAllTrees.py -x res \   
-d /afs/cern.ch/work/r/rateixei/work/DiHiggs/flashgg_Moriond17/CMSSW_8_0_25/src/flashgg/bbggTools/test/RunJobs/newData_HHTagger/Hadd \   
-s /afs/cern.ch/work/r/rateixei/work/DiHiggs/flashgg_Moriond17/CMSSW_8_0_25/src/flashgg/bbggTools/test/RunJobs/Signal_HHTagger400/Hadd/ \   
-f LT_ResHM --highMassRes     
```    


### Details 
The *C++ Loop* code to produce the Limit Trees is located at
*src/bbggLTMaker.cc*. In order to run it over a single tree use the
python script *scripts/pyLimitTreeMaker.py*, which exists in the
*$PATH* after scram build. To run it just do:
```
pyLimitTreeMaker.py -f fileName.root -o outDir
```

where `fileName.root` is a an input Flat tree to be run over, and
`outDir` is where the output trees will be created. The
`makeAllTrees.py` mentioned in the beginning utilizes the
`pyLimitTreeMaker.py` and runs it over many input  files.


More options for the `pyLimitTreeMaker.py` can be specified:
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

Example to make all resonant limit trees:   
```
makeAllTrees.py -x res \   
-d /afs/cern.ch/work/r/rateixei/work/DiHiggs/flashgg_Moriond17/CMSSW_8_0_25/src/flashgg/bbggTools/test/RunJobs/DataJan20/Hadd \   
-s /afs/cern.ch/work/r/rateixei/work/DiHiggs/flashgg_Moriond17/CMSSW_8_0_25/src/flashgg/bbggTools/test/RunJobs/SignalJan20/Hadd \   
-f LT_Jan26_Res_
```   

### Set the Limits
For setting the Resonant limits, the instructions will come soon.

For Non-Resonant limit, stay here and run:
```
pyNonResLimits.py -f NonRes.json --nodes 2 3 SM --points 0-10
```
This example will run all fitting and limits for *Nodes 2,3,SM* and
the re-weighting to *points 0-10* out of 0-1506 avaialable. The results
of the combined limit should appear in *outDir/CombinedCard_Node_X/*
subdirectories.

* `-o <dir>` - Output directory (will be created)
* `--overwrite` - If the output directory exists - overwite it. By default, the script will exit
* `-v <integer>` - Verbosity level: 0 - Minimal or no messages; 1 - INFO; 2 - DEBUG; 3 - Go crazy.
If the verbosity level is greater than zero, log files are created at `/tmp/logs/` per process.
* `-j <ncpu>`  - Number of CPU cores to use. Default is 2.
* `-t <sec>` - Per job timeout (in seconds) for multiprocessing. Jobs will be killed if run longer than this.

Good luck!
