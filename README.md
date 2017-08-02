# Package for computing limits for the Run II analyses

## Instalation
First, setup the environment with the Higgs Combination tools: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHiggsAnalysisCombinedLimit  
Currently working with 74X (check latest on HiggsCombine twiki).   


#### Step 1: Get Combine   
(Taken from combine twiki on June 27)   

```
export SCRAM_ARCH=slc6_amd64_gcc491
cmsrel CMSSW_7_4_7
cd CMSSW_7_4_7/src 
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v6.3.1
scramv1 b clean; scramv1 b # always make a clean build, as scram doesn't always see updates to src/LinkDef.h
```       

#### Step 2: Get HH support stuff    
(Needed for analytical reweighting. When trying to reproduce EPS results, check with Konstantin and Alexandra about which tag to use for the following repositories.)    

```
cd ${CMSSW_BASE}/src/
git clone git@github.com:cms-hh/HHStatAnalysis.git
scramv1 b HHStatAnalysis/AnalyticalModels #only build whats needed, no need Harvester (will keep complaining, though. If you want to compile the full HHStatAnalysis, also get CombineHarvester package.)
cd ${CMSSW_BASE}/src/HHStatAnalysis
git clone git@github.com:cms-hh/Support.git
```    

#### Step 3: Get bbggTools    
```
cd ${CMSSW_BASE}/src/HiggsAnalysis/
git clone -b dev-rafael-May8 git@github.com:ResonantHbbHgg/bbggLimits.git
cd ${CMSSW_BASE}/src/HiggsAnalysis/bbggLimits/
scramv1 b # a lot of complaints about bbggHighMassFitter.cc (this is not used anymore, needs to be deleted)
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


### Set the Limits
In order to reporduce EPS17 results, follow instructions here:
https://github.com/ResonantHbbHgg/bbggLimits/blob/master/SmartScripts/README.md

Good luck!
