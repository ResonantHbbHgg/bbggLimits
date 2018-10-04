# Package for computing limits for the Run II analyses

## Installation
First, setup the environment with the Higgs Combination tools: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHiggsAnalysisCombinedLimit  
Currently working with 74X (check latest on HiggsCombine twiki).   


#### Step 1: Get Combine   
(Taken from combine twiki on June 27)   

```
export SCRAM_ARCH=slc6_amd64_gcc491
cmsrel CMSSW_7_4_7
cd CMSSW_7_4_7/src 
cmsenv
mkdir ${CMSSW_BASE}/src/HiggsAnalysis/
cd ${CMSSW_BASE}/src/HiggsAnalysis/
git clone -b tth-bdt-tagger git@github.com:ResonantHbbHgg/bbggLimits.git
mv bbggLimits/CombinedLimit .
mv bbggLimits/HHStatAnalysis ../
scramv1 b clean
scramv1 b
```       
