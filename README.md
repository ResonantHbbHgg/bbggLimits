# bbggLimits
Package for computing limits for the Run II analyses

## Instalation
First, setup the environment with the Higgs Combine tools: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHiggsAnalysisCombinedLimit#For_end_users_that_don_t_need_to   
Currently recommended CMSSW version: 71X
```
setenv SCRAM_ARCH slc6_amd64_gcc481
cmsrel CMSSW_7_1_5 ### must be a 7_1_X release  >= 7_1_5;  (7.0.X and 7.2.X are NOT supported either) 
cd CMSSW_7_1_5/src 
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v5.0.2   # try v5.0.1 if any issues occur
scramv1 b clean; scramv1 b # always make a clean build, as scram doesn't always see updates to src/LinkDef.h

```    
Get bbggLimits:   
```
cd ${CMSSW_BASE}/src/HiggsAnalysis/
git clone git@github.com:ResonantHbbHgg/bbggLimits.git
./Compile.sh
```

How to run it :
```
		
1) Edit your .json ( example in LimitSetting/jest.json )
if you want to change the value of "minMggMassFit" "maxMggMassFit" ... etc for one Mass in particular (Mx=300 for example ), just add the line :	
"param_300" :[minMggMassFit,maxMggMassFit,minMjjMassFit,maxMjjMassFit,minSigFitMgg,maxSigFitMgg,minSigFitMjj,maxSigFitMjj,minHigMggFit,maxHigMggFit,minHigMjjFit,maxHigMjj],
	in "signal"
if you want to run runCombine and BrazilianFlag at the same time than bbgg2DFit, just put runCombine and doBrazilianFlag accordingly.

2) Run bbgg2DFit MyJsonFile MyFolder
	with MyJsonFile is the Json file you have created and MyFolder is the name of the folder in witch bbgg2DFit  runCombine and BrazilianFlag will put all their output. The name of the directory will be MyFolder_v{version} with {version} the number provided in MyJsonFile
	If you don't provide MyFolder argument bbgg2DFit will create bbggToolsResults_v{version} by default.

3) If you want to run RunCombine and/or BrazilianFlag alone run :
	runCombine MyJsonFile MyFolder
	BrazilianFlag MyJsonFile MyFolder
```         
