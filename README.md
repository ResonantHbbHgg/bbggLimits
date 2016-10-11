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
#scramv1 b -j 10
./Compile.sh
```
   
### Make Limits Trees
See LimitSetting folder.   

### How to run it :
```
1) Edit your .json ( example in LimitSetting/jest.json ) :
	If you want to change the value of "minMggMassFit" "maxMggMassFit" ... etc for one Mass in particular just add the line in "signal":	
"param_Mass" :[minMggMassFit,maxMggMassFit,minMjjMassFit,maxMjjMassFit,minSigFitMgg,maxSigFitMgg,minSigFitMjj,maxSigFitMjj,minHigMggFit,maxHigMggFit,minHigMjjFit,maxHigMjj],
If you want to run runCombine and BrazilianFlag at the same time than bbgg2DFit, just put runCombine and doBrazilianFlag accordingly.

2) Run bbgg2DFit MyJsonFile MyFolder
	MyJsonFile is the Json file you have created 
	MyFolder is the name of the folder in witch bbgg2DFit runCombine and BrazilianFlag will put all their outputs. 
	The name of the directory will be MyFolder_v{version} with {version} the number provided in MyJsonFile
	If you don't provide MyFolder argument bbgg2DFit will create bbggToolsResults_v{version} by default.

3) If you want to run RunCombine and/or BrazilianFlag alone run :
	runCombine MyJsonFile MyFolder
	BrazilianFlag MyJsonFile MyFolder
```         
