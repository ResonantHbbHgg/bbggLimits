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
scramv1 b clean
scramv1 b
```       

#### Step 2: Get HH support stuff    
(Needed for analytical reweighting. When trying to reproduce EPS results, check with Konstantin and Alexandra about which tag to use for the following repositories.)    

```
cd ${CMSSW_BASE}/src/
git clone git@github.com:cms-hh/HHStatAnalysis.git
scramv1 b HHStatAnalysis/AnalyticalModels
cd ${CMSSW_BASE}/src/HHStatAnalysis
git clone git@github.com:cms-hh/Support.git
```    

#### Step 3: Get bbggLimits code    
```
cd ${CMSSW_BASE}/src/HiggsAnalysis/
git clone  git@github.com:ResonantHbbHgg/bbggLimits.git
cd ${CMSSW_BASE}/src/HiggsAnalysis/bbggLimits/
scramv1 b # a lot of complaints about bbggHighMassFitter.cc (this is not used anymore, needs to be deleted)
```


## Making Limit Trees

In order to make limit trees from all samples use these script:
```
makeAllTrees.py -x nonres [--NRW]
```

### Working examples

Make non-res shape benchmark points trees (MVA based with 350 M(HH) threshold):
```
makeAllTrees.py -x nonres -f LT_OutDir \   
--doCatMVA --MVAHMC0 0.970 --MVAHMC1 0.600 --MVALMC0 0.985 --MVALMC1 0.600 --massNR 350 --LMLJBTC 0.55 --LMSJBTC 0.55
```   
You can aslo provide the locations of the flat tree ntuples via `-s`, `-d`, `-b` options.

Make non-res shape benchmark points trees (cut based with 400 Mhh threshold and cut on cos theta star):   
```
TBD
```   
   
Make resonant limit trees with low mass categorization:   
```
TBD
```   
   
Make resonant limit trees with high mass categorization:   
```
TBD
```    

In order to re-produce the limit trees used for EPS17 results, follow instructions in [SmartScripts/README.md](SmartScripts/README.md)


#### Details 
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


More options for the `pyLimitTreeMaker.py` can be specified, 
as can be seen [directly in the code](https://github.com/ResonantHbbHgg/bbggLimits/blob/10c319b013134e5bb15a561557f960dc2f1ea6b2/scripts/pyLimitTreeMaker.py#L11-L85)

### Set the Limits 
Once the limit trees are produced, we would like to make the 2D fits in
(m(gg), m(bb)) plane, for each category, and then run the limits.

The main functions to do the fits are implemented in `src/bbgg2DFitter.cc`.  The python
scripts are needed to handle many different situations (resonant, non-resonant,
re-weighting to various benchmark points, etc.). In order to run just one limit you need
`scripts/pyLimits.py`. Minimal options for the *SM* point are:  
``` 
./pyLimits.py -fconf_NonRes_EPS17.json -o outputDirName --nodes SM 
```

The above command must be run on _lxplus_, because the input root files are located on EOS
(the path is specified in json config file).  
The `pyLimits.py` script would call _runFullChain()_ method which is implemented in
`python/LimitsUtil.py`.  So in fact, the [LimitsUtil.py](python/LimitsUtil.py) script is
the base code which interacts with the functions in `bbgg2DFitter.cc`.  
Using the `--nodes SM` option tells it to use the Limit Tree produced from a single SM MC
sample.  Alternatively, one can do the re-weighting of all existing non-resonant
samples and therefore increase the statistics of the SM signal (number of events in a single
sample is only 50K). Analytical re-weighting was used for EPS17 results of 2016 data.  
Run it like so: 
```
./pyLimits.py -f conf_NonRes_EPS17.json -o outputDirName --analyticalRW 
```

The above command should give you the limits identical to
[the ones on SVN](https://svnweb.cern.ch/cern/wsvn/cmshcg/trunk/cadi/HIG-17-008/NonResonant/Benchmarks/CombinedCard_Node_SMkl1p0_kt1p0_cg0p0_c20p0_c2g0p0/result_2_L_CombinedCard_Node_SMkl1p0_kt1p0_cg0p0_c20p0_c2g0p0.log).
In order to reporduce the rest of _EPS17_ results, follow the instructions here:
[SmartScripts/README.md](SmartScripts/README.md)

Good luck!
