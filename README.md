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
You can aslo provide the locations of the flat trees if they are not the ones hardcoded in
the script, via `-s`, `-d`, `-b` options. For example, to make the trees from single H, use: 
```
makeAllTrees.py -x nonres -f LT_OutDir -s FlatT_SignalDir -d 0 \  
--doCatMVA --MVAHMC0 0.970 --MVAHMC1 0.600 --MVALMC0 0.985 --MVALMC1 0.600 --massNR 350 --doSMHiggs --LMLJBTC 0.55 --LMSJBTC 0.55 --genDiPhotonFilte
```  

 
In order to re-produce the limit trees used for EPS17 results, follow instructions in
[SmartScripts/README.md](SmartScripts/README.md).


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

### Non-resonant reweighted trees 
In the non-resonant search, to get the limits at any parameter point of 5D space, we need
to reweight the signal sample to that point.  To do that we have a script,
`scripts/MakeARWTree.py`. Have a look at it to understand what it does. Then, to simplify
the production of those reweighted trees, we have another script which does everything on
batch, `scripts/ArwTreesOnLSF.py`. For example, to make trees for *kl* scan run:
```
python scripts/ArwTreesOnLSF.py -t KL
```


### Set the Limits 
Once the limit trees are produced, we would like to make the 2D fits in
(m(gg), m(bb)) plane, for each category, and then run the limits.

The main functions to do the fits are implemented in `src/bbgg2DFitter.cc`.  The python
scripts are needed to handle many different situations (resonant, non-resonant,
re-weighting to various benchmark points, etc.). In order to run just one limit you need
`scripts/pyLimits.py`. Minimal options for the *SM* point are:  
``` 
pyLimits.py -f conf_default.json -o outputDirName --nodes SM 
```  
*Important*: one has to specify the location of the inpot limit trees in
`conf_default.json` file.  The above command must be run on _lxplus_, if the input root
files are located on EOS.  
The `pyLimits.py` script would call _runFullChain()_ method which is implemented in
`python/LimitsUtil.py`.  So in fact, the [LimitsUtil.py](python/LimitsUtil.py) script is
the base code which interacts with the functions in `bbgg2DFitter.cc`.  

Using the `--nodes SM` option tells it to use the Limit Tree produced from a single SM MC
sample.  
Alternatively, one can do the limit os the re-weighted samples of the merged non-resonant
samples. (For the SM point this allows to increase the statistics of the signal
sample. Such re-weighting was used for EPS17 results of 2016 data.)  
Run it like so:  
```
pyLimits.py -f conf_default.json -o outputDirName --analyticalRW
```
In case of problems it's useful to increase verbosity level with `-v 1(2,3)` option. In
this case the logs should be found in your `/tmp/username/logs` and in the _master_ log,
`outputDirName/mainLog_date-time.log`


We have another script to facilitate runing the limit for _benchmarks_,_kl_ and _kl-kt_ scans:  
```
python scripts/runLimitsOnLSF.py -t [JHEP, KL, KLKT]
```

The above command should give you the limits identical to
[the ones on SVN](https://svnweb.cern.ch/cern/wsvn/cmshcg/trunk/cadi/HIG-17-008/NonResonant/Benchmarks/).

Good luck!

### Scripts for making plots
To make the non-resonant benchmark limit plot:  
```
python scripts/MakeBenchmarksPlot.py -f LIMSDIR
```
where _LIMSDIR_ is a directory with the limits output.

To get the *kl* scan plot:
```
python scripts/MakeKLambdaScan.py -f LIMSDIR
```

PS. In order to reporduce the _EPS17_ results, follow the instructions here:
[SmartScripts/README.md](SmartScripts/README.md)

