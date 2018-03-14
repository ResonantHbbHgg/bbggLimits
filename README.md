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
scramv1 b
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
You can also provide the locations of the flat trees if they are not the ones hard-coded in
the script, via `-s`, `-d` options.  
In order to make the trees from the single Higgs samples, use `--doSMHiggs` option, and don't run over data (`-d 0`): 
```
makeAllTrees.py -x nonres -f LT_OutDir -d 0 --doSMHiggs --genDiPhotonFilter \  
--doCatMVA --MVAHMC0 0.970 --MVAHMC1 0.600 --MVALMC0 0.985 --MVALMC1 0.600 --massNR 350  --LMLJBTC 0.55 --LMSJBTC 0.55
```  

Before we can proceed furthed, some `hadd`ing needs to be done:
```
for m in LowMass HighMass; do hadd -f LT_OutDir_${m}/LT_output_bbHToGG_M-125_13TeV_amcatnlo.root LT_OutDir_${m}/LT_output_bbHToGG_M-125_4FS_yb*.root; done
for m in LowMass HighMass; do hadd -f LT_OutDir_${m}/LT_output_GluGluToHHTo2B2G_AllNodes.root LT_OutDir_${m}/LT_output_GluGluToHHTo2B2G_node_*.root; done
```  

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

More options for the `pyLimitTreeMaker.py` can be specified. To see all of them look
[directly in the code](https://github.com/ResonantHbbHgg/bbggLimits/blob/cc11d25a97392ee55116bac9d08b77f5f4128998/scripts/pyLimitTreeMaker.py#L11-L85).


### Non-resonant reweighted trees 
In the non-resonant search, to get the limits at any parameter point of 5D space, we need
to reweigh the signal sample to that point.  To do that we have a script,
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
pyLimits.py -f conf_default.json -o outDir --nodes SM 
```  
*Important*: one has to specify the location of the input limit trees in
`conf_default.json` file.  The above command must be run on _lxplus_, if the input root
files are located on EOS.  
The `pyLimits.py` script would call _runFullChain()_ method which is implemented in
`python/LimitsUtil.py`.  So in fact, the [LimitsUtil.py](python/LimitsUtil.py) script is
the base code which interacts with the functions in `bbgg2DFitter.cc`.  

Using the `--nodes SM` option tells it to use the Limit Tree produced from a single SM MC
sample.  
Alternatively, one can do the limit on the re-weighted samples of the merged non-resonant
samples. (For the SM point this allows to increase the statistics of the signal
sample. Such re-weighting was used for EPS17 results of 2016 data.)  
Run it like so:  
```
pyLimits.py -f conf_default.json -o outDir --analyticalRW
```
In case of problems it's useful to increase verbosity level with `-v 1(2,3)` option. In
this case the logs should be found in your `outDir/logs` and in the _master_ log,
`outDir/mainLog_date-time.log`

We have another script to facilitate running the limit for _benchmarks_, _kl_ and _kl-kt_ scans:  
```
python scripts/runLimitsOnLSF.py -f conf_default.json -t [JHEP, KL, KLKT] [-o outDir]
```

The above command should give you the limits identical to
[the ones on SVN](https://svnweb.cern.ch/cern/wsvn/cmshcg/trunk/cadi/HIG-17-008/NonResonant/Paper_v14/).

Good luck!

### Scripts for making plots
Make background fit plots for m(gg) and m(jj) in all categories:
```  
source scripts/MakeSMHHFullBkgPlots.sh LIMSDIR
```  
where _LIMSDIR_ is a directory with the limits output. Similarly, for the signal shape (SM
point), run:  
```
python scripts/MakeSigPlotSimple.py LIMSDIR
```  

In order make the non-resonant benchmark limit plot:  
```
python scripts/MakeBenchmarksPlot.py LIMSDIR
```

To get the *kl* scan plot:
```
python scripts/MakeKLambdaScan.py LIMSDIR
```

For *kl-kt* scan plot, we first need to gather the results of all limits in a text file,
and then run the plotting script:  
```
python scripts/MakeKLKTScanTxtList.py LIMSDIR [-s]
python scripts/MakeKLKTplot.py -l LIMSDIR/KLKT_Scan_List.txt
```  
Here, `-s` option can be used if only limits for `kt>0` are produced. In this case the
plot is simply drawn symmetrically over (0,0) point in (kl,kt) coordinates.


PS. In order to reproduce the _EPS17_ results, follow the instructions here:
[SmartScripts/README.md](SmartScripts/README.md)

