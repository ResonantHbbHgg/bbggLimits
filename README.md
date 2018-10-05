# Package for computing limits for the Run II analyses



### Step 1: Installation


```
export SCRAM_ARCH=slc6_amd64_gcc491
cmsrel CMSSW_7_4_7
cd CMSSW_7_4_7/src 
cmsenv
mkdir ${CMSSW_BASE}/src/HiggsAnalysis/
cd ${CMSSW_BASE}/src/HiggsAnalysis/
git clone -b ProjectionECFA git@github.com:ResonantHbbHgg/bbggLimits.git
mv bbggLimits/CombinedLimit .
mv bbggLimits/HHStatAnalysis ../
scramv1 b clean
scramv1 b
```       

### Step2: To create a datacard and analyse

```
pyLimits.py -f json/conf_MILANO.json -o outDir_MILANO_v3_split480 --nodes SM
```

This job may take forever. Need to fix that. In practice you simply need to do ctrl+z once you see that this file exist

outDir_MILANO_v3_split480/CombinedCard_Node_SM/hhbbgg_13TeV_DataCard.txt

Shall take 2-3 mn.

Then to produce Significance and errors:

sh SmartScripts/Analyzer.sh outDir_MILANO_v3_split480

You can find in outDir_MILANO_v3_split480/CombinedCard_Node_SM/ 2 files that contains significance, uncertainty and number of events per category

ProfileLikelihood.txt
MaxLikelihood.txt

### Step3

To produce background plots you need to add by hand the signal normalisation to  SmartScripts/ForApproval_MakeSMHHFullBkgPlots_UPGRADE.sh. Once I have time I would write an automatic script

HMHP_NORM=13
HMMP_NORM=9.6
LMHP_NORM=10.2
LMMP_NORM=14

you can also tune the y-axis in MakeFullBackgroundFit.py (need to automatise)

yLimits = {'mgg': [60, 700, 60, 400], 'mjj': [110, 1000, 110, 700]}


```
sh SmartScripts/ForApproval_MakeSMHHFullBkgPlots_UPGRADE.sh outDir_MILANO_v3_split480
```

the plots are in  outDir_MILANO_v3_split480/CombinedCard_Node_SM/Background

To produce signal plots

```
sh SmartScripts/ForApproval_MakeSMHHSignalFits_UPGRADE.sh outDir_MILANO_v3_split480
```

the plots are in  outDir_MILANO_v3_split480/CombinedCard_Node_SM/SignalShapes

To produce initial fit plot for background you need a bit of gymnastics

```
cp -r outDir_MILANO_v3_split480/CombinedCard_Node_SM outDir_MILANO_v3_split480/CombinedCard_Node_SM_Prefit
```

go to  outDir_MILANO_v3_split480/CombinedCard_Node_SM_Prefit/hhbbgg_13TeV_DataCard.txt  and remove everything from 

yield_norm    rateParam ch1_cat2 vbf 3000  

till the end. Then run

```
cd ${CMSSW_BASE}/src/HiggsAnalysis/bbggLimits
combine -M MaxLikelihoodFit -d outDir_MILANO_v3_split480/CombinedCard_Node_SM_Prefit/hhbbgg_13TeV_DataCard.txt --saveWorkspace --saveShapes --saveNormalization -S 0
mv MaxLikelihoodFitResult.root outDir_MILANO_v3_split480/CombinedCard_Node_SM_Prefit/
sh SmartScripts/ForApproval_MakeSMHHFullBkgPlots_UPGRADE_DATA.sh outDir_MILANO_v3_split480_Exp_Hist_5
```

the plots are in CombinedCard_Node_SM_Prefit/Background

you can also tune the proper scale of y-axis in MakeFullBackgroundFit.py

yLimits = {'mgg': [0.03, 0.1, 0.03, 0.1], 'mjj': [0.03, 0.1, 0.03, 0.1]}

======================

Enjoy...

