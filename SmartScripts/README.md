# Steps for full analysis

In order to reproduce exactly the results put out for EPS-2017, checkout the *v2.1-goodbye-rafael* tag of this repository.
Set up the environment from the main [README.md](README.md), then follow instructions below.

## Nonresonant Analysis

#### 0) Repeat HIG-17-008

svn co -N svn+ssh://mgouzevi@svn.cern.ch/reps/cmshcg/trunk/ SVNREP

cd SVNREP/

svn update -N cadi

svn update cadi/HIG-17-008

cd cadi/HIG-17-008/

combine -M Asymptotic LIMS_Paper_v14/CombinedCard_ARW_kl_1p0_kt_1p0_cg_0p0_c2_0p0_c2g_0p0/hhbbgg_13TeV_DataCard.txt

#### 1) Make limit trees   

```
cp SmartScripts/ForApproval_MakeNonResonantLimitTrees.py .
python ForApproval_MakeNonResonantLimitTrees.py
```

**Note**: The script above can also be used to optimize the MVA and Mhh categorization, by adding other working points to the for loops.

#### 2) Make JHEP benchmark limit trees   
NOTE: Run more than once to make sure all limit trees have been produced (needs to investigate why some fail).   

```
cp SmartScripts/ForApproval_MakeJHEPLimitTrees.py .
python ForApproval_MakeJHEPLimitTrees.py
```

#### 3) Make limit trees for kappa lambda and kl x kt scans   

NOTE: Run more than once to make sure all limit trees have been produced (needs to investigate why some fail).   

```
cp SmartScripts/ForApproval_MakeARWBSMLimitTrees.py .
python ForApproval_MakeARWBSMLimitTrees.py
```

#### 4) Run JHEP benchmark limits
```
cp SmartScripts/ForApproval_MakeJHEPLimits.py .
python ForApproval_MakeJHEPLimits.py
```
##### 4.1) Check that all the limits were computed   

NOTE: Some limits fail, not sure why. Need to run again (needs to investigate why some fail, is it a problem with the script?).

```
cp SmartScripts/ForApproval_MakeJHEPTxtList.py .
python ForApproval_MakeJHEPTxtList.py
```

If there are points in the "missing" text file, run:   

```
cp SmartScripts/ForApproval_MakeJHEPLimits_Missing.py .
python ForApproval_MakeJHEPLimits_Missing.py
```

#### 5) Make JHEP benchmarks plot

```
cp SmartScripts/ForApproval_MakeJHEPplot.sh .
source ForApproval_MakeJHEPplot.sh
```

#### 6) Run kappa lambda scan limits:

```
cp SmartScripts/ForApproval_MakeKLScanLimits.py .
python ForApproval_MakeKLScanLimits.py
```

##### 6.1) Check that all the limits were computed

NOTE: Some limits fail, not sure why. Need to run again (needs to investigate why some fail, is it a problem with the script?).

```
cp SmartScripts/ForApproval_MakeKLScanTxtList.py .
python ForApproval_MakeKLScanTxtList.py
```

If there are points in the "missing" text file, run:

```
cp SmartScripts/ForApproval_MakeKLScanLimits_Missing.py .
python ForApproval_MakeKLScanLimits_Missing.py
```

#### 7) Make kappa lambda scan plot:

```
cp SmartScripts/ForApproval_MakeKLplot.sh .
source ForApproval_MakeKLplot.sh
```

#### 8) Make kl x kt scan limits:

```
cp SmartScripts/ForApproval_MakeKLKTScanLimits.py .
python ForApproval_MakeKLKTScanLimits.py
```

##### 8.1) Check that all the limits were computed

NOTE: Some limits fail, not sure why. Need to run again (needs to investigate why some fail, is it a problem with the script?).

```
cp SmartScripts/ForApproval_MakeKLKTScanTxtList.py .
python ForApproval_MakeKLKTScanTxtList.py
``` 

If there are points in the "missing" text file, run:

```
cp SmartScripts/ForApproval_MakeKLKTScanLimits_Missing.py .
python ForApproval_MakeKLKTScanLimits_Missing.py
```

#### 9) Make kl x kt plot:

You need to run 6.1 first, even if all the limits are done. It will create a txt file with all points and limits that will be used by the plotting macro.

```
cp SmartScripts/ForApproval_MakeKLKTplot.sh .
source ForApproval_MakeKLKTplot.sh
```

## Resonant Analysis

#### 1) Make resonant limit trees and run limits

```
cp SmartScripts/ForApproval_MakeResonantAnalysis.py .
python ForApproval_MakeResonantAnalysis.py
```

This script will produce the limit trees and run the limits for both the low mass and high mass categorizations, for radion and graviton signal hypotheses. 

#### 2) Make resonant limit plot

```
cp SmartScripts/ForApproval_MakeResPlot.sh
source ForApproval_MakeResPlot.sh
```

## Support plots

#### 1) Make nonresonant background plots (background+signal+single Higgs):

```
cp SmartScripts/ForApproval_MakeSMHHFullBkgPlots.sh .
source ForApproval_MakeSMHHFullBkgPlots.sh
```

#### 2) Make nonresonant signal-like single Higgs fits:

```
cp SmartScripts/ForApproval_MakeSMHHHiggsSignalFits.sh .
source ForApproval_MakeSMHHHiggsSignalFits.sh
```

#### 3) Make SMHH signal fits:

```
cp SmartScripts/ForApproval_MakeSMHHSignalFits.sh .
source ForApproval_MakeSMHHSignalFits.sh
```

#### 3) Make resonant signal fits:

```
cp SmartScripts/ForApproval_MakeResonantSignalFits.sh .
source ForApproval_MakeResonantSignalFits.sh
```

#### 4) Make resonant background fits (only background and fit error):

```
cp SmartScripts/ForApproval_MakeResonantBackgroundFits.sh .
source ForApproval_MakeResonantBackgroundFits.sh
```

#### 4) Make SM HH limit category breakdown

```
cp SmartScripts/ForApproval_MakeSMHHCatsPlot.sh .
source ForApproval_MakeSMHHCatsPlot.sh
```
