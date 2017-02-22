### Running Resonant Limits


#### Make Low mass+High mass plot:   
```
python scripts/LowHighResPlotter.py \
-L <Folder with low mass results> \
-H <Folder with high mass results> \
-l <Int. Lumi.> \
-G/-R (Graviton or Radion) \
-o (if plot observed)
```   

#### Make signal region fits with errors:   
```
python scripts/MakeBkgPlot.py \
-w DataWorkspace.root \
-c 0 [cat0 or cat1] \
-o mgg,mjj [observables] \
-l 2.70 [lumi] \
-a "pp#rightarrowX#rightarrowHH#rightarrowb#bar{b}#gamma#gamma|M_{X} = $i GeV Selection" [legend] \
-b 80,40 [bins for respective observables]   
```     

The script that reads the fitting results (expected number of signal events, which order of pol was used, etc) is in LimitSetting/scripts/DataCardMaker.py.   
The code is still creating the old version of datacards (Run1), this should be removed.


