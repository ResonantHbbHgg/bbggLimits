#!/bin/bash

for i in 250 260 270 280 300 320 350 400 450 500 550 600 ; do
echo $i;
python scripts/MakeBkgPlot.py -w LIMS_LT_Res_500_HMHPC_960_HMMPC_700_LMHPC_960_LMMPC_700_v66/BulkGraviton_Node_${i}LT_Res_500_HMHPC_960_HMMPC_700_LMHPC_960_LMMPC_700/workspaces/hhbbgg.inputbkg_13TeV.root -c 0 -o "mjj,mgg" -l 35.87 -a "#font[61]{pp#rightarrowX#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|m_{X} = ${i} GeV Selection" -b 24,80 -L "Res${i}" #-B
python scripts/MakeBkgPlot.py -w LIMS_LT_Res_500_HMHPC_960_HMMPC_700_LMHPC_960_LMMPC_700_v66/BulkGraviton_Node_${i}LT_Res_500_HMHPC_960_HMMPC_700_LMHPC_960_LMMPC_700/workspaces/hhbbgg.inputbkg_13TeV.root -c 1 -o "mjj,mgg" -l 35.87 -a "#font[61]{pp#rightarrowX#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|m_{X} = ${i} GeV Selection" -b 24,80 -L "Res${i}" #-B
mv Res* ~/www/HHBBGG/ForApproval/ResonantBkgFits
done

for i in 600 650 700 750 800 900; do
echo $i;
python scripts/MakeBkgPlot.py -w LIMS_LT_Res_500_HMHPC_500_HMMPC_000_LMHPC_500_LMMPC_000_v66/BulkGraviton_Node_${i}LT_Res_500_HMHPC_500_HMMPC_000_LMHPC_500_LMMPC_000/workspaces/hhbbgg.inputbkg_13TeV.root -c 0 -o "mjj,mgg" -l 35.87 -a "#font[61]{pp#rightarrowX#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|m_{X} = ${i} GeV Selection" -b 24,80 -L "Res${i}" #-B
python scripts/MakeBkgPlot.py -w LIMS_LT_Res_500_HMHPC_500_HMMPC_000_LMHPC_500_LMMPC_000_v66/BulkGraviton_Node_${i}LT_Res_500_HMHPC_500_HMMPC_000_LMHPC_500_LMMPC_000/workspaces/hhbbgg.inputbkg_13TeV.root -c 1 -o "mjj,mgg" -l 35.87 -a "#font[61]{pp#rightarrowX#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|m_{X} = ${i} GeV Selection" -b 24,80 -L "Res${i}" #-B
mv Res* ~/www/HHBBGG/ForApproval/ResonantBkgFits
done


#python scripts/MakeBkgPlot.py -w LIMS_LT_Res_500_HMHPC_500_HMMPC_000_LMHPC_500_LMMPC_000_v66/BulkGraviton_Node_500LT_Res_500_HMHPC_500_HMMPC_000_LMHPC_500_LMMPC_000/workspaces/hhbbgg.inputbkg_13TeV.root -c 0 -o "mjj,mgg" -l 35.87 -a "#font[61]{pp#rightarrowX#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|m_{X} = 500 GeV Selection" -b 24,80 -L "Grav500" #-B
