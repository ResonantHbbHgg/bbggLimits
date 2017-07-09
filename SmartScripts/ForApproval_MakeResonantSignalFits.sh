#!/bin/bash

for i in 250 260 270 280 300 320 350 400 450 500 550 600 ; do
echo $i;
python scripts/MakeSigPlot.py -w LIMS_LT_Res_500_HMHPC_960_HMMPC_700_LMHPC_960_LMMPC_700_v66/BulkGraviton_Node_${i}LT_Res_500_HMHPC_960_HMMPC_700_LMHPC_960_LMMPC_700/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 0 -o "mjj,mgg" -l 35.87 -a "#font[61]{pp#rightarrowX#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|m_{X} = ${i} GeV, Spin-2" -b 24,80 -L "SigFit_Res${i}_Grav_" --DSCB
python scripts/MakeSigPlot.py -w LIMS_LT_Res_500_HMHPC_960_HMMPC_700_LMHPC_960_LMMPC_700_v66/BulkGraviton_Node_${i}LT_Res_500_HMHPC_960_HMMPC_700_LMHPC_960_LMMPC_700/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 1 -o "mjj,mgg" -l 35.87 -a "#font[61]{pp#rightarrowX#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|m_{X} = ${i} GeV, Spin-2" -b 24,80 -L "SigFit_Res${i}_Grav_" --DSCB
python scripts/MakeSigPlot.py -w LIMS_LT_Res_500_HMHPC_960_HMMPC_700_LMHPC_960_LMMPC_700_v66/Radion_Node_${i}LT_Res_500_HMHPC_960_HMMPC_700_LMHPC_960_LMMPC_700/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 0 -o "mjj,mgg" -l 35.87 -a "#font[61]{pp#rightarrowX#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|m_{X} = ${i} GeV, Spin-0" -b 24,80 -L "SigFit_Res${i}_Rad_" --DSCB
python scripts/MakeSigPlot.py -w LIMS_LT_Res_500_HMHPC_960_HMMPC_700_LMHPC_960_LMMPC_700_v66/Radion_Node_${i}LT_Res_500_HMHPC_960_HMMPC_700_LMHPC_960_LMMPC_700/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 1 -o "mjj,mgg" -l 35.87 -a "#font[61]{pp#rightarrowX#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|m_{X} = ${i} GeV, Spin-0" -b 24,80 -L "SigFit_Res${i}_Rad_" --DSCB
mv SigFit_Res* ~/www/HHBBGG/ForApproval/ResonantSigFits
done

for i in 600 650 700 750 800 900; do
echo $i;
python scripts/MakeSigPlot.py -w LIMS_LT_Res_500_HMHPC_500_HMMPC_000_LMHPC_500_LMMPC_000_v66/BulkGraviton_Node_${i}LT_Res_500_HMHPC_500_HMMPC_000_LMHPC_500_LMMPC_000/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 0 -o "mjj,mgg" -l 35.87 -a "#font[61]{pp#rightarrowX#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|m_{X} = ${i} GeV, Spin-2" -b 24,80 -L "SigFit_Res${i}_Grav_" --DSCB
python scripts/MakeSigPlot.py -w LIMS_LT_Res_500_HMHPC_500_HMMPC_000_LMHPC_500_LMMPC_000_v66/BulkGraviton_Node_${i}LT_Res_500_HMHPC_500_HMMPC_000_LMHPC_500_LMMPC_000/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 1 -o "mjj,mgg" -l 35.87 -a "#font[61]{pp#rightarrowX#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|m_{X} = ${i} GeV, Spin-2" -b 24,80 -L "SigFit_Res${i}_Grav_" --DSCB
python scripts/MakeSigPlot.py -w LIMS_LT_Res_500_HMHPC_500_HMMPC_000_LMHPC_500_LMMPC_000_v66/Radion_Node_${i}LT_Res_500_HMHPC_500_HMMPC_000_LMHPC_500_LMMPC_000/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 0 -o "mjj,mgg" -l 35.87 -a "#font[61]{pp#rightarrowX#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|m_{X} = ${i} GeV, Spin-0" -b 24,80 -L "SigFit_Res${i}_Rad_" --DSCB
python scripts/MakeSigPlot.py -w LIMS_LT_Res_500_HMHPC_500_HMMPC_000_LMHPC_500_LMMPC_000_v66/Radion_Node_${i}LT_Res_500_HMHPC_500_HMMPC_000_LMHPC_500_LMMPC_000/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 1 -o "mjj,mgg" -l 35.87 -a "#font[61]{pp#rightarrowX#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|m_{X} = ${i} GeV, Spin-0" -b 24,80 -L "SigFit_Res${i}_Rad_" --DSCB
mv SigFit_Res* ~/www/HHBBGG/ForApproval/ResonantSigFits
done

