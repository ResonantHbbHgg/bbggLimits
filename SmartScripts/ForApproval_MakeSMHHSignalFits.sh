#!/bin/bash

#PF=LT_MJJ70NewMVA_350_HMHPC_970_HMMPC_600_LMHPC_985_LMMPC_600
#KLSCAN_LIMS_LT_350_HMHPC_970_HMMPC_600_LMHPC_985_LMMPC_600_v66
#PF=LT_MJJ70NewMVA_350_HMHPC_970_HMMPC_600_LMHPC_985_LMMPC_600_allNodes
PF=LT_350_HMHPC_970_HMMPC_600_LMHPC_985_LMMPC_600
PF2=kl1p0_kt1p0_cg0p0_c20p0_c2g0p0
DIR=LIMS_${PF}_v66

DIR_LM=${DIR}/LowMass_Node_SM${PF2}
DIR_HM=${DIR}/HighMass_Node_SM${PF2}

OUT=~/www/HHBBGG/ForApproval/SignalShapes_${PF2}

mkdir ${OUT}

cp ~/www/HHBBGG/index.php ${OUT}

#python ../scripts/MakeSigPlot.py -w ${DIR_LM}/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 0 -o "mjj,mgg" -l 35.87 -a "Nonresonant Analysis|SM-like HH|Low mass region" -b 24,160 -L "SMHHLM" --DSCB
python scripts/MakeSigPlot.py -w ${DIR_LM}/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 0 -o "mjj,mgg" -l 35.87 -a "#font[61]{pp#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|SM-like HH|Low mass region" -b 24,160 -L "SMHHLM" --DSCB
mv SMHH* ${OUT}
#python ../scripts/MakeSigPlot.py -w ${DIR_LM}/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 1 -o "mjj,mgg" -l 35.87 -a "Nonresonant Analysis|SM-like HH|Low mass region" -b 24,160 -L "SMHHLM" --DSCB
python scripts/MakeSigPlot.py -w ${DIR_LM}/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 1 -o "mjj,mgg" -l 35.87 -a "#font[61]{pp#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|SM-like HH|Low mass region" -b 24,160 -L "SMHHLM" --DSCB
mv SMHH* ${OUT}
#python ../scripts/MakeSigPlot.py -w ${DIR_HM}/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 2 -o "mjj,mgg" -l 35.87 -a "Nonresonant Analysis|SM-like HH|High mass region" -b 24,160 -L "SMHHHM" --DSCB
python scripts/MakeSigPlot.py -w ${DIR_HM}/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 2 -o "mjj,mgg" -l 35.87 -a "#font[61]{pp#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|SM-like HH|High mass region" -b 24,160 -L "SMHHHM" --DSCB
mv SMHH* ${OUT}
#python ../scripts/MakeSigPlot.py -w ${DIR_HM}/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 3 -o "mjj,mgg" -l 35.87 -a "Nonresonant Analysis|SM-like HH|High mass region" -b 24,160 -L "SMHHHM" --DSCB
python scripts/MakeSigPlot.py -w ${DIR_HM}/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 3 -o "mjj,mgg" -l 35.87 -a "#font[61]{pp#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|SM-like HH|High mass region" -b 24,160 -L "SMHHHM" --DSCB
mv SMHH* ${OUT}
