#!/bin/bash

DIR=outDir_MILANO_v3_split480_Exp_Hist_2

DIR_LM=${DIR}/LowMass_Node_SM
DIR_HM=${DIR}/HighMass_Node_SM

OUT=${DIR}/CombinedCard_Node_SM/SignalShapes

mkdir ${OUT}


#python ../scripts/MakeSigPlot.py -w ${DIR_LM}/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 0 -o "mjj,mgg" -l 35.9 -a "Nonresonant Analysis|SM-like HH|Low mass region" -b 24,160 -L "SMHHLM" --DSCB
python scripts/MakeSigPlot.py -w ${DIR_LM}/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 0 -o "mjj,mgg" -l 35.9 -a "#font[61]{pp#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|SM-like HH|Low mass region" -b 24,160 -L "SMHHLM" --DSCB
mv SMHH* ${OUT}
#python ../scripts/MakeSigPlot.py -w ${DIR_LM}/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 1 -o "mjj,mgg" -l 35.9 -a "Nonresonant Analysis|SM-like HH|Low mass region" -b 24,160 -L "SMHHLM" --DSCB
python scripts/MakeSigPlot.py -w ${DIR_LM}/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 1 -o "mjj,mgg" -l 35.9 -a "#font[61]{pp#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|SM-like HH|Low mass region" -b 24,160 -L "SMHHLM" --DSCB
mv SMHH* ${OUT}
#python ../scripts/MakeSigPlot.py -w ${DIR_HM}/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 2 -o "mjj,mgg" -l 35.9 -a "Nonresonant Analysis|SM-like HH|High mass region" -b 24,160 -L "SMHHHM" --DSCB
python scripts/MakeSigPlot.py -w ${DIR_HM}/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 2 -o "mjj,mgg" -l 35.9 -a "#font[61]{pp#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|SM-like HH|High mass region" -b 24,160 -L "SMHHHM" --DSCB
mv SMHH* ${OUT}
#python ../scripts/MakeSigPlot.py -w ${DIR_HM}/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 3 -o "mjj,mgg" -l 35.9 -a "Nonresonant Analysis|SM-like HH|High mass region" -b 24,160 -L "SMHHHM" --DSCB
python scripts/MakeSigPlot.py -w ${DIR_HM}/workspaces/hhbbgg.mH125_13TeV.inputsig.root -c 3 -o "mjj,mgg" -l 35.9 -a "#font[61]{pp#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|SM-like HH|High mass region" -b 24,160 -L "SMHHHM" --DSCB
mv SMHH* ${OUT}
