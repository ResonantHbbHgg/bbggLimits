#!/bin/bash

PF=LT_MJJ70NewMVA_350_HMHPC_970_HMMPC_600_LMHPC_985_LMMPC_600

DIR=../LIMS_${PF}_v66

DIR_LM=${DIR}/LowMass_Node_SM${PF}
DIR_HM=${DIR}/HighMass_Node_SM${PF}

OUT=~/www/HHBBGG/BackgroundShapes_${PF}_UB

mkdir ${OUT}

cp ~/www/HHBBGG/index.php ${OUT}

python ../scripts/MakeBkgPlot.py -w ${DIR_LM}/workspaces/hhbbgg.inputbkg_13TeV.root -c 0 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|Low Mass Region|SM-like HH" -b 24,80 -L "SMHHLM" #-B
mv SMHH* ${OUT}
python ../scripts/MakeBkgPlot.py -w ${DIR_LM}/workspaces/hhbbgg.inputbkg_13TeV.root -c 1 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|Low Mass Region|SM-like HH" -b 24,80 -L "SMHHLM" #-B
mv SMHH* ${OUT}
python ../scripts/MakeBkgPlot.py -w ${DIR_HM}/workspaces/hhbbgg.inputbkg_13TeV.root -c 2 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|High Mass Region|SM-like HH" -b 24,80 -L "SMHHHM" #-B
mv SMHH* ${OUT}
python ../scripts/MakeBkgPlot.py -w ${DIR_HM}/workspaces/hhbbgg.inputbkg_13TeV.root -c 3 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|High Mass Region|SM-like HH" -b 24,80 -L "SMHHHM" #-B
mv SMHH* ${OUT}
