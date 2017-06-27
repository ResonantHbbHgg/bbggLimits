#!/bin/bash

#BASE=LT_MJJ70NewMVA_350_HMHPC_970_HMMPC_925_LMHPC_975_LMMPC_925
#BASE=LT_MJJ70NewMVA_350_HMHPC_970_HMMPC_600_LMHPC_985_LMMPC_600
BASE=LT_testPY_350_HMHPC_970_HMMPC_600_LMHPC_985_LMMPC_600
BASE2=${BASE}
DIR=../LIMS_${BASE}_v66
DIR_LM=${DIR}/LowMass_Node_SM${BASE2}
DIR_HM=${DIR}/HighMass_Node_SM${BASE2}
OUT=~/www/HHBBGG/SingleHiggsShapes_${BASE2}/

mkdir ${OUT}
cp ~/www/HHBBGG/index.php ${OUT}

python ../scripts/MakeHigPlot.py -w ${DIR_HM}/workspaces/hhbbgg.bbh.inputhig.root -c 2 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|High Mass Region|SM Higgs, b#bar{b}H" -b 12,160 -i bbh -L "bbh_HM" --DSCB
python ../scripts/MakeHigPlot.py -w ${DIR_HM}/workspaces/hhbbgg.bbh.inputhig.root -c 3 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|High Mass Region|SM Higgs, b#bar{b}H" -b 12,160 -i bbh -L "bbh_HM" --DSCB
python ../scripts/MakeHigPlot.py -w ${DIR_LM}/workspaces/hhbbgg.bbh.inputhig.root -c 1 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|Low Mass Region|SM Higgs, b#bar{b}H" -b 12,160 -i bbh -L "bbh_LM" --DSCB
python ../scripts/MakeHigPlot.py -w ${DIR_LM}/workspaces/hhbbgg.bbh.inputhig.root -c 0 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|Low Mass Region|SM Higgs, b#bar{b}H" -b 12,160 -i bbh -L "bbh_LM" --DSCB
mv *signal_fit* ${OUT}
python ../scripts/MakeHigPlot.py -w ${DIR_HM}/workspaces/hhbbgg.vbf.inputhig.root -c 2 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|High Mass Region|SM Higgs, VBF H" -b 12,160 -i vbf -L "vbf_HM" --DSCB
python ../scripts/MakeHigPlot.py -w ${DIR_HM}/workspaces/hhbbgg.vbf.inputhig.root -c 3 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|High Mass Region|SM Higgs, VBF H" -b 12,160 -i vbf -L "vbf_HM" --DSCB
python ../scripts/MakeHigPlot.py -w ${DIR_LM}/workspaces/hhbbgg.vbf.inputhig.root -c 1 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|Low Mass Region|SM Higgs, VBF H" -b 12,160 -i vbf -L "vbf_LM" --DSCB
python ../scripts/MakeHigPlot.py -w ${DIR_LM}/workspaces/hhbbgg.vbf.inputhig.root -c 0 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|Low Mass Region|SM Higgs, VBF H" -b 12,160 -i vbf -L "vbf_LM" --DSCB
mv *signal_fit* ${OUT}
python ../scripts/MakeHigPlot.py -w ${DIR_HM}/workspaces/hhbbgg.vh.inputhig.root -c 2 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|High Mass Region|SM Higgs, VH" -b 12,160 -i vh -L "vh_HM" --DSCB
python ../scripts/MakeHigPlot.py -w ${DIR_HM}/workspaces/hhbbgg.vh.inputhig.root -c 3 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|High Mass Region|SM Higgs, VH" -b 12,160 -i vh -L "vh_HM" --DSCB
python ../scripts/MakeHigPlot.py -w ${DIR_LM}/workspaces/hhbbgg.vh.inputhig.root -c 1 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|Low Mass Region|SM Higgs, VH" -b 12,160 -i vh -L "vh_LM" --DSCB
python ../scripts/MakeHigPlot.py -w ${DIR_LM}/workspaces/hhbbgg.vh.inputhig.root -c 0 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|Low Mass Region|SM Higgs, VH" -b 12,160 -i vh -L "vh_LM" --DSCB
mv *signal_fit* ${OUT}
python ../scripts/MakeHigPlot.py -w ${DIR_HM}/workspaces/hhbbgg.tth.inputhig.root -c 2 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|High Mass Region|SM Higgs, t#bar{t}H" -b 12,160 -i tth -L "tth_HM" --DSCB
python ../scripts/MakeHigPlot.py -w ${DIR_HM}/workspaces/hhbbgg.tth.inputhig.root -c 3 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|High Mass Region|SM Higgs, t#bar{t}H" -b 12,160 -i tth -L "tth_HM" --DSCB
python ../scripts/MakeHigPlot.py -w ${DIR_LM}/workspaces/hhbbgg.tth.inputhig.root -c 1 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|Low Mass Region|SM Higgs, t#bar{t}H" -b 12,160 -i tth -L "tth_LM" --DSCB
python ../scripts/MakeHigPlot.py -w ${DIR_LM}/workspaces/hhbbgg.tth.inputhig.root -c 0 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|Low Mass Region|SM Higgs, t#bar{t}H" -b 12,160 -i tth -L "tth_LM" --DSCB
mv *signal_fit* ${OUT}
python ../scripts/MakeHigPlot.py -w ${DIR_HM}/workspaces/hhbbgg.ggh.inputhig.root -c 2 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|High Mass Region|SM Higgs, ggH" -b 12,160 -i ggh -L "ggh_HM" --DSCB
python ../scripts/MakeHigPlot.py -w ${DIR_HM}/workspaces/hhbbgg.ggh.inputhig.root -c 3 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|High Mass Region|SM Higgs, ggH" -b 12,160 -i ggh -L "ggh_HM" --DSCB
python ../scripts/MakeHigPlot.py -w ${DIR_LM}/workspaces/hhbbgg.ggh.inputhig.root -c 1 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|Low Mass Region|SM Higgs, ggH" -b 12,160 -i ggh -L "ggh_LM" --DSCB
python ../scripts/MakeHigPlot.py -w ${DIR_LM}/workspaces/hhbbgg.ggh.inputhig.root -c 0 -o "mjj,mgg" -l 35.87 -a "Non-Resonant Analysis|Low Mass Region|SM Higgs, ggH" -b 12,160 -i ggh -L "ggh_LM" --DSCB
mv *signal_fit* ${OUT}

