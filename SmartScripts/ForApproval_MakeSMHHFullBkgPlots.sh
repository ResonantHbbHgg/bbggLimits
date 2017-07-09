#!/bin/bash

## 1) Make MaxLikelihood fit

LIMFOLDER=LIMS_LT_350_HMHPC_970_HMMPC_600_LMHPC_985_LMMPC_600_v66
POINT=kl1p0_kt1p0_cg0p0_c20p0_c2g0p0 #this is SM

combine --datacard ${LIMFOLDER}/CombinedCard_Node_SM${POINT}/hhbbgg_13TeV_DataCard.txt -M MaxLikelihoodFit --saveWorkspace --saveShapes --saveNormalization --X-rtd TMCSO_AdaptivePseudoAsimov=50 -n SMHHForBkgPlots

mv mlfitSMHHForBkgPlots.root ${LIMFOLDER}/CombinedCard_Node_SM${POINT}/
mv higgsCombineSMHHForBkgPlots.MaxLikelihoodFit.mH120.root ${LIMFOLDER}/CombinedCard_Node_SM${POINT}/
mv MaxLikelihoodFitResult.root ${LIMFOLDER}/CombinedCard_Node_SM${POINT}/

## 2) Make plots with fit

INFILE=${LIMFOLDER}/CombinedCard_Node_SM${POINT}/MaxLikelihoodFitResult.root
OUTFILE_HM=FullBkgPlot_HM
OUTFILE_LM=FullBkgPlot_LM
#all four below to be taken from the datacard (signal rate)*33.49*0.0026
HMHP_NORM=0.37697
HMMP_NORM=0.35948
LMHP_NORM=0.02473
LMMP_NORM=0.04940
#these are multiplicative factors to make the signal show up
HMHP_FACT=20
HMMP_FACT=100
LMHP_FACT=100
LMMP_FACT=1000

HMTEXT="#font[61]{pp#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|High mass region"
LMTEXT="#font[61]{pp#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|Low mass region"

python scripts/MakeFullBackgroundFit.py -i ${INFILE} -o ${OUTFILE_HM} \
--signalNormalization ${HMHP_NORM} ${HMMP_NORM} \
--signalFactor ${HMHP_FACT} ${HMMP_FACT} --addHiggs \
--text "${HMTEXT}" --unblind

python scripts/MakeFullBackgroundFit.py -i ${INFILE} -o ${OUTFILE_LM} \
--signalNormalization ${LMHP_NORM} ${LMMP_NORM} \
--signalFactor ${LMHP_FACT} ${LMMP_FACT} --addHiggs \
--text "${LMTEXT}" --unblind

