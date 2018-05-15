#!/bin/bash

## 1) Make MaxLikelihood fit

LIMFOLDER=$1
POINT="kl_1p0_kt_1p0_cg_0p0_c2_0p0_c2g_0p0" #this is SM
CombDir="CombinedCard_ARW_${POINT}"

#CombDir="CombinedCard_Node_SM"

combine --datacard ${LIMFOLDER}/${CombDir}/hhbbgg_13TeV_DataCard.txt -M MaxLikelihoodFit --saveWorkspace --saveShapes --saveNormalization --X-rtd TMCSO_AdaptivePseudoAsimov=50 -n SMHHForBkgPlots

mv mlfitSMHHForBkgPlots.root ${LIMFOLDER}/${CombDir}/
mv higgsCombineSMHHForBkgPlots.MaxLikelihoodFit.mH120.root ${LIMFOLDER}/${CombDir}/
mv MaxLikelihoodFitResult.root ${LIMFOLDER}/${CombDir}/

## 2) Make plots with fit

INFILE=${LIMFOLDER}/${CombDir}/MaxLikelihoodFitResult.root
OUTFILE_HM=${LIMFOLDER}/FullBkgPlot_HM
OUTFILE_LM=${LIMFOLDER}/FullBkgPlot_LM
#all four below to be taken from the datacard (signal rate)*33.49*0.0026
HMHP_NORM=0.37144
HMMP_NORM=0.35764
LMHP_NORM=0.02331
LMMP_NORM=0.04631
#these are multiplicative factors to make the signal show up
HMHP_FACT=20
HMMP_FACT=100
LMHP_FACT=100
LMMP_FACT=1000

HMTEXT="#font[61]{pp#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|High-mass region"
LMTEXT="#font[61]{pp#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|Low-mass region"

python scripts/MakeFullBackgroundFit.py -i ${INFILE} -o ${OUTFILE_HM} \
       --signalNormalization ${HMHP_NORM} ${HMMP_NORM} \
       --signalFactor ${HMHP_FACT} ${HMMP_FACT} \
       --addHiggs \
       --text "${HMTEXT}" --unblind

python scripts/MakeFullBackgroundFit.py -i ${INFILE} -o ${OUTFILE_LM} \
       --signalNormalization ${LMHP_NORM} ${LMMP_NORM} \
       --signalFactor ${LMHP_FACT} ${LMMP_FACT} \
       --addHiggs \
       --text "${LMTEXT}" --unblind

