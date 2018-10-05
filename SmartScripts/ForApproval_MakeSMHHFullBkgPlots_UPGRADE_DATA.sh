#!/bin/bash

## 1) Make MaxLikelihood fit

LIMFOLDER=$1/CombinedCard_Node_SM_Prefit/

#combine --datacard ${LIMFOLDER}/CombinedCard_Node_SM${POINT}/hhbbgg_13TeV_DataCard.txt -M MaxLikelihoodFit --saveWorkspace --saveShapes --saveNormalization --X-rtd TMCSO_AdaptivePseudoAsimov=50 -n SMHHForBkgPlots

#mv mlfitSMHHForBkgPlots.root ${LIMFOLDER}
#mv higgsCombineSMHHForBkgPlots.MaxLikelihoodFit.mH120.root ${LIMFOLDER}
#mv MaxLikelihoodFitResult.root ${LIMFOLDER}

## 2) Make plots with fit

INFILE=${LIMFOLDER}/MaxLikelihoodFitResult.root

mkdir ${LIMFOLDER}/Background

echo $INFILE


OUTFILE_HM=FullBkgPlot_HM
OUTFILE_LM=FullBkgPlot_LM
#all four below to be taken from the datacard (signal rate)*33.49*0.0026
HMHP_NORM=13
HMMP_NORM=9.6
LMHP_NORM=10.2
LMMP_NORM=14
#these are multiplicative factors to make the signal show up
HMHP_FACT=1
HMMP_FACT=1
LMHP_FACT=1
LMMP_FACT=1

HMTEXT="#font[61]{pp#rightarrowHH#rightarrow#gamma#gammab#bar{b}}|High-mass region"
LMTEXT="#font[61]{pp#rightarrowHH#rightarrow#gamma#gammab#bar{b}}|Low-mass region"

python scripts/MakeFullBackgroundFit_Data.py -i ${INFILE} -o ${OUTFILE_HM} \
--signalNormalization ${HMHP_NORM} ${HMMP_NORM} \
--signalFactor ${HMHP_FACT} ${HMMP_FACT} --addHiggs \
--text "${HMTEXT}" --unblind

python scripts/MakeFullBackgroundFit_Data.py -i ${INFILE} -o ${OUTFILE_LM} \
--signalNormalization ${LMHP_NORM} ${LMMP_NORM} \
--signalFactor ${LMHP_FACT} ${LMMP_FACT} --addHiggs \
--text "${LMTEXT}" --unblind

mv FullBkgPlot* $LIMFOLDER/Background/
