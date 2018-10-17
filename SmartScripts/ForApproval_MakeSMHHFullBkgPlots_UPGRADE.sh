#!/bin/bash

## 1) Make MaxLikelihood fit

LIMFOLDER=$1/CombinedCard_Node_SM

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
#these are multiplicative factors to make the signal show up
HMHP_FACT=1
HMMP_FACT=1
LMHP_FACT=1
LMMP_FACT=1

HMTEXT="#font[61]{pp#rightarrowHH#rightarrow#gamma#gammab#bar{b}}|High-mass region"
LMTEXT="#font[61]{pp#rightarrowHH#rightarrow#gamma#gammab#bar{b}}|Low-mass region"

python scripts/MakeFullBackgroundFit.py -i ${INFILE} -o ${OUTFILE_HM} \
--signalFactor ${HMHP_FACT} ${HMMP_FACT} --addHiggs \
--text "${HMTEXT}" --unblind

python scripts/MakeFullBackgroundFit.py -i ${INFILE} -o ${OUTFILE_LM} \
--signalFactor ${LMHP_FACT} ${LMMP_FACT} --addHiggs \
--text "${LMTEXT}" --unblind


mv FullBkgPlot* $LIMFOLDER/Background/
