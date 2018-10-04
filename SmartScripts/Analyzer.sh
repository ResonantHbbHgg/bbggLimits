OUT=0_Exp

text2workspace.py outDir_MILANO_v3_split${OUT}/CombinedCard_Node_SM/hhbbgg_13TeV_DataCard.txt outDir_MILANO_v3_split${OUT}/CombinedCard_Node_SM/hhbbgg_13TeV_DataCard.root -m 125 --X-nuisance-group-function 'theory' '0.5'

DIR=outDir_MILANO_v3_split${OUT}/CombinedCard_Node_SM/
echo $DIR
cd $DIR

combine -M MaxLikelihoodFit -t -1 --expectSignal 1 -d hhbbgg_13TeV_DataCard.root --saveWorkspace --saveShapes --saveNormalization &>MaxLikelihood.txt

combine -M ProfileLikelihood -d hhbbgg_13TeV_DataCard.root --significance -t -1 --expectSignal 1 -m 125 -n SM_13TeV_3ab &>ProfileLikelihood.txt

combine -M ProfileLikelihood -d hhbbgg_13TeV_DataCard.root --significance -t -1 --expectSignal 1 -m 125 -n SM_13TeV_3ab

cd -
