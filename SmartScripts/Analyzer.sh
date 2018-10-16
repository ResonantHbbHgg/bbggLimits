#!/bin/bash

DIR=$1/CombinedCard_Node_SM
echo $DIR

text2workspace.py $DIR/hhbbgg_13TeV_DataCard.txt $DIR/hhbbgg_13TeV_DataCard.root -m 125 --X-nuisance-group-function 'theory' '0.5'



cd $DIR

combine -M MaxLikelihoodFit -t -1 --expectSignal 1 -d hhbbgg_13TeV_DataCard.root --saveWorkspace --saveShapes --saveNormalization &>MaxLikelihood.txt

combine -M ProfileLikelihood -d hhbbgg_13TeV_DataCard.root --significance -t -1 --expectSignal 1 -m 125 -n SM_13TeV_3ab &>ProfileLikelihood.txt

combine -M ProfileLikelihood -d hhbbgg_13TeV_DataCard.root --significance -t -1 --expectSignal 1 -m 125 -n SM_13TeV_3ab

combine -M Asymptotic -d hhbbgg_13TeV_DataCard.root -t -1  --run blind --expectSignal 1 -m 125 -n SM_13TeV_3ab -S 0 &>Limit_stat.txt

combine -M Asymptotic -d hhbbgg_13TeV_DataCard.root -t -1  --run blind --expectSignal 1 -m 125 -n SM_13TeV_3ab &>Limit.txt


cd -
