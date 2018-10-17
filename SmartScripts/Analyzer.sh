#!/bin/bash

DIR=$1/CombinedCard_Node_SM
echo $DIR

text2workspace.py $DIR/hhbbgg_13TeV_DataCard.txt $DIR/hhbbgg_13TeV_DataCard.root -m 125 --X-nuisance-group-function 'theory' '0.5'

cd $DIR

combine -M MaxLikelihoodFit -t -1 --expectSignal 1 -d hhbbgg_13TeV_DataCard.root -S 0 &>MaxLikelihood_stat.txt

echo 'Finished Max Likelihood Fit Stat'

tail MaxLikelihood_stat.txt

printf '\n=====\n\n'

combine -M MaxLikelihoodFit -t -1 --expectSignal 1 -d hhbbgg_13TeV_DataCard.root --saveWorkspace --saveShapes --saveNormalization &>MaxLikelihood.txt

echo 'Finished Max Likelihood Fit'

tail MaxLikelihood.txt

printf '\n=====\n\n'

combine -M ProfileLikelihood -d hhbbgg_13TeV_DataCard.root --significance -t -1 --expectSignal 1 -m 125 -n SM_13TeV_3ab -S 0 &>ProfileLikelihood_stat.txt

echo 'Finished Profile Likelihood Fit stat'

head ProfileLikelihood_stat.txt

printf '\n=====\n\n'

combine -M ProfileLikelihood -d hhbbgg_13TeV_DataCard.root --significance -t -1 --expectSignal 1 -m 125 -n SM_13TeV_3ab &>ProfileLikelihood.txt

echo 'Finished Profile Likelihood Fit'

head ProfileLikelihood.txt

printf '\n=====\n\n'

combine -M Asymptotic -d hhbbgg_13TeV_DataCard.root -t -1  --run blind --expectSignal 1 -m 125 -n SM_13TeV_3ab -S 0 &>Limit_stat.txt

echo 'Finished stat only limits'


head Limit_stat.txt

printf '\n=====\n\n'

combine -M Asymptotic -d hhbbgg_13TeV_DataCard.root -t -1  --run blind --expectSignal 1 -m 125 -n SM_13TeV_3ab &>Limit.txt



echo 'Finished full limits'

head Limit.txt

printf '\n=====\n\n'

cd -
