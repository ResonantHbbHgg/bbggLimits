#!/bin/bash

echo '=== Preparing datacards in root format ==='

DIR=$1/CombinedCard_Node_SM
echo $DIR

text2workspace.py $DIR/hhbbgg_13TeV_DataCard.txt $DIR/hhbbgg_13TeV_DataCard.root -m 125 --X-nuisance-group-function 'theory' '0.5'

DIR_Prefit=$1/CombinedCard_Node_SM_Prefit
echo ${DIR_Prefit}

cp -r ${DIR} ${DIR_Prefit}

text2workspace.py $DIR_Prefit/hhbbgg_13TeV_DataCard.txt $DIR_Prefit/hhbbgg_13TeV_DataCard.root -m 125 --X-nuisance-group-function 'theory' '0.5'

cd $DIR

printf '\n=== Start fitting ==\n\n'

combine -M MaxLikelihoodFit -t -1 --expectSignal 1 -d hhbbgg_13TeV_DataCard.root -S 0 &>MaxLikelihood_stat.txt


echo '== 1.1) Finished Max Likelihood Fit Stat'

tail MaxLikelihood_stat.txt

printf '\n=====\n\n'

combine -M MaxLikelihoodFit -t -1 --expectSignal 1 -d hhbbgg_13TeV_DataCard.root --saveWorkspace --saveShapes --saveNormalization &>MaxLikelihood.txt

echo '== 1.2) Finished Max Likelihood Fit'

tail -36 MaxLikelihood.txt

printf '\n=====\n\n'

combine -M ProfileLikelihood -d hhbbgg_13TeV_DataCard.root --significance -t -1 --expectSignal 1 -m 125 -n SM_13TeV_3ab -S 0 &>ProfileLikelihood_stat.txt

echo '== 2.1) Finished Profile Likelihood Fit stat'

head ProfileLikelihood_stat.txt

printf '\n=====\n\n'

combine -M ProfileLikelihood -d hhbbgg_13TeV_DataCard.root --significance -t -1 --expectSignal 1 -m 125 -n SM_13TeV_3ab &>ProfileLikelihood.txt

echo '== 2.2) Finished Profile Likelihood Fit'

head ProfileLikelihood.txt

printf '\n=====\n\n'

combine -M Asymptotic -d hhbbgg_13TeV_DataCard.root -t -1  --run blind --expectSignal 1 -m 125 -n SM_13TeV_3ab -S 0 &>Limit_stat.txt

echo '== 3.1) Finished stat only limits'


head Limit_stat.txt

printf '\n=====\n\n'

combine -M Asymptotic -d hhbbgg_13TeV_DataCard.root -t -1  --run blind --expectSignal 1 -m 125 -n SM_13TeV_3ab &>Limit.txt



echo '== 3.2) Finished full limits'

head Limit.txt

printf '\n=====\n\n'

cd -

echo ' == 4.1) Making prefit'

cd ${DIR_Prefit}

combine -M MaxLikelihoodFit -d hhbbgg_13TeV_DataCard.root --saveWorkspace --saveShapes --saveNormalization --setPhysicsModelParameters yield_norm=1 -S 0

printf '\n== END ===\n\n'

cd -