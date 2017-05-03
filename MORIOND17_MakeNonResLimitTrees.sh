#!/bin/bash

#high mass mpc
for M1 in 780 ;
do 

#low mass mpc
for M2 in 475 ;
do

#low mass hpc
for C1 in 980 ;
do 

#high mass hpc
for C0 in 990 ;
do

#Low Mass leading jet b-tag (loose)
LMLJBTC=0.56
#Low Mass subleading jet b-tag (medium)
LMSJBTC=0.86

#low mass/high mass boundary
for MM in 350 ;
do

#samples on EOS
    
    
    # L_SIGNAL=root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/April26_Mjj70to190/EGML_Signal_Mjj70to190/Hadd
    L_SIGNAL=./FlatTrees-May02-Sig
    L_BACKGROUND=root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/April26_Mjj70to190/EGML_Background_Mjj70to190/Hadd
    L_DATA=root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/April26_Mjj70to190/EGML_Data_Mjj70to190/Hadd

HERE=${PWD}

L_OUTDIR=LT_${MM}_HMHPC_${C0}_HMMPC_${M1}_LMHPC_${C1}_LMMPC_${M2}

# For Non-Res samples, including the weights:
python scripts/makeAllTrees.py -x nonres -d 0 -s ${L_SIGNAL} -f ${L_OUTDIR} --NRW --doCatMVA --MVAHMC0 0.${C0} --MVAHMC1 0.${M1} --MVALMC0 0.${C1} --MVALMC1 0.${M2} --massNR ${MM} --LMLJBTC ${LMLJBTC} --LMSJBTC ${LMSJBTC};

#python scripts/makeAllTrees.py -x nonres -d ${L_DATA} -s ${L_SIGNAL} -f ${L_OUTDIR} --doCatMVA --MVAHMC0 0.${C0} --MVAHMC1 0.${M1} --MVALMC0 0.${C1} --MVALMC1 0.${M2} --massNR ${MM} --onlySMHH --LMLJBTC ${LMLJBTC} --LMSJBTC ${LMSJBTC};

#python scripts/makeAllTrees.py -x nonres -s ${L_BACKGROUND} -d 0 -f ${L_OUTDIR} --doCatMVA --MVAHMC0 0.${C0} --MVAHMC1 0.${M1} --MVALMC0 0.${C1} --MVALMC1 0.${M2} --massNR ${MM} --doSMHiggs --LMLJBTC ${LMLJBTC} --LMSJBTC ${LMSJBTC} --genDiPhotonFilter ;


#cd ${L_OUTDIR}_HighMass
#hadd LT_output_bbHToGG_M-125_13TeV_amcatnlo.root LT_output_bbHToGG_M-125_4FS_y*
#cd ..
#cd ${L_OUTDIR}_LowMass
#hadd LT_output_bbHToGG_M-125_13TeV_amcatnlo.root LT_output_bbHToGG_M-125_4FS_y*
#cd ..


#python scripts/RunnerOfLimitsCustom.py --dirname ${L_OUTDIR}


done
done
done
done
done
