from ROOT import *
import os,sys

#Low Mass leading jet b-tag (loose)
LMLJBTC=0.00
#Low Mass subleading jet b-tag (medium)
LMSJBTC=0.00

#L_SIGNAL='/afs/cern.ch/work/r/rateixei/work/DiHiggs/bbggTools_flashgg_tag-Moriond17-v10/CMSSW_8_0_26_patch1/src/flashgg/bbggTools/test/RunJobs/EGML_Signal_GEN/Hadd/'
L_BACKGROUND='/tmp/rateixei/eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/May2_Mjj70to190_NewCatMVA/EGML_Background_Mjj70_NewMVA/Hadd/'
L_DATA='/tmp/rateixei/eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/May2_Mjj70to190_NewCatMVA/EGML_Data_Mjj70_NewMVA/Hadd/'
L_SIGNAL='/tmp/rateixei/eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/May2_Mjj70to190_NewCatMVA/EGML_Signal_Mjj70_NewMVA/Hadd/'

bashFile_org = '''
#!/bin/bash

cd HERE

eval `scramv1 runtime -sh`

python scripts/RunnerOfLimitsCustom_mjj70.py --dirname L_OUTDIR --mass RMASS --jsonName JSONNAME --resType RTYPE

'''

HERE = os.environ['PWD']

LM_Cat1 = ['600', '650', '700', '750', '800', '850']

MM = '500'

LowMasses = ['250', '260', '270', '280', '300', '320', '350', '400', '450', '500', '550', '600']
HighMasses = ['500', '550', '600', '650', '700', '750', '800', '900']

resTypes = ['BulkGraviton', 'Radion']
#resType = 'Radion'
#--resType", dest="resType", choices=['Radion', 'BulkGraviton']

HighMass = {'M1': ['000'], 'C1': ['500'], 'Mass': HighMasses}
LowMass = {'M1': ['700'], 'C1': ['960'], 'Mass': LowMasses}

Case = [LowMass, HighMass]

for resType in resTypes:
 for cc in Case:
  for M1 in cc['M1']:
    print M1
    M2=M1
    print M2

    for C1 in cc['C1']:
        print C1
        C0=C1
        print C0

        for Mass in cc['Mass']:
          print MM

          L_OUTDIR = 'LT_Res_'+str(MM)+'_HMHPC_'+str(C0)+'_HMMPC_'+str(M1)+'_LMHPC_'+str(C1)+'_LMMPC_'+str(M2)

          command = 'python scripts/makeAllTrees.py -x res -d '+ L_DATA +' -s ' + L_SIGNAL + ' -f ' + L_OUTDIR + '_ --doCatMVA --MVAHMC0 0.' + C0 + ' --MVAHMC1 0.' + M1 + ' --MVALMC0 0.' + C1 + ' --MVALMC1 0.' + M2 + ' --massNR ' + MM + ' --LMLJBTC ' + str(LMLJBTC) + ' --LMSJBTC ' + str(LMSJBTC) + ' --resMass ' + Mass + ' --resType ' + resType

          print command
          os.system(command)

          jsonname = 'json_'+L_OUTDIR+'_Mass_'+Mass+'_'+resType+'.json'

          bashFile = bashFile_org.replace('HERE', HERE).replace('L_OUTDIR', L_OUTDIR).replace('RMASS', Mass).replace('JSONNAME', jsonname).replace('RTYPE', resType)
          bFile = open('batch_'+L_OUTDIR+'_mass_'+Mass+'_'+resType+'.sh', "w+")
          bFile.write(bashFile)
          bFile.close()
          command = 'chmod a+rwx ' + 'batch_'+L_OUTDIR+'_mass_'+Mass+'_'+resType+'.sh'
          os.system(command)
          command = "bsub -q 1nd  -J batch_" + L_OUTDIR+'_mass_'+Mass+'_'+resType  + " < batch_" + L_OUTDIR+'_mass_'+Mass+'_'+resType+'.sh'
          os.system(command)


