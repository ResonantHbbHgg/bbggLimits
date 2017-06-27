from ROOT import *
import os,sys
import getpass
username = getpass.getuser()

#Low Mass leading jet b-tag (loose)
LMLJBTC=0.55
#Low Mass subleading jet b-tag (medium)
LMSJBTC=0.55

if os.path.isdir("/tmp/"+username+"/eos/") == False :
  print "Mounting eos under /tmp/"+username+"/eos/ ..."
  command = 'mkdir /tmp/'+username+'/eos/'
  os.system(command)
  command = "eosmount /tmp/"+username+"/eos/"

L_BACKGROUND='/tmp/'+username+'/eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/May2_Mjj70to190_NewCatMVA/EGML_Background_Mjj70_NewMVA/Hadd/'
L_SIGNAL='/tmp/'+username+'/eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/May2_Mjj70to190_NewCatMVA/EGML_Signal_GEN/Hadd/'
L_DATA='/tmp/'+username+'/eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/May2_Mjj70to190_NewCatMVA/EGML_Data_Mjj70_NewMVA/Hadd/'

bashFile = '''
#!/bin/bash

cd HERE

eval `scramv1 runtime -sh`

python scripts/RunnerOfLimitsCustom_mjj70.py --dirname L_OUTDIR

'''

HERE = os.environ['PWD']

for M1 in ['600']:
  print M1

  for M2 in ['600']:
    print M2

    for C1 in ['985']:
      print C1

      for C0 in ['970']:
        print C0

        for MM in ['350']:

          L_OUTDIR = 'LT_'+str(MM)+'_HMHPC_'+str(C0)+'_HMMPC_'+str(M1)+'_LMHPC_'+str(C1)+'_LMMPC_'+str(M2)

          command = 'python scripts/makeAllTrees.py -x nonres -d ' + L_DATA + ' -s ' + L_SIGNAL + ' -f ' + L_OUTDIR + ' --doCatMVA --MVAHMC0 0.' + C0 + ' --MVAHMC1 0.' + M1 + ' --MVALMC0 0.' + C1 + ' --MVALMC1 0.' + M2 + ' --massNR ' + MM + ' --LMLJBTC ' + str(LMLJBTC) + ' --LMSJBTC ' + str(LMSJBTC) + ' '
          print command
          os.system(command)
          break

          command = 'python scripts/makeAllTrees.py -x nonres -s ' + L_BACKGROUND + ' -d 0 -f ' + L_OUTDIR + ' --doCatMVA --MVAHMC0 0.' + C0 + ' --MVAHMC1 0.' + M1 + ' --MVALMC0 0.' + C1 + ' --MVALMC1 0.' + M2 + ' --massNR ' + MM + ' --doSMHiggs --LMLJBTC ' + str(LMLJBTC) + ' --LMSJBTC ' + str(LMSJBTC) + ' --genDiPhotonFilter ' + ' '
          print command
          os.system(command)


          command = 'hadd ' + L_OUTDIR + '_HighMass/LT_output_bbHToGG_M-125_13TeV_amcatnlo.root ' + L_OUTDIR + '_HighMass/LT_output_bbHToGG_M-125_4FS_y*'
          os.system(command)
          command = 'hadd ' + L_OUTDIR + '_HighMass/LT_output_GluGluToHHTo2B2G_AllNodes.root ' + L_OUTDIR + '_HighMass/LT_output_GluGluToHHTo2B2G_nod*'
          os.system(command)
          command = 'hadd ' + L_OUTDIR + '_LowMass/LT_output_bbHToGG_M-125_13TeV_amcatnlo.root ' + L_OUTDIR + '_LowMass/LT_output_bbHToGG_M-125_4FS_y*'
          os.system(command)
          command = 'hadd ' + L_OUTDIR + '_LowMass/LT_output_GluGluToHHTo2B2G_AllNodes.root ' + L_OUTDIR + '_LowMass/LT_output_GluGluToHHTo2B2G_nod*'
          os.system(command)

'''
NOTE: UNCOMMENT THIS SECTION IF YOU ALSO WANT TO RUN THE LIMITS ON THE V0 BENCHMARKS WHILE MAKING THE LIMIT TREES
          bashFile = bashFile.replace('HERE', HERE).replace('L_OUTDIR', L_OUTDIR)
          bFile = open('batch_'+L_OUTDIR+'.sh', "w+")
          bFile.write(bashFile)
          bFile.close()
          command = 'chmod a+rwx ' + 'batch_'+L_OUTDIR+'.sh'
          os.system(command)
          command = "bsub -q 1nd  -J batch_" + L_OUTDIR  + " < batch_" + L_OUTDIR+'.sh'
          os.system(command)
'''
