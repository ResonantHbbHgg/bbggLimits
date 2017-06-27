from ROOT import *
import os,sys
import getpass
from HiggsAnalysis.bbggLimits.DefineScans import *
username = getpass.getuser()

org_bashFile = '''
#!/bin/bash

cd HERE

eval `scramv1 runtime -sh`

if [ ! -d "L_OUTLOC/L_OUTDIR_HighMass" ]; then
  echo 'MAKING DIRS... L_OUTLOC/ and L_OUTLOC/L_OUTDIR_HighMass and L_OUTLOC/L_OUTDIR_LowMass'
  mkdir L_OUTLOC/
  mkdir L_OUTLOC/L_OUTDIR_HighMass
  mkdir L_OUTLOC/L_OUTDIR_LowMass
fi

if [ ! -d "L_OUTLOC" ]; then
  echo "COULDNT CREATE DIR!!!! L_OUTLOC "
  exit 1
fi
if [ ! -d "L_OUTLOC/L_OUTDIR_HighMass" ]; then
  echo "COULDNT CREATE DIR!!!! L_OUTLOC/L_OUTDIR_HighMass"
  exit 1
fi
if [ ! -d "L_OUTLOC/L_OUTDIR_LowMass" ]; then
  echo "COULDNT CREATE DIR!!!! L_OUTLOC/L_OUTDIR_LowMass"
  exit 1
fi

if [ ! -f L_OUTLOC/L_OUTDIR_HighMass/ARWFILE ]; then
  echo 'Copying EOSLOC/L_OUTDIR_HighMass/ARWFILE to L_OUTLOC/L_OUTDIR_HighMass'
  xrdcp EOSLOC/L_OUTDIR_HighMass/ARWFILE L_OUTLOC/L_OUTDIR_HighMass
  xrdcp EOSLOC/L_OUTDIR_LowMass/ARWFILE L_OUTLOC/L_OUTDIR_LowMass
#  ls L_OUTLOC/L_OUTDIR_LowMass
#  xrdcp EOSLOC/L_OUTDIR_HighMass/ARWFILE L_OUTLOC/L_OUTDIR_HighMass
#  xrdcp EOSLOC/L_OUTDIR_LowMass/ARWFILE L_OUTLOC/L_OUTDIR_LowMass
fi

if [ ! -f L_OUTLOC/L_OUTDIR_HighMass/ARWFILE ]; then
  echo "COULDNT GET SIGNAL FILE FROM EOS!!!!!! L_OUTLOC/L_OUTDIR_HighMass/ARWFILE"
  exit 1
fi

if [ ! -f L_OUTLOC/L_OUTDIR_HighMass/LT_DoubleEG.root ]; then

  xrdcp EOSLOC/L_OUTDIR_HighMass/LT_DoubleEG.root L_OUTLOC/L_OUTDIR_HighMass ;

  xrdcp EOSLOC/L_OUTDIR_HighMass/LT_output_GluGluHToGG_M-125_13TeV_powheg_pythia8.root L_OUTLOC/L_OUTDIR_HighMass ;

  xrdcp EOSLOC/L_OUTDIR_HighMass/LT_output_VBFHToGG_M-125_13TeV_powheg_pythia8.root L_OUTLOC/L_OUTDIR_HighMass ;

  xrdcp EOSLOC/L_OUTDIR_HighMass/LT_output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root L_OUTLOC/L_OUTDIR_HighMass ;

  xrdcp EOSLOC/L_OUTDIR_HighMass/LT_output_bbHToGG_M-125_13TeV_amcatnlo.root L_OUTLOC/L_OUTDIR_HighMass ;

  xrdcp EOSLOC/L_OUTDIR_HighMass/LT_output_bbHToGG_M-125_4FS_yb2_13TeV_amcatnlo.root L_OUTLOC/L_OUTDIR_HighMass ;

  xrdcp EOSLOC/L_OUTDIR_HighMass/LT_output_bbHToGG_M-125_4FS_ybyt_13TeV_amcatnlo.root L_OUTLOC/L_OUTDIR_HighMass ;

  xrdcp EOSLOC/L_OUTDIR_HighMass/LT_output_ttHToGG_M125_13TeV_powheg_pythia8_v2.root L_OUTLOC/L_OUTDIR_HighMass ;

  xrdcp EOSLOC/L_OUTDIR_LowMass/LT_DoubleEG.root L_OUTLOC/L_OUTDIR_LowMass ;

  xrdcp EOSLOC/L_OUTDIR_LowMass/LT_output_GluGluHToGG_M-125_13TeV_powheg_pythia8.root L_OUTLOC/L_OUTDIR_LowMass ;

  xrdcp EOSLOC/L_OUTDIR_LowMass/LT_output_VBFHToGG_M-125_13TeV_powheg_pythia8.root L_OUTLOC/L_OUTDIR_LowMass ;

  xrdcp EOSLOC/L_OUTDIR_LowMass/LT_output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root L_OUTLOC/L_OUTDIR_LowMass ;

  xrdcp EOSLOC/L_OUTDIR_LowMass/LT_output_bbHToGG_M-125_13TeV_amcatnlo.root L_OUTLOC/L_OUTDIR_LowMass ;

  xrdcp EOSLOC/L_OUTDIR_LowMass/LT_output_bbHToGG_M-125_4FS_yb2_13TeV_amcatnlo.root L_OUTLOC/L_OUTDIR_LowMass ;

  xrdcp EOSLOC/L_OUTDIR_LowMass/LT_output_bbHToGG_M-125_4FS_ybyt_13TeV_amcatnlo.root L_OUTLOC/L_OUTDIR_LowMass ;

  xrdcp EOSLOC/L_OUTDIR_LowMass/LT_output_ttHToGG_M125_13TeV_powheg_pythia8_v2.root L_OUTLOC/L_OUTDIR_LowMass ;

fi

echo 'CONTENT of LOWMASS:'
ls L_OUTLOC/L_OUTDIR_LowMass

echo 'CONTENT of HIGHMASS:'
ls L_OUTLOC/L_OUTDIR_HighMass

if [ ! -f L_OUTLOC/L_OUTDIR_HighMass/LT_DoubleEG.root ]; then
  echo "COULDNT GET DATA FILE FROM EOS!!!!!!"
  exit 1
fi

python HERE/scripts/RunnerOfLimitsCustom_mjj70.py --dirname L_OUTDIR --dirloc L_OUTLOC --extra 'EXTRA' --jsonName JSONNAME

#exit 1
'''

HERE = os.environ['PWD']

scan = scan_2d

#FOLDER = 'LT_MJJ70NewMVA_350_HMHPC_970_HMMPC_600_LMHPC_985_LMMPC_600_allNodes'
FOLDER = 'LT_350_HMHPC_970_HMMPC_600_LMHPC_985_LMMPC_600'
EOSLOC = 'root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/resonant_HH/RunII/LimitTrees/2016/'
L_OUTLOC_org = '/tmp/'+username+'/'

counter = 0
for kl in scan['kl']:
 for kt in scan['kt']:
  for cg in scan['cg']:
   for c2 in scan['c2']:
    for c2g in scan['c2g']:
     print counter
     counter += 1
     print kl, kt, cg, c2, c2g
#     if kl > 0.99 and kl < 1.01 and kt > 0.99 and 1.01 and cg == 0 and c2g == 0 and c2 == 0: continue
#     if kl > 9.99 and kl < 10.01 and kt > 0.99 and 1.01 and cg == 0 and c2g == 0 and c2 == 0: continue
#     if kl < -9.99 and kl > -10.01 and kt > 0.99 and 1.01 and cg == 0 and c2g == 0 and c2 == 0: continue

     extra = ' --analyticalRW '
     extra += ' --kl ' + str(kl)
     extra += ' --kt ' + str(kt)
     extra += ' --cg ' + str(cg)
     extra += ' --c2 ' + str(c2)
     extra += ' --c2g ' + str(c2g)
     pointStr = 'kl'+str(kl).replace('.', 'p').replace('-', 'm')+'_kt'+str(kt).replace('.', 'p').replace('-', 'm')+'_cg'+str(cg).replace('.', 'p').replace('-', 'm')+'_c2'+str(c2).replace('.', 'p').replace('-', 'm')+'_c2g'+str(c2g).replace('.', 'p').replace('-', 'm')
     extra += ' --extraLabel ' + pointStr + ' '

     L_OUTLOC = L_OUTLOC_org + '/' + pointStr + '/'

     ARWFILE = "/LT_NR_Nodes_All_merged_kl_"+str(kl).replace(".", "p")+"_kt_"+str(kt).replace(".", "p")+"_cg_"+str(cg).replace(".", "p")+"_c2_"+str(c2).replace(".", "p")+"_c2g_"+str(c2g).replace(".", "p")+".root"

     JSONNAME = '/tmp/'+username+'/json_'+pointStr+'.json'

     bashFile = org_bashFile.replace('HERE', HERE).replace('L_OUTDIR', FOLDER).replace('EXTRA', extra).replace('L_OUTLOC', L_OUTLOC).replace("EOSLOC", EOSLOC).replace("ARWFILE", ARWFILE).replace("JSONNAME", JSONNAME)

     bFile = open('/tmp/'+username+'/batch_'+pointStr+'.sh', "w+")
     bFile.write(bashFile)
     bFile.close()
     command = 'chmod a+rwx ' + '/tmp/'+username+'/batch_'+pointStr+'.sh'
     os.system(command)
#     command = 'source ' + '/tmp/'+username+'/batch_'+pointStr+'.sh'
#     os.system(command)
     command = "bsub -q 1nh -o /tmp/"+username+"/"+pointStr+".log -J batch_" + pointStr  + " < /tmp/"+username+"/batch_" + pointStr + '.sh'
#     command = "bsub -q 1nh  -J batch_" + pointStr  + " < /tmp/"+username+"/batch_" + pointStr + '.sh'
     print command
     os.system(command)
#     sys.exit()

