import time
import argparse, os, sys
import getpass
import os.path
from HiggsAnalysis.bbggLimits.DefineScans import *
username = getpass.getuser()
cwd = os.getcwd()

'''
NOTE:
Unfortunately there's something wrong with this script.
During the batch step, several files will fail because you won't be able to mount EOS on that specific lxbatch node.
Therefore, you will need to run this script several times until no file is missing.
The script will only run again the files that don't exist in the folder, so you won't run everything again if you run the script twice (only the missing files).
'''

if os.path.isdir("/tmp/"+username+"/eos/") == False :
  print "Mounting eos under /tmp/"+username+"/eos/ ..."
  command = 'mkdir /tmp/'+username+'/eos/'
  os.system(command)
  command = "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select -b fuse mount /tmp/"+username+"/eos/"
  os.system(command)


D_HM = "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/resonant_HH/RunII/LimitTrees/2016/LT_350_HMHPC_970_HMMPC_600_LMHPC_985_LMMPC_600_HighMass/"
D_LM = "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/resonant_HH/RunII/LimitTrees/2016/LT_350_HMHPC_970_HMMPC_600_LMHPC_985_LMMPC_600_LowMass/"
FILE = "/LT_output_GluGluToHHTo2B2G_AllNodes.root"

case = [D_HM, D_LM]

org_batch = '''
#!/bin/bash

cd CWD

eval `scramv1 runtime -sh`

if [ ! -d /tmp/USER/POINT/ ] ; then 
mkdir /tmp/USER/POINT/
/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select cp FOLD/LT_output_GluGluToHHTo2B2G_AllNodes.root /tmp/USER/POINT/
fi

if [ ! -f /tmp/USER/POINT/LT_output_GluGluToHHTo2B2G_AllNodes.root ] ; then
echo "Couldnt download LT_output_GluGluToHHTo2B2G_AllNodes.root"
exit 1
fi

COMMAND

/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select cp /tmp/USER/POINT/FILE FOLD

#exit 1
'''

scans = [ scan_2d, scan_kl ]
mkount = 0
for scan in scans:
 for cc in case:
  print cc
  for kl in scan['kl']:
   for kt in scan['kt']:
    for cg in scan['cg']:
     for c2 in scan['c2']:
      for c2g in scan['c2g']:
  #     print kl, kt, cg, c2, c2g
  
       fname = 'LT_NR_Nodes_All_merged_kl_'+ str(kl).replace('.', 'p') + '_kt_' + str(kt).replace('.', 'p') + '_cg_' + str(cg).replace('.', 'p') + '_c2_' + str(c2).replace('.', 'p') + '_c2g_' + str(c2g).replace('.', 'p') + '.root'
  #     my_file = Path(opt.out + '/' + fname)
  #     if my_file.is_file():
  #     print opt.out + '/' + fname
       locPath = "/tmp/"+username+"/eos/" + cc.replace("root://eoscms.cern.ch//eos/","")
       if os.path.isfile(locPath + '/' + fname):
         print 'FILE ALREADY EXISTS, JUMPING..', kl, kt, cg, c2, c2g
         continue
       else:
         print 'FILE DOESNT EXIST!!', kl, kt, cg, c2, c2g
 
       pointStr = 'kl'+str(kl).replace('.', 'p').replace('-', 'm')+'_kt'+str(kt).replace('.', 'p').replace('-', 'm')+'_cg'+str(cg).replace('.', 'p').replace('-', 'm')+'_c2'+str(c2).replace('.', 'p').replace('-', 'm')+'_c2g'+str(c2g).replace('.', 'p').replace('-', 'm')

       dtemp = pointStr+str(mkount)

       command = 'python scripts/MakeARWTree.py -f ' + cc+FILE + ' -o /tmp/' + username + '/'+ dtemp + ' --kl ' + str(kl) + ' --kt ' + str(kt) + ' --cg ' + str(cg) + ' --c2 ' + str(c2) + ' --c2g ' + str(c2g) + '  '
  #     print command


       batch = org_batch.replace("COMMAND", command).replace("CWD", cwd).replace("USER", username).replace('POINT', dtemp ).replace("FOLD", cc).replace("FILE", fname)
   
       bFile = open('/tmp/'+username+'/batch_LT_'+pointStr+'_'+str(mkount)+'.sh', "w+")
       bFile.write(batch)
       bFile.close()
       command = 'chmod a+rwx ' + '/tmp/'+username+'/batch_LT_'+pointStr+'_'+str(mkount)+'.sh'
       os.system(command)
       command = "bsub -q 8nm -o /tmp/"+username+"/LT_"+pointStr+".log -J batch_LT_" + pointStr+'_'+str(mkount)  + " < /tmp/"+username+"/batch_LT_" + pointStr+'_'+str(mkount) + '.sh'
  #     command = "bsub -q 8nm -J batch_LT_" + pointStr  + " < /tmp/"+username+"/batch_LT_" + pointStr + '.sh'
       os.system(command)
       mkount += 1
