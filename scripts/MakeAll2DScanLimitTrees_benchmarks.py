import time
import argparse, os
import getpass
username = getpass.getuser()

import os.path


parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument("-f", "--file", dest="f", type=str)
parser.add_argument("-o", "--outfolder", dest="out", type=str)
opt = parser.parse_args()

org_batch = '''
#!/bin/bash

cd /afs/cern.ch/work/r/rateixei/work/DiHiggs/bbggLimits_May8/CMSSW_7_4_7/src/HiggsAnalysis/bbggLimits/

eval `scramv1 runtime -sh`

if [ ! -d /tmp/rateixei/eos/ ] ; then 
mkdir /tmp/rateixei/eos/
fi

if [ ! -d /tmp/rateixei/eos/cms/ ] ; then
/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select -b fuse mount /tmp/rateixei/eos/
fi

COMMAND

#exit 1
'''

scan_2d = {
'kl': [float(i)*2.0 for i in range(-10,11)], #-20, -17.5, -15, -12.5, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20],
'kt': [float(i)/4.0 for i in range(-10,11)],#-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5],
'cg': [0.0],
'c2': [0.0],
'c2g': [0.0]
}

scan_2d_fix = {
'kl': [float(i)*2.0 for i in range(-10,11)], #-20, -17.5, -15, -12.5, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20],
'kt': [0.0],#-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5],
'cg': [0.0],
'c2': [0.0],
'c2g': [0.0]
}

scan_kl = {
'kl': [float(i)/10. for i in range(-200,201)],
'kt': [1.0],
'cg': [0.0],
'c2': [0.0],
'c2g': [0.0]
}

scan = scan_2d

counter = 0

klJHEP=[1.0,  7.5,  1.0,  1.0,  -3.5, 1.0, 2.4, 5.0, 15.0, 1.0, 10.0, 2.4, 15.0]
ktJHEP=[1.0,  1.0,  1.0,  1.0,  1.5,  1.0, 1.0, 1.0, 1.0,  1.0, 1.5,  1.0, 1.0]
c2JHEP=[0.0,  -1.0, 0.5, -1.5, -3.0,  0.0, 0.0, 0.0, 0.0,  1.0, -1.0, 0.0, 1.0]
cgJHEP=[0.0,  0.0, -0.8,  0.0, 0.0,   0.8, 0.2, 0.2, -1.0, -0.6, 0.0, 1.0, 0.0]
c2gJHEP=[0.0, 0.0, 0.6, -0.8, 0.0, -1.0, -0.2,-0.2,  1.0,  0.6, 0.0, -1.0, 0.0]

for ii in xrange(0,13):
# for ikt in scan['kt']:
#  for cg in scan['cg']:
#   for c2 in scan['c2']:
#    for c2g in scan['c2g']:

     kl = klJHEP[ii]
     kt = ktJHEP[ii]
     c2 = c2JHEP[ii]
     cg = cgJHEP[ii]
     c2g = c2gJHEP[ii]
#     print kl, kt, cg, c2, c2g

     fname = 'LT_NR_Nodes_All_merged_kl_'+ str(kl).replace('.', 'p') + '_kt_' + str(kt).replace('.', 'p') + '_cg_' + str(cg).replace('.', 'p') + '_c2_' + str(c2).replace('.', 'p') + '_c2g_' + str(c2g).replace('.', 'p') + '.root'
#     my_file = Path(opt.out + '/' + fname)
#     if my_file.is_file():
#     print opt.out + '/' + fname
     if os.path.isfile(opt.out + '/' + fname):
       print 'FILE ALREADY EXISTS, JUMPING..', kl, kt, cg, c2, c2g
       continue
     else:
       print 'FILE DOESNT EXIST!!', kl, kt, cg, c2, c2g

     command = 'python scripts/MakeARWTree.py -f ' + opt.f + ' -o ' + opt.out + ' --kl ' + str(kl) + ' --kt ' + str(kt) + ' --cg ' + str(cg) + ' --c2 ' + str(c2) + ' --c2g ' + str(c2g) + '  '
#     print command
     batch = org_batch.replace("COMMAND", command)
     pointStr = 'kl'+str(kl).replace('.', 'p').replace('-', 'm')+'_kt'+str(kt).replace('.', 'p').replace('-', 'm')+'_cg'+str(cg).replace('.', 'p').replace('-', 'm')+'_c2'+str(c2).replace('.', 'p').replace('-', 'm')+'_c2g'+str(c2g).replace('.', 'p').replace('-', 'm')
 
     bFile = open('/tmp/'+username+'/batch_LT_'+pointStr+'.sh', "w+")
     bFile.write(batch)
     bFile.close()
     command = 'chmod a+rwx ' + '/tmp/'+username+'/batch_LT_'+pointStr+'.sh'
     os.system(command)
     command = "bsub -q 8nm -o /tmp/"+username+"/LT_"+pointStr+".log -J batch_LT_" + pointStr  + " < /tmp/"+username+"/batch_LT_" + pointStr + '.sh'
#     command = "bsub -q 8nm -J batch_LT_" + pointStr  + " < /tmp/"+username+"/batch_LT_" + pointStr + '.sh'
     os.system(command)

