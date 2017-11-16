from ROOT import *
import os,sys
import getpass
from HiggsAnalysis.bbggLimits.DefineScans import *
username = getpass.getuser()

# This one id for Analytical reweighting:
org_bashFile = '''
#!/bin/bash

cd HERE

eval `scramv1 runtime -sh`

## This script is ugly, need to get rid of it:
python HERE/scripts/RunnerOfLimitsCustom_mjj70.py --dirname L_OUTDIR --dirloc L_OUTLOC --extra 'EXTRA' --jsonName JSONNAME

'''

HERE = os.environ['PWD']
FOLDER = 'MXFIX_LT_350_HMHPC_970_HMMPC_600_LMHPC_985_LMMPC_600'
L_OUTLOC_org = '/afs/cern.ch/user/a/andrey/work/hh/'

import argparse

parser =  argparse.ArgumentParser(description='submit the limit to the batch')
parser.add_argument('-t', '--type', dest="limtype", default=None, type=str, nargs='+',
                    choices=['JHEP', 'KL', 'KLKT', 'grid'],
                    help = "Choose the type of limits to run")
                    
opt = parser.parse_args()


def submitPoint(kl, kt, cg, c2, c2g):
  extra = ' --analyticalRW '
  extra += ' --kl ' + str(kl)
  extra += ' --kt ' + str(kt)
  extra += ' --cg ' + str(cg)
  extra += ' --c2 ' + str(c2)
  extra += ' --c2g ' + str(c2g)
  pointStr = "_".join(['kl'+str(kl), 'kt' + str(kt), 'cg'+ str(cg), 'c2' + str(c2), 'c2g' + str(c2g)]).replace('.', 'p').replace('-', 'm')

  extra += ' --extraLabel ' + pointStr + ' '

  L_OUTLOC = L_OUTLOC_org
  
  JSONNAME = '/tmp/'+username+'/json_'+pointStr+'.json'

  bashFile = org_bashFile.replace('HERE', HERE).replace('L_OUTDIR', FOLDER).replace('EXTRA', extra).replace('L_OUTLOC', L_OUTLOC).replace("JSONNAME", JSONNAME)

  bFile = open('/tmp/'+username+'/batch_'+pointStr+'.sh', "w+")
  bFile.write(bashFile)
  bFile.close()
  command = 'chmod a+rwx ' + '/tmp/'+username+'/batch_'+pointStr+'.sh'
  os.system(command)
  
  command = "bsub -q 1nh -J batch_" + pointStr  + " < /tmp/"+username+"/batch_" + pointStr + '.sh'

  print command
  os.system(command)
     

if __name__ == "__main__":
  print "This is the __main__ part"


  if "JHEP" in opt.limtype:
    counter = 0
    for ii in range(0, len(klJHEP)):
      kl = klJHEP[ii]
      kt = ktJHEP[ii]
      cg = cgJHEP[ii]
      c2 = c2JHEP[ii]
      c2g =c2gJHEP[ii]
      
      print counter
      counter += 1
      print kl, kt, cg, c2, c2g
      submitPoint(kl, kt, cg, c2, c2g)   


  if "KL" in opt.limtype:
    counter = 0
    for kl in [float(i)/(2.5) for i in range(-50,51)]:
      kt = 1.0
      cg = 0.0
      c2 = 0.0
      c2g = 0.0
            
      print counter
      counter += 1
      print kl, kt, cg, c2, c2g
      submitPoint(kl, kt, cg, c2, c2g)   

  if "grid" in opt.limtype:
    counter = 0
    for j in xrange(0,5):
      command = ' '.join(['bsub -q 1nh','./scripts/batch_job_for_grid.sh',HERE,'GridLimits', str(j)]) 
      print command
      os.system(command)
