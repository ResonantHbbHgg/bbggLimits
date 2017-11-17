#!/usr/bin/env python

import os.path
import getpass
username = getpass.getuser()
cwd = os.getcwd()

D_HM = "/afs/cern.ch/user/a/andrey/work/hh/NewHope_HighMass/"
D_LM = "/afs/cern.ch/user/a/andrey/work/hh/NewHope_LowMass"

FILE = "/LT_output_GluGluToHHTo2B2G_AllNodes.root"

case = [D_HM, D_LM]

org_batch = '''
#!/bin/bash
cd CWD
eval `scramv1 runtime -sh`
COMMAND

'''

import argparse

parser =  argparse.ArgumentParser(description='submit the tree maker to the batch')
parser.add_argument('-t', '--type', dest="scanType", default=None, type=str,
                    choices=['JHEP', 'KL', 'KLKT', 'manual'],
                    help = "Choose the type Trees to produce")

opt = parser.parse_args()


def makeMyTreeGood(inPath, kl=1, kt=1, cg=0, c2=0, c2g=0):
  pointStr = "_".join(['kl',str(float(kl)),'kt',str(float(kt)),'cg',str(float(cg)),'c2',str(float(c2)),'c2g',str(float(c2g))]).replace('.', 'p').replace('-', 'm')
  fname = 'LT_NR_Nodes_All_merged_'+pointStr+'.root'

  outputFile = inPath + '/' + fname
  if os.path.isfile(outputFile):
    print 'FILE ALREADY EXISTS, JUMPING..', kl, kt, cg, c2, c2g
    return
  else:
    print 'FILE DOESNT EXIST:', kl, kt, cg, c2, c2g

  command = 'python scripts/MakeARWTree.py -f ' + inPath+FILE + ' -o ' + outputFile + ' --kl ' + str(kl) + ' --kt ' + str(kt) + ' --cg ' + str(cg) + ' --c2 ' + str(c2) + ' --c2g ' + str(c2g) + '  '
  #     print command

  batch = org_batch.replace("COMMAND", command).replace("CWD", cwd)

  bFile = open('/tmp/'+username+'/batch_LT_'+pointStr+'.sh', "w+")
  bFile.write(batch)
  bFile.close()
  command = 'chmod a+rwx ' + '/tmp/'+username+'/batch_LT_'+pointStr+'.sh'
  os.system(command)
  command = "bsub -q 8nm -J batch_LT_" + pointStr+ " < /tmp/"+username+"/batch_LT_" + pointStr+'.sh'

  os.system(command)


if __name__ == "__main__":
  print "This is the __main__ part"

  from HiggsAnalysis.bbggLimits.DefineScans import *

  if 'JHEP' in opt.scanType:
    for cc in case:
      print cc
      for ii in range(0, len(klJHEP)):
        kl = klJHEP[ii]
        kt = ktJHEP[ii]
        cg = cgJHEP[ii]
        c2 = c2JHEP[ii]
        c2g =c2gJHEP[ii]

        makeMyTreeGood(cc, kl,kt,cg,c2,c2g)

  if 'KL' in opt.scanType:
    for cc in case:
      print 'Making trees for KL scan, ', cc
      for kl in scan_kl['kl']:
        print '\t kl = ', kl
        kt = 1.0
        cg = 0.0
        c2 = 0.0
        c2g = 0.0
        makeMyTreeGood(cc, kl,kt,cg,c2,c2g)
        
  if 'KLKT' in opt.scanType:
    for cc in case:
      print 'Making trees for KLKT scan, ', cc
      for kl in scan_2d['kl']:
        for kt in scan_2d['kt']:
          print '\t kl = ', kl, '\t kt = ', kt
          cg = 0.0
          c2 = 0.0
          c2g = 0.0
          makeMyTreeGood(cc, kl,kt,cg,c2,c2g)
        
  if 'manual' in opt.scanType:
    # Here make trees for specific points
    for cc in case:
      kl = 1.0
      kt = 1.0
      cg = 0.0
      c2 = 0.0
      c2g = 0.0
      
      makeMyTreeGood(cc, kl,kt,cg,c2,c2g)
