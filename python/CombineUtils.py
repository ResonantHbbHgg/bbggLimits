from ROOT import *
import os,sys,json,time,re
import logging
from shutil import copy
from pprint import pformat
import getpass
username = getpass.getuser()


def runCombine(inDir, doBlind, log, combineOpt = 1, Combinelxbatch = 0, Label = None, scaleSingleHiggs=False):
  log.info('Running combine tool.  Dir: %s Blinded: %r', inDir, doBlind)
  log.debug('inDir should be the immediate directory where the card is located')
  if Combinelxbatch:
    runCombineOnLXBatch(inDir, doBlind, log, combineOpt, Label)
    return

  if doBlind and combineOpt!=3:
    # In HybridNew this option does not work
    blinded = "--run blind"
  else:
    blinded = ''

  if combineOpt==1:
    combineMethod = 'Asymptotic'
  elif combineOpt==2:
    combineMethod = 'Asymptotic --X-rtd TMCSO_AdaptivePseudoAsimov=50'
  elif combineOpt==3:
    combineMethod = 'HybridNew --testStat=LHC --frequentist'
  else:
    log.error('This option is not supported: %r', combineOpt)
    return __BAD__

  if scaleSingleHiggs:
    cardName = inDir+"/kt_scaled_hhbbgg_13TeV_DataCard.txt"
    resFile  = inDir+"/kt_scaled_result_"+str(combineOpt)+".log"
  else:
    cardName = inDir+"/hhbbgg_13TeV_DataCard.txt"
    resFile  = inDir+"/result_"+str(combineOpt)+".log"


  command1 = ' '.join(['combine -M', combineMethod,'-m 125 -n',Label,blinded,cardName,">",resFile,"2>&1"])
  log.info('Combine command we run:\n%s', command1)

  combExitCode = os.system(command1)

  if combineOpt in [1,2]:
    fName = 'higgsCombine'+Label+'.Asymptotic.mH125.root'
  elif combineOpt in [3]:
    fName = 'higgsCombine'+Label+'.HybridNew.mH125.root'

  outDir = inDir
  os.rename(fName, outDir+'/'+fName.replace('.root', '_%i.root'%combineOpt))

  log.info("Combine is done. %r should be produced", fName)
  return combExitCode

#########################

def runCombineOnLXBatch(inDir, doBlind, log, combineOpt=1, Label=None):
  log.info('Running combine tool.  Dir: %s Blinded: %r', inDir, doBlind)
  log.debug('inDir should be the immediate directory where the card is located')
  print "im here8"

  if doBlind and combineOpt!=3:
    # In HybridNew this option does not work
    blinded = "--run blind"
  else:
    blinded = ''

  if combineOpt==1:
    combineMethod = 'Asymptotic'
  elif combineOpt==2:
    combineMethod = 'Asymptotic --X-rtd TMCSO_AdaptivePseudoAsimov=50'
  elif combineOpt==3:
    combineMethod = 'HybridNew --testStat=LHC --frequentist '
  else:
    log.error('This option is not supported: %r', combineOpt)
    return __BAD__

  print "im here9"

  cardName = inDir+"/hhbbgg_13TeV_DataCard.txt"
  resFile  = inDir+"/result_"+str(combineOpt)
  batchFileName = inDir+"/batch_"+str(combineOpt)
  if Label!=None:
    batchFileName += "_L_"+str(Label)
    resFile += "_L_"+str(Label)
  batchFileName += ".sh"
  resFile += ".log"

  if Label == None:
    thisLabel = ''
  else:
    thisLabel = "-n " + Label

  cmssw_base =  os.getenv("CMSSW_BASE")
  outputFileStringTmp0 = '''
#!/bin/bash

cd CMSSWBASE/src/HiggsAnalysis/bbggLimits/
eval `scramv1 runtime -sh`
mkdir -p /tmp/USER/LABELP/
cd /tmp/USER/LABELP/
combine -M COMBINEMETHOD -m 125 LABEL BLINDED --datacard CARDNAME > RESFILE 2>&1

mv /tmp/USER/LABELP/OUTFILE OUTDIR

'''

  print "im here10"

  outputFileStringTmp0 = outputFileStringTmp0.replace("USER", username)
  outputFileStringTmp0 = outputFileStringTmp0.replace("LABELP", Label)
  outputFileStringTmp0 = outputFileStringTmp0.replace("BLINDED", blinded)
  outputFileStringTmp1 = outputFileStringTmp0.replace("CMSSWBASE", cmssw_base)
  outputFileStringTmp2 = outputFileStringTmp1.replace("CARDNAME", cardName)
  outputFileStringTmp4 = outputFileStringTmp2.replace("OUTDIR", inDir)

  if combineOpt < 3:
    print "im here11"
    outputFileStringTmp5 = outputFileStringTmp4.replace("LABEL", thisLabel)
    outputFileStringTmp7 = outputFileStringTmp5.replace("COMBINEMETHOD", combineMethod)
    outputFileStringTmp8 = outputFileStringTmp7.replace("RESFILE", resFile)
    outCombineFileName = 'higgsCombine'+thisLabel.replace("-n ", "") +'.Asymptotic.mH125.root'
    outputFileStringTmp9 = outputFileStringTmp8.replace("OUTFILE", outCombineFileName)
    print "im here 19"
    batchFile = open(batchFileName, "w+")
    batchFile.write(outputFileStringTmp9)
    batchFile.close()
    os.system("chmod a+rwx " + batchFileName)
#    os.system("bsub -q 1nd -o "+ resFile.replace(".log", "_batch.out") + " < " + batchFileName)
    command = "bsub -q 1nd  -J " + batchFileName.split('/')[len(batchFileName.split('/'))-1].replace('.sh', '')  + " < " + batchFileName
    print command
    os.system(command)

  if combineOpt == 3:
    print "im here12"
    #do expected bands
    quantiles = ["0.025","0.160","0.500","0.840","0.975"]
    qNames = ['m2s', 'm1s', 'central', 'p1s', 'p2s']
    for iqt,qt in enumerate(quantiles):
      myLabel = thisLabel + "_qt_"+str(qNames[iqt])
      if thisLabel == '':
        myLabel = "-n " + "_qt_"+str(qNames[iqt])
      myMethod = combineMethod + " --expectedFromGrid " + str(qt)
      myName = batchFileName.replace(".sh", "_qt_"+str(qNames[iqt])+".sh")
      myresFile = resFile.replace(".log", "_qt_"+str(qNames[iqt])+".log")
      outputFileStringTmp5 = outputFileStringTmp4.replace("LABEL", myLabel)
      outputFileStringTmp6 = outputFileStringTmp5.replace("COMBINEMETHOD", myMethod)
      outputFileStringTmp7 = outputFileStringTmp6.replace("RESFILE", myresFile)
      outCombineFileName = 'higgsCombine'+myLabel.replace("-n ", "")+".HybridNew.mH125.quant"+str(qt)+".root"
      outputFileStringTmp8 = outputFileStringTmp7.replace("OUTFILE", outCombineFileName)
      print "im here13"
      batchFile = open(myName, "w+")
      batchFile.write(outputFileStringTmp8)
      batchFile.close()
      os.system("chmod a+rwx " + myName)
      os.system("bsub -q 1nd -o "+ myresFile.replace(".log", "_batch.out") + " < " + myName)
    #if not blinded, do observed
    if not doBlind:
      print "im here14"
      myName = batchFileName.replace(".sh", "_qt_observed.sh")
      myresFile = resFile.replace(".log", "_qt_observed.log")
      myLabel = thisLabel + "_observed"
      if thisLabel == '':
        myLabel = "-n _observed"
      outputFileStringTmp5 = outputFileStringTmp4.replace("LABEL", myLabel)
      outputFileStringTmp6 = outputFileStringTmp5.replace("COMBINEMETHOD", combineMethod)
      outputFileStringTmp7 = outputFileStringTmp6.replace("RESFILE", myresFile)
      outCombineFileName = 'higgsCombine'+myLabel.replace("-n ", "")+'.HybridNew.mH125.root'
      outputFileStringTmp8 = outputFileStringTmp7.replace("OUTFILE", outCombineFileName)
      print "im here 15"
      batchFile = open(myName, "w+")
      batchFile.write(outputFileStringTmp8)
      batchFile.close()
      os.system("chmod a+rwx " + myName)
      os.system("bsub -q 1nd -o "+ myresFile.replace(".log", "_batch.out") + " < " + myName)

