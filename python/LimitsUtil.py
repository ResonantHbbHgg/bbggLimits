from ROOT import *
import os,sys,json,time,re
import logging
from shutil import copy
from pprint import pformat
# import pebble as pb
from multiprocessing import Pool, TimeoutError, current_process

__BAD__ = 666

#def DataCardMaker(newFolder, NCAT, sigExpStr, bkgObsStr):
def DataCardMaker(Folder, nCats, signalExp, observed, isRes = 0):
  if isRes == 0 and nCats == 1:
    print 'Resonant needs two cats!'
    sys.exit(2)

  if nCats == 2:
    inputDatacardName = os.getenv("CMSSW_BASE")+'/src/HiggsAnalysis/bbggLimits/LimitSetting/Models/LowMassResDatacardModel.txt'
    if isRes == 0:
      inputDatacardName = os.getenv("CMSSW_BASE")+'/src/HiggsAnalysis/bbggLimits/LimitSetting/Models/NonResDatacardModel.txt'

    inputDatacard = open(inputDatacardName, 'r')
    outputDatacard = open(Folder+'/datacards/hhbbgg_13TeV_DataCard.txt', 'w')
    outToWrite = ''
    for line in inputDatacard:
      outTemp = line.replace("INPUTBKGLOC", Folder+'/workspaces/hhbbgg.inputbkg_13TeV.root')
      outTemp2 = outTemp.replace("INPUTSIGLOC", Folder+'/workspaces/hhbbgg.mH125_13TeV.inputsig.root')
      outTemp3 = outTemp2.replace("OBSCAT0", '{:.0f}'.format(float(str(observed.split(',')[0]))))
      outTemp4 = outTemp3.replace("OBSCAT1", '{:.0f}'.format(float(str(observed.split(',')[1]))))
      outTemp5 = outTemp4.replace("SIGCAT0", str(signalExp.split(',')[0]))
      outTemp6 = outTemp5.replace("SIGCAT1", str(signalExp.split(',')[1]))
      if float(observed.split(',')[0]) < 11:
        newTemp1 = outTemp6
        newTemp2 = newTemp1.replace('CMS_hhbbgg_13TeV_mgg_bkg_slope3_cat0', '### CMS_hhbbgg_13TeV_mgg_bkg_slope3_cat0')
        outTemp6 = newTemp2.replace('CMS_hhbbgg_13TeV_mjj_bkg_slope3_cat0', '### CMS_hhbbgg_13TeV_mjj_bkg_slope3_cat0')
      if float(observed.split(',')[1]) < 11:
        newTemp1 = outTemp6
        newTemp2 = newTemp1.replace('CMS_hhbbgg_13TeV_mgg_bkg_slope3_cat1', '### CMS_hhbbgg_13TeV_mgg_bkg_slope3_cat1')
        outTemp6 = newTemp2.replace('CMS_hhbbgg_13TeV_mjj_bkg_slope3_cat1', '### CMS_hhbbgg_13TeV_mjj_bkg_slope3_cat1')
      outToWrite += outTemp6
    outputDatacard.write(outToWrite)
    outputDatacard.close()

  if nCats == 1 and isRes == 1:
    inputDatacardName = 'Models/HighMassResDatacardModel.txt'
    inputDatacard = open(inputDatacardName, 'r')
    outputDatacard = open(Folder+'/datacards/hhbbgg_13TeV_DataCard.txt', 'w')
    outToWrite = ''
    for line in inputDatacard:
      outTemp = line.replace("INPUTBKGLOC", Folder+'/workspaces/hhbbgg.inputbkg_13TeV.root')
      outTemp2 = outTemp.replace("INPUTSIGLOC", Folder+'/workspaces/hhbbgg.mH125_13TeV.inputsig.root')
      outTemp3 = outTemp2.replace("OBSCAT0", str(observed))
      outTemp5 = outTemp3.replace("SIGCAT0", str(signalExp))
      if float(observed.split(',')[0]) < 11:
        newTemp1 = outTemp5
        newTemp2 = newTemp1.replace('CMS_hhbbgg_13TeV_mgg_bkg_slope3_cat0', '### CMS_hhbbgg_13TeV_mgg_bkg_slope3_cat0')
        outTemp5 = newTemp2.replace('CMS_hhbbgg_13TeV_mjj_bkg_slope3_cat0', '### CMS_hhbbgg_13TeV_mjj_bkg_slope3_cat0')
      outToWrite += outTemp5
    outputDatacard.write(outToWrite)
    outputDatacard.close()

######
######
def createDir(myDir, log=None, over=True):
  if log!=None:
    log.info('Creating a new directory: %s', myDir)
  if os.path.exists(myDir):
    if log!=None:
      log.warning("\t This directory already exists: %s", myDir)
    if over:
      # Overwrite it...
      if log!=None:
        log.warning("But we will continue anyway (I will --overwrite it!)")
    else:
      if log!=None:
        log.error(' And so I exit this place...')
      print 'The directory exist and we exit. Dir = ', myDir
      sys.exit(1)
  else:
    try: os.makedirs(myDir)
    except OSError:
      if os.path.isdir(myDir): pass
      else: raise
######
######

def printTime(t1, t2, log):
  tNew = time.time()
  log.debug('Time since start of worker: %.2f sec; since previous point: %.2f sec' %(tNew-t2,tNew-t1))
  return tNew

######
######
def runCombineOnLXBatch(inDir, doBlind, log, combineOpt=1, Label=None):
  log.info('Running combine tool.  Dir: %s Blinded: %r', inDir, doBlind)
  log.debug('inDir should be the immediate directory where the card is located')

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
combine -M COMBINEMETHOD -m 125 LABEL BLINDED CARDNAME > RESFILE 

mv OUTFILE OUTDIR

'''

  outputFileStringTmp1 = outputFileStringTmp0.replace("CMSSWBASE", cmssw_base)
  outputFileStringTmp2 = outputFileStringTmp1.replace("CARDNAME", cardName)
  outputFileStringTmp4 = outputFileStringTmp2.replace("OUTDIR", inDir)

  if combineOpt < 3:
    outputFileStringTmp5 = outputFileStringTmp4.replace("LABEL", thisLabel)
    outputFileStringTmp6 = outputFileStringTmp5.replace("BLINDED", blinded)
    outputFileStringTmp7 = outputFileStringTmp6.replace("COMBINEMETHOD", combineMethod)
    outputFileStringTmp8 = outputFileStringTmp7.replace("RESFILE", resFile)
    outCombineFileName = 'higgsCombine'+thisLabel.replace("-n ", "") +'.Asymptotic.mH125.root'
    outputFileStringTmp9 = outputFileStringTmp8.replace("OUTFILE", outCombineFileName)

    batchFile = open(batchFileName, "w+")
    batchFile.write(outputFileStringTmp9)
    batchFile.close()
    os.system("chmod a+rwx " + batchFileName)
    os.system("bsub -q 1nh < " + batchFileName)

  if combineOpt == 3:
    #do expected bands
    quantiles = [0.025,0.16,0.5,0.84,0.975]
    qNames = ['m2s', 'm1s', 'central', 'p1s', 'p2s']
    for iqt,qt in enumerate(quantiles):
      thisLabel += "_qt_"+str(qNames[iqt])
      combineMethod += " --expectedFromGrid " + str(qt)
      batchFileName.replace(".sh", "_qt_"+str(qNames[iqt])+".sh")
      resFile.replace(".log", "_qt_"+str(qNames[iqt])+".log")
      outputFileStringTmp5 = outputFileStringTmp4.replace("LABEL", thisLabel)
      outputFileStringTmp6 = outputFileStringTmp5.replace("COMBINEMETHOD", combineMethod)
      outputFileStringTmp7 = outputFileStringTmp6.replace("RESFILE", resFile)
      outCombineFileName = 'higgsCombine'+thisLabel.replace("-n ", "")+'.HybridNew.mH125.root'
      outputFileStringTmp8 = outputFileStringTmp7.replace("OUTFILE", outCombineFileName)

      batchFile = open(batchFileName, "w+")
      batchFile.write(outputFileStringTmp8)
      batchFile.close()
      os.system("chmod a+rwx " + batchFileName)
      os.system("bsub -q 1nh < " + batchFileName)
    #if not blinded, do observed
    if not doBlind:
      batchFileName.replace(".sh", "_qt_observed.sh")
      resFile.replace(".log", "_qt_observed.log")
      thisLabel += "_observed"
      outputFileStringTmp5 = outputFileStringTmp4.replace("LABEL", thisLabel)
      outputFileStringTmp6 = outputFileStringTmp5.replace("COMBINEMETHOD", combineMethod)
      outputFileStringTmp7 = outputFileStringTmp6.replace("RESFILE", resFile)
      outCombineFileName = 'higgsCombine'+thisLabel.replace("-n ", "")+'.HybridNew.mH125.root'
      outputFileStringTmp8 = outputFileStringTmp7.replace("OUTFILE", outCombineFileName)

      batchFile = open(batchFileName, "w+")
      batchFile.write(outputFileStringTmp8)
      batchFile.close()
      os.system("chmod a+rwx " + batchFileName)
      os.system("bsub -q 1nh < " + batchFileName)
      
######
######

def runCombine(inDir, doBlind, log, combineOpt = 1, Label = None):
  log.info('Running combine tool.  Dir: %s Blinded: %r', inDir, doBlind)
  log.debug('inDir should be the immediate directory where the card is located')

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

  return combExitCode

######
######


def runFullChain(opt, Params, point=None, NRgridPoint=-1):
  #print 'Running: ', sys._getframe().f_code.co_name, " Node=",point
  # print sys._getframe().f_code
  PID = os.getpid()

  if opt.verb==0:
    logLvl = logging.ERROR
  elif opt.verb==1:
    logLvl = logging.INFO
  else:
    logLvl = logging.DEBUG

  LTDir_type  = os.getenv("CMSSW_BASE")+Params['LTDIR']
  signalModelCard = os.getenv("CMSSW_BASE")+Params['signal']['signalModelCard']
  lumi = Params['other']["integratedLumi"];
  energy = str(Params['other']["energy"])
  mass   = Params['other']["higgsMass"]
  addHiggs   = Params['other']["addHiggs"]
  doBlinding = Params['other']["doBlinding"]
  doBands = Params['other']["doBands"]
  NCAT    = Params['other']["ncat"]
  doBrazilianFlag = Params['other']["doBrazilianFlag"]
  Combinelxbatch = Params['other']['Combinelxbatch']
  doSingleLimit = Params['other']['doSingleLimit']
  drawSignalFit = Params['other']['drawSignalFit']
  doCombine       = Params['other']["runCombine"]
  useSigTheoryUnc = Params['other']["useSigTheoryUnc"]
  analysisType = Params['other']["analysisType"]
  HH   = Params['other']["HH"]
  base = Params['other']["base"]
  low  = Params['other']["low"]
  obs  = Params['other']["obs"]
  twotag=Params['other']["twotag"]
  dataName = Params['data']['name']
  combineOpt = Params['other']['combineOption']
  doBias = Params['other']['doBias']
  biasConfig = Params['other']['biasConfig']

  massCuts = [Params['other']["minMggMassFit"], Params['other']["maxMggMassFit"],
              Params['other']["minMjjMassFit"], Params['other']["maxMjjMassFit"],
              Params['other']["minSigFitMgg"],  Params['other']["maxSigFitMgg"],
              Params['other']["minSigFitMjj"],  Params['other']["maxSigFitMjj"],
              Params['other']["minHigMggFit"],  Params['other']["maxHigMggFit"],
              Params['other']["minHigMjjFit"],  Params['other']["maxHigMjjFit"]]

  if NCAT > 3:
    procLog.error("Error NCAT>3!")
    return __BAD__
  
  signalTypes = Params['signal']['types']

  if point!=None and NRgridPoint!=-1:
    print 'WARning: cannot have both the Node and grid Point. Chose one and try again'
    return __BAD__
  elif point!=None:
    Label = "_Node_"+str(point)
  elif NRgridPoint!=-1:
    Label = "_gridPoint_"+str(NRgridPoint)
  else:
    print 'WARning: using list of nodes from the json input file'
    return __BAD__

  sigCat = 0
  if point==None:
    sigCat = -1
  elif point == 'SM':
    sigCat = 0
  elif point == 'box':
    sigCat = 1
  elif point > 15:
    sigCat = int(point)
    isRes = 1
    Label.replace("Node", "Mass")
  else:
    sigCat = int(point)


  if opt.outDir:
    baseFolder="./"+opt.outDir+"_v"+str(Params['other']["version"])
  else:
    baseFolder="./bbggToolsResults_v"+str(Params['other']["version"])

  # Create PID file to track the job:
  pidfile = "/tmp/PIDs/PoolWorker"+Label+".pid"
  file(pidfile, 'w').write(str(PID))

  procName = current_process().name
  try:
    logging.basicConfig(level=logLvl,
                        format='%(asctime)s PID:%(process)d %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename='/tmp/logs/processLog_'+str(procName)+'.log',
                        filemode='w')
  except:
    print 'I got excepted!'
    return __BAD__

  procLog = logging.getLogger('Process.Log')

  procLog.info('\n\n New process Log started. PID = %d,  job label: %r\n', PID, Label)
  procLog.info("This log filename = "+logging.getLoggerClass().root.handlers[0].baseFilename)
  procLog.info('Node or Mass=%r  gridPoint=%r  PID=%r \n Options: %s',point, NRgridPoint, PID, pformat(opt))

  start = time.time()


  # For now the mass cuts are all the same, but can be changed in future.
  # ParamsForFits = {'SM': massCuts, 'box': massCuts}

  SignalFile = "/LT_output_GluGluToHHTo2B2G_node_"+str(point)+"_13TeV-madgraph.root"
  if isRes:
    SignalFile = "/LT_output_GluGluToTYPEToHHTo2B2G_M-"+str(point)+"_narrow_13TeV-madgraph.root"

  if NRgridPoint >= 0:
    SignalFile = "/LT_NR_Nodes_2to13_merged.root"

  procLog.debug('%s, %s', SignalFile, pformat(signalTypes))


  for t in signalTypes:
    newFolder = baseFolder+ str('/'+t+Label)
    thisSignalFile = SignalFile.replace("TYPE", t)

    procLog.info('Type = %s, %s', t, newFolder)

    createDir(newFolder,procLog)

    HLFactoryname= str(t+Label)
    hlf = RooStats.HLFactory(HLFactoryname, signalModelCard, False)
    w = hlf.GetWs()

    theFitter = bbgg2DFitter()
    theStyle = theFitter.style()
    gROOT.SetStyle('hggPaperStyle')

    theFitter.Initialize( w, sigCat, lumi, newFolder, energy, doBlinding, NCAT, addHiggs,
                          massCuts[0],massCuts[1],massCuts[2],
                          massCuts[3],massCuts[4],massCuts[5],
                          massCuts[6],massCuts[7],massCuts[8],
                          massCuts[9],massCuts[10],massCuts[11], NRgridPoint,
                          logging.getLoggerClass().root.handlers[0].baseFilename+'.bbgg2D')

    theFitter.SetVerbosityLevel(opt.verb)
    LTDir = LTDir_type.replace('TYPE', t)
    mass = 125.0

    openStatus = theFitter.AddSigData( mass, str(LTDir+thisSignalFile))
    if openStatus==-1:
      procLog.error('There is a problem with openStatus')
      return __BAD__
    procLog.info("\t SIGNAL ADDED. Node=%r, GridPoint=%r, type=%r", point,NRgridPoint,t)
    if opt.verb>0: p1 = printTime(start, start, procLog)

    createDir(newFolder+'/workspaces',procLog)
    createDir(newFolder+'/datacards',procLog)

    theFitter.SigModelFit( mass)
    procLog.info("\t SIGNAL FITTED. Node=%r, GridPoint=%r, type=%r", point,NRgridPoint,t)
    if opt.verb>0: p2 = printTime(p1,start, procLog)

    fileBaseName = "hhbbgg.mH"+str(mass)[0:3]+"_13TeV"
    theFitter.MakeSigWS( fileBaseName)
    procLog.info("\t SIGNAL'S WORKSPACE DONE. Node=%r, GridPoint=%r, type=%r", point,NRgridPoint,t)
    if opt.verb>0: p3 = printTime(p2,start,procLog)

    if drawSignalFit: 
      theFitter.MakePlots( mass)
      procLog.info("\t SIGNAL'S PLOT DONE. Node=%r, GridPoint=%r, type=%r", point,NRgridPoint,t)
      if opt.verb>0: p4 = printTime(p3,start,procLog)

    if addHiggs:
      procLog.debug('Here will add SM Higgs contributions')
      # theFitter.AddHigData( mass,direc,1)

    ddata = str(LTDir + '/LT_'+dataName+'.root')
    ddata = ddata.replace("MASS", str(point))

    theFitter.AddBkgData(ddata)
    procLog.info("\t BKG ADDED. Node=%r, GridPoint=%r, type=%r, data file=%s", point,NRgridPoint,t,ddata)
    if opt.verb>0: p4 = printTime(p3,start, procLog)

    if opt.verb>1:
      theFitter.PrintWorkspace();

    fitresults = theFitter.BkgModelFit( doBands, addHiggs)
    procLog.info("\t BKG FITTED. Node=%r, GridPoint=%r, type=%r", point,NRgridPoint,t)
    if opt.verb>0: p5 = printTime(p4,start,procLog)
    if fitresults==None:
      procLog.error("PROBLEM with fitresults !!")
      return __BAD__

    if opt.verb>1:
      fitresults.Print()

    wsFileBkgName = "hhbbgg.inputbkg_13TeV"
    theFitter.MakeBkgWS( wsFileBkgName);
    procLog.info("\t BKG'S WORKSPACE DONE. Node=%r, GridPoint=%r, type=%r", point,NRgridPoint,t)
    if opt.verb>0: p6 = printTime(p5,start,procLog)

    ##do fits for bias study, if needed
    if doBias:
      createDir(newFolder+'/bias',procLog)
      theFitter.MakeFitsForBias(str(os.getenv("CMSSW_BASE")+'/src/HiggsAnalysis/bbggLimits/'+biasConfig), str(newFolder+'/bias/biasWorkspace.root'))

    # This is making cards ala 8 TeV. We don't need this for now
    #theFitter.MakeDataCard( fileBaseName, wsFileBkgName, useSigTheoryUnc)
    #print "\t 8TeV DATACARD DONE"

    sigExp = []
    bkgObs = []
    for cc in xrange(NCAT):
      sigExp.append(-1)
      bkgObs.append(-1)

    sigExpStr = ''
    bkgObsStr = ''
    for cc in xrange(NCAT):
      sigExp[cc] = theFitter.GetSigExpectedCats(cc);
      if not doBlinding:
        bkgObs[cc] = theFitter.GetObservedCats(cc);

      sigExpStr += "%f" % sigExp[cc]
      bkgObsStr += "%f" % bkgObs[cc]
      if cc < NCAT-1:
        sigExpStr += ","
        bkgObsStr += ","

    # Make datacards:
    DataCardMaker(str(newFolder), NCAT, sigExpStr, bkgObsStr, isRes)
    procLog.info("\t DATACARD DONE. Node/Mass=%r, GridPoint=%r, type=%r", point,NRgridPoint,t)
    if opt.verb>0: p7 = printTime(p6,start,procLog)

    # Limits by type:
    if doSingleLimit or isRes:
      if doCombine:
        if Combinelxbatch:
          runCombineOnLXBatch(newFolder+"/datacards/", doBlinding, procLog, combineOpt, t+Label)
        else:
          runCombine(newFolder+"/datacards/", doBlinding, procLog, combineOpt, t+Label)

    

    # End of loop over Types
  ## <-- indent

  #Nonresonant data card massaging...
  if not isRes:
    # Here we merge datacars of all categories (in this case two)
    cardsToMerge = ''
    for t in signalTypes:
      cardsToMerge += baseFolder+'/'+t+Label+'/datacards/hhbbgg_13TeV_DataCard.txt '

    newDir = baseFolder+'/CombinedCard'+Label
    createDir(newDir,procLog)

    combCard = newDir+'/hhbbgg_13TeV_DataCard.txt'
    os.system("combineCards.py "+ cardsToMerge + " > " + combCard+' ')

    # Now we actually need to fix the combined card
    for t in signalTypes:
      strReplace = baseFolder+'/'+t+Label+'/datacards/'
      os.system("sed -i 's|"+strReplace+"./|./|g' "+combCard)

    if doCombine:
      for method in [1,2,3]:
        # If options 1,2,3 are provided - run the corresponding limits:
        # 1 - asymptotic, 2 - asymptotoc with adaptive azimov option; 3 - hybridnew
        # If combineOpt==4: run all of them at once
        if combineOpt!=4 and method!=combineOpt: continue
        try:
          combStatus = runCombine(newDir, doBlinding, procLog, method, Combinelxbatch, Label=Label)
        except:
          return __BAD__
        procLog.info("\t COMBINE with Option=%r is DONE. Node=%r, GridPoint=%r, type=%r \n \t Status = %r",
                    method, point,NRgridPoint,t, combStatus)
        if combStatus!=0:
          procLog.error('Combine failed...')
          # return __BAD__


  if opt.verb>0: p8 = printTime(p7,start,procLog)
  os.remove(pidfile)
    # procLog.handlers = []
  procLog.info('This process has ended. Label=%r', Label)
  return 42
