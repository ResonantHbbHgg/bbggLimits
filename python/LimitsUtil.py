from ROOT import *
import os,sys,json,time,re
import logging
from shutil import copy,copy2
from pprint import pformat
from multiprocessing import Pool, TimeoutError, current_process
from HiggsAnalysis.bbggLimits.DataCardUtils import *
from HiggsAnalysis.bbggLimits.IOUtils import *
from HiggsAnalysis.bbggLimits.CombineUtils import *
from HiggsAnalysis.bbggLimits.AnalyticalReweighting import *
from os import listdir
from os.path import isfile, join
import getpass

__username__ = getpass.getuser()
__BAD__ = 666

def printTime(t1, t2, log):
  tNew = time.time()
  log.debug('Time since start of worker: %.2f sec; since previous point: %.2f sec' %(tNew-t2,tNew-t1))
  return tNew


def runFullChain(opt, Params, point=None, NRgridPoint=-1, extraLabel=''):
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
  if '/tmp' in Params['LTDIR'] or '/store' in Params['LTDIR'] or '/afs' in Params['LTDIR']:
    LTDir_type = Params['LTDIR']
    if '/store' in Params['LTDIR']:
      LTDir_type = '/eos/cms'+Params['LTDIR']

  signalModelCard = os.getenv("CMSSW_BASE")+Params['signal']['signalModelCard']
  lumi = 35.87 # Only used for plot produced by bbgg2Dfitter
  energy = str(Params['other']["energy"])
  mass   = Params['other']["higgsMass"]
  addHiggs   = Params['other']["addHiggs"]
  scaleSingleHiggs = Params['other']["scaleSingleHiggs"]
  doBlinding = Params['other']["doBlinding"]
  doBands = Params['other']["doBands"]
  NCAT    = Params['other']["ncat"]
  doBrazilianFlag = Params['other']["doBrazilianFlag"]
  Combinelxbatch = Params['other']['Combinelxbatch']
  doSingleLimit = Params['other']['doSingleLimit']
  drawSignalFit = Params['other']['drawSignalFit']
  doCombine       = Params['other']["runCombine"]
  useSigTheoryUnc = Params['other']["useSigTheoryUnc"]
  HH   = Params['other']["HH"]
  base = Params['other']["base"]
  low  = Params['other']["low"]
  obs  = Params['other']["obs"]
  twotag=Params['other']["twotag"]
  dataName = Params['data']['name']
  combineOpt = Params['other']['combineOption']
  doBias = Params['other']['doBias']
  biasConfig = Params['other']['biasConfig']
  doDoubleSidedCB = Params['other']['doDoubleSidedCB']

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
  elif opt.analyticalRW==True:
    Label = "_ARW_"
  elif point!=None:
    Label = "_Node_"+str(point)
  elif opt.Mjj_VBFCut!=None and opt.dEta_VBFCut!=None:
    Label = "VBFHH"
  elif NRgridPoint!=-1:
    Label = "_gridPoint_"+str(NRgridPoint)
  
  else:
    print 'WARning: using list of nodes from the json input file'
    return __BAD__

  sigCat = 0
  isRes = 0
  if point==None:
    sigCat = -1
  elif point == 'SM':
    sigCat = 0
  elif point == 'box':
    sigCat = 1
  elif int(point) > 15:
    sigCat = int(point)
    isRes = 1
    Label.replace("Node", "Mass")
  else:
    sigCat = int(point)

  Label +=  extraLabel

  print "Label=",Label


  if opt.outDir:
    baseFolder="./"+opt.outDir
  else:
    baseFolder="./bbggToolsResults"

  # Create PID file to track the job:
  pidfile = "/tmp/"+__username__+"/PIDs/PoolWorker"+Label+".pid"
  file(pidfile, 'w').write(str(PID))

  procName = current_process().name
  try:
    logging.basicConfig(level=logLvl,
                        format='%(asctime)s PID:%(process)d %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename=baseFolder+'/logs/processLog_'+str(procName)+Label+'.log',
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
  if "LT_StrikeBack" in LTDir_type or "MadMax" in LTDir_type or "ttH" in LTDir_type:
      SignalFile = "/LT_output_GluGluToHHTo2B2G_node_"+str(point)+"_13TeV-madgraph_0.root"
  if isRes:
    SignalFile = "/LT_output_GluGluToTYPEToHHTo2B2G_M-"+str(point)+"_narrow_13TeV-madgraph.root"
    if "RES_Mar21" in LTDir_type:
      SignalFile = "/LT_output_GluGluToTYPEToHHTo2B2G_M-"+str(point)+"_narrow_13TeV-madgraph_0.root"
    
  if NRgridPoint >= 0:
    SignalFile = "/LT_NR_Nodes_2to13_merged.root"

  if opt.analyticalRW == True:
    pointStr = "_".join(['kl',str(opt.ARW_kl),'kt',str(opt.ARW_kt),'cg',str(opt.ARW_cg),'c2',str(opt.ARW_c2),'c2g',str(opt.ARW_c2g)]).replace('.', 'p').replace('-', 'm')
    SignalFile="/LT_NR_Nodes_All_merged_"+pointStr+".root"

  if opt.Mjj_VBFCut!=None and opt.dEta_VBFCut!=None:
    SignalFile="/LT_output_VBFHHTo2B2G_CV_1_C2V_1_C3_1_13TeV-madgraph_0.root"


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
                          logging.getLoggerClass().root.handlers[0].baseFilename+'.bbgg2D', opt.analyticalRW)

    theFitter.SetVerbosityLevel(opt.verb)

    if opt.ttHTaggerCut!=None:
      theFitter.SetCut("ttHTagger > "+str(opt.ttHTaggerCut))
      if opt.verb>0:
        procLog.info('Apply the cut on ttHTagger: ' + str(opt.ttHTaggerCut))
    if opt.Mjj_VBFCut!=None and opt.dEta_VBFCut!=None:
      if not opt.removeHHVBF:
        print "mjj_vbf > "+str(opt.Mjj_VBFCut)+" && deta_vbf > "+str(opt.dEta_VBFCut)
        theFitter.SetCut("mjj_vbf > "+str(opt.Mjj_VBFCut)+" && deta_vbf > "+str(opt.dEta_VBFCut))
        procLog.info("mjj_vbf > "+str(opt.Mjj_VBFCut)+" && deta_vbf > "+str(opt.dEta_VBFCut))
      else:
        theFitter.SetCut("mjj_vbf < "+str(opt.Mjj_VBFCut)+" || deta_vbf < "+str(opt.dEta_VBFCut))
        procLog.info("mjj_vbf < "+str(opt.Mjj_VBFCut)+" || deta_vbf < "+str(opt.dEta_VBFCut))


        
    if 'HighMass' in t:
      theFitter.SetNCat0(2)
    else:
      theFitter.SetNCat0(0)

    if opt.verb>0:
      procLog.info('Using Double Sided Crystal Ball as Signal Model: %r', doDoubleSidedCB)
    if doDoubleSidedCB: theFitter.UseDoubleSidedCB()

    LTDir = LTDir_type.replace('TYPE', t)

    mass = 125.0
    if opt.verb>0:
      procLog.info('Signal File:\n'+LTDir+thisSignalFile)

    if not os.path.isfile(LTDir+thisSignalFile):
      print 'File does not exist: ', LTDir+thisSignalFile
      return __BAD__

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
      higTypes = Params['higgs']['type']
      if opt.verb>1:
        procLog.debug('Here will add SM Higgs contributions \n Higgs types: '+ pformat(higTypes))
      higgsExp = {}
      for iht,HT in enumerate(higTypes):
        higgsExp[HT] = [0,0]
        ht = higTypes[HT]
        if opt.verb>1:
          procLog.debug('iht = %r, ht = %r, HT = %r' % (iht,ht,HT))
        higFileName = str(LTDir)+"/LT_output_"+str(ht)+".root"

        exphig = theFitter.AddHigData( mass,higFileName,iht, str(HT))
        theFitter.HigModelFit(mass,iht, str(HT) )
        theFitter.MakeHigWS(str('hhbbgg.')+str(HT), iht, str(HT))

        higgsExp[HT] = [exphig[0], exphig[1]]

      if opt.verb>1:
        procLog.debug("Done SM Higgs bzz")

    ddata = str(LTDir + '/LT_'+dataName+'.root')
    ddata = ddata.replace("%MASS%", str(point))

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
      procLog.debug("\n Fit Results: \n\n"+pformat(fitresults.Print()))

    wsFileBkgName = "hhbbgg.inputbkg_13TeV"
    theFitter.MakeBkgWS( wsFileBkgName);
    procLog.info("\t BKG'S WORKSPACE DONE. Node=%r, GridPoint=%r, type=%r", point,NRgridPoint,t)
    if opt.verb>0: p6 = printTime(p5,start,procLog)

    ##do fits for bias study, if needed

    procLog.info("\t Making Fits and WS for Bias Study? * %r *  Node=%r, GridPoint=%r, type=%r", doBias, point,NRgridPoint,t)
    if doBias:
      createDir(newFolder+'/bias',procLog)
      theFitter.MakeFitsForBias(str(os.getenv("CMSSW_BASE")+'/src/HiggsAnalysis/bbggLimits/'+biasConfig), str(newFolder+'/bias/biasWorkspace.root'))

    procLog.info("\t BIAS FITS DONE. Node=%r, GridPoint=%r, type=%r", point,NRgridPoint,t)
    if opt.verb>0: p7 = printTime(p6,start,procLog)

    # print PID, "IM HERE"

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

    # print PID, "IM HERE2"

    # Make datacards:
    myLoc = os.getenv("CMSSW_BASE") + '/src/HiggsAnalysis/bbggLimits/'+newFolder
    if isRes==1:
      DataCardMaker(str(myLoc), NCAT, sigExpStr, bkgObsStr, isRes)
    elif addHiggs == 0:
      DataCardMaker(str(myLoc), NCAT, sigExpStr, bkgObsStr, isRes, t)
    else:
      DataCardMaker_wHiggs(str(myLoc), NCAT, sigExpStr, bkgObsStr, higgsExp, t)

    procLog.info("\t DATACARD DONE. Node/Mass=%r, GridPoint=%r, type=%r", point,NRgridPoint,t)
    if opt.verb>0: p8 = printTime(p7,start,procLog)

    # print PID, "IM HERE3"
    # Limits by type:
    if doSingleLimit or isRes:
      # print PID, "IM HERE4"
      if doCombine:
        # print PID, "IM HERE5"
        if Combinelxbatch:
          # print PID, "IM HERE6"
          runCombineOnLXBatch(myLoc+"/datacards/", doBlinding, procLog, combineOpt, t+Label)
        else:
          # print PID, "IM HERE7"
          runCombine(newFolder+"/datacards/", doBlinding, procLog, combineOpt, Combinelxbatch, t+Label)



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
    os.system("combineCards.py "+ cardsToMerge + " > " + combCard)

    # Now we actually need to fix the combined card
    for t in signalTypes:
      strReplace = baseFolder+'/'+t+Label+'/datacards/'+os.getenv("CMSSW_BASE")+'/src/HiggsAnalysis/bbggLimits/'
      os.system("sed -i 's|"+strReplace+"||g' "+combCard)
      procLog.info("String to replace: "+ strReplace)

    if opt.analyticalRW:
      # Make another datacard, with entries for kt-reweihgting of single Higgs background
      ktScaled = newDir+'/kt_scaled_hhbbgg_13TeV_DataCard.txt'
      copy2(combCard, ktScaled)
      with open(ktScaled, "a") as myfile:
        HigScale = opt.ARW_kt*opt.ARW_kt
        appendString = '\n \n'
        appendString+= 'h_norm_ggh rateParam * ggh %.4f \n' % HigScale
        appendString+= 'nuisance edit freeze h_norm_ggh \n'
        appendString+= 'h_norm_tth rateParam * tth %.4f \n' % HigScale
        appendString+= 'nuisance edit freeze h_norm_tth \n'
        myfile.write(appendString)
        
    if doCombine:
      if Combinelxbatch:
        # print PID, "IM HERE6"
        myLoc = os.getenv("CMSSW_BASE") + '/src/HiggsAnalysis/bbggLimits/' + newDir
        runCombineOnLXBatch(myLoc+"/", doBlinding, procLog, combineOpt, "CombinedCard"+Label)
      else:
        for method in [1,2,3]:
          # If options 1,2,3 are provided - run the corresponding limits:
          # 1 - asymptotic, 2 - asymptotoc with adaptive azimov option; 3 - hybridnew
          # If combineOpt==4: run all of them at once
          if combineOpt!=4 and method!=combineOpt: continue
          try:
            combStatus = runCombine(newDir, doBlinding, procLog, method, Combinelxbatch, Label, scaleSingleHiggs)
          except:
            return __BAD__
          procLog.info("\t COMBINE with Option=%r is DONE. Node=%r, GridPoint=%r, type=%r \n \t Status = %r",
                       method, point,NRgridPoint,t, combStatus)
          if combStatus!=0:
            procLog.error('Combine failed...')
            # return __BAD__


  if opt.verb>0: p9 = printTime(p8,start,procLog)
  os.remove(pidfile)
  # procLog.handlers = []
  procLog.info('This process has ended. Label=%r', Label)
  return 42
