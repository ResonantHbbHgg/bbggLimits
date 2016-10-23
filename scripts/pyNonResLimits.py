#!/usr/bin/env python

from ROOT import *
import os,sys,json,time
gROOT.SetBatch()

__author__ = 'Andrey Pozdnyakov'
import argparse
parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument('-f', '--inputFile', dest="fname", type=str, default=None, required=True,
                    help="Json config file")
parser.add_argument('-o', '--outDir', dest="outDir", type=str, default=None,
                    help="Output directory (will be created).")
parser.add_argument('--NRW', dest="NRW", action="store_true", default=False,
                    help="Use non-resonant weights")
parser.add_argument('--overwrite', dest="overwrite", action="store_true", default=False,
                    help="Overwrite the results into the same directory")
parser.add_argument("-v", dest="verb", type=int, default=0,
                    help="Verbosity level: 0 - Minimal; 1 - Talk to me; 2 - talk more, I like to listen; 3 - not used; 4 - Debug messages only (minimal of other messages)")
parser.add_argument('-j', '--ncpu',dest="ncpu", type=int, default=2,
                    help="Number of cores to run on.")

opt = parser.parse_args()
#opt.func()

#if opt.hhelp:
#  parser.print_help()
#  sys.exit(1)

if opt.verb in [0,4]:
  gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
  RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)
  RooMsgService.instance().setSilentMode(True)

begin = time.time()
  
def createDir(myDir):
  if opt.verb>0 and opt.verb<4: print 'Creating a new directory: ', myDir
  if os.path.exists(myDir):
    if opt.verb<4:
      print "\t This directory already exists:", myDir
    if opt.overwrite:
      if opt.verb<4:
        print ' But we will continue anyway (--rewrite it!)'
    else:
      if opt.verb<4:
        print ' And so I exit this place...'
      sys.exit(1)
  else:
    try: os.makedirs(myDir)
    except OSError:
      if os.path.isdir(myDir): pass
      else: raise

def printTime(t1, t2):
  tNew = time.clock()
  print 'Time since start of worker: %.2f sec; since previous point: %.2f sec' %(tNew-t2,tNew-t1)
  return tNew

def runCombine(inDir, doBlind, Label = None, asimov=False):
  print 'Running combine tool.  Dir:', inDir, 'Blind:', doBlind
  if opt.verb==4:
    print 'inDir should be the immediate directory where the card is located'
  
  if doBlind:
    blinded = "--run blind"
  else:
    blinded = ''

  asimovOpt = ''
  if asimov:
    asimovOpt = '--X-rtd TMCSO_AdaptivePseudoAsimov=50'
  
  cardName = inDir+"/hhbbgg_13TeV_DataCard.txt"
  logFile  = inDir+"/result.log"

  command1 = "combine -M Asymptotic -m 125 -n "+Label+" "+blinded+" "+asimovOpt+" "+cardName+" > "+logFile+" 2>&1"
  os.system(command1)

  fName = 'higgsCombine'+Label+'.Asymptotic.mH125.root'
  
  outDir = inDir
  os.rename(fName, outDir+'/'+fName)


def runFullChain(Params, NRnode=None, NRgridPoint=None):
  #print 'Running: ', sys._getframe().f_code.co_name, " Node=",NRnode
  # print sys._getframe().f_code

  if opt.verb>1:
    print 'Node=',NRnode, 'gridPoint=',NRgridPoint,'Options:\n ', opt, 
  
  if NRnode!=None and NRgridPoint!=None:
    print 'WARning: cannot have both the Node and grid Point. Chose one and try again'
    sys.exit(1)
  elif NRnode!=None:
    Label = "_Node_"+NRnode
  elif NRgridPoint!=None:
    Label = "_gridPoint_"+NRgridPoint
  else:
    print 'WARning: must provide one of these: NRnode or NRgridPoint'
    sys.exit(1)

  start = time.clock()

  signalDir   = os.getenv("CMSSW_BASE")+Params['signal']['dir']
  signalTypes = Params['signal']['types']
  signalModelCard = os.getenv("CMSSW_BASE")+Params['signal']['signalModelCard']

  lumi = Params['other']["integratedLumi"];
  energy = str(Params['other']["energy"])
  mass   = Params['other']["higgsMass"]
  addHiggs   = Params['other']["addHiggs"]
  doBlinding = Params['other']["doBlinding"]
  doBands = Params['other']["doBands"]
  NCAT    = Params['other']["ncat"]
  doBrazilianFlag = Params['other']["doBrazilianFlag"]

  if NCAT > 3:
    print "Error NCAT>3!"
    sys.exit(1)

  doCombine       = Params['other']["runCombine"]
  useSigTheoryUnc = Params['other']["useSigTheoryUnc"]
  analysisType = Params['other']["analysisType"]
  HH   = Params['other']["HH"]
  base = Params['other']["base"]
  low  = Params['other']["low"]
  obs  = Params['other']["obs"]
  twotag=Params['other']["twotag"]

  dataDir  = os.getenv("CMSSW_BASE")+Params['data']['dir']
  dataName = Params['data']['name']

  sigCat = 1
  sigMas = '125'

  if opt.outDir:
    baseFolder=outDir+"_v"+str(Params['other']["version"])
  else:
    baseFolder="./bbggToolsResults_v"+str(Params['other']["version"])

  massCuts = [Params['other']["minMggMassFit"], Params['other']["maxMggMassFit"],
              Params['other']["minMjjMassFit"], Params['other']["maxMjjMassFit"],
              Params['other']["minSigFitMgg"],  Params['other']["maxSigFitMgg"],
              Params['other']["minSigFitMjj"],  Params['other']["maxSigFitMjj"],
              Params['other']["minHigMggFit"],  Params['other']["maxHigMggFit"],
              Params['other']["minHigMjjFit"],  Params['other']["maxHigMjjFit"]]

  # For now the mass cuts are all the same, but can be changed in future.
  # ParamsForFits = {'SM': massCuts, 'box': massCuts}

  NonResSignalFile = "/LT_output_GluGluToHHTo2B2G_node_"+NRnode+"_13TeV-madgraph.root"


  for t in signalTypes:
    newFolder = baseFolder+ str('/'+t+Label)
    if opt.verb>1:
      print 'Type = ', t, newFolder

    createDir(newFolder)

    HLFactoryname= str(t+Label)
    hlf = RooStats.HLFactory(HLFactoryname, signalModelCard, False)
    w = hlf.GetWs()

    theFitter = bbgg2DFitter()
    theFitter.Initialize( w, sigCat, lumi, newFolder, energy, doBlinding, NCAT, addHiggs,
                          massCuts[0],massCuts[1],massCuts[2],
                          massCuts[3],massCuts[4],massCuts[5],
                          massCuts[6],massCuts[7],massCuts[8],
                          massCuts[9],massCuts[10],massCuts[11])

    theFitter.SetVerbosityLevel(opt.verb)
    theFitter.style()

    signalDir = signalDir.replace('TYPE', t)
    dataDir   = dataDir.replace('TYPE', t)
    mass = 125.0

    openStatus = theFitter.AddSigData( mass, str(signalDir+NonResSignalFile))
    if openStatus==-1:
      print 'There is a problem with openStatus'
      sys.exit(1)
    print "\t SIGNAL ADDED. Node=",NRnode, 'type=',t
    if opt.verb>0: p1 = printTime(start, start)
    
    createDir(newFolder+'/workspaces')
    createDir(newFolder+'/datacards')

    theFitter.SigModelFit( mass)
    print "\t SIGNAL FITTED. Node=",NRnode, 'type=',t
    if opt.verb>0: p2 = printTime(p1,start)

    fileBaseName = "hhbbgg.mH"+str(mass)[0:3]+"_13TeV"
    theFitter.MakeSigWS( fileBaseName)
    print "\t SIGNAL'S WORKSPACE DONE. Node=",NRnode, 'type=',t
    if opt.verb>0: p3 = printTime(p2,start)

    theFitter.MakePlots( mass)
    print "\t SIGNAL'S PLOT DONE.. Node=",NRnode, 'type=',t
    if opt.verb>0: p4 = printTime(p3,start)

    if addHiggs:
      print 'Here will add SM Higgs contributions'
      # theFitter.AddHigData( mass,direc,1)

    ddata = str(dataDir + '/LT_'+dataName+'.root')

    theFitter.AddBkgData(ddata)
    print "\t BKG ADDED. Node=",NRnode, 'type=',t
    if opt.verb>0: p4 = printTime(p3,start)

    if opt.verb==3:
      TheFitter.PrintWorkspace();

    fitresults = theFitter.BkgModelFit( doBands, addHiggs)
    print "\t BKG FITTED. Node=",NRnode, 'type=',t
    if opt.verb>0: p5 = printTime(p4,start)
    if fitresults==None:
      print "PROBLEM with fitresults !!"
      sys.exit(1)

    if opt.verb in [1,2,3]:
      fitresults.Print()

    wsFileBkgName = "hhbbgg.inputbkg_13TeV"
    theFitter.MakeBkgWS( wsFileBkgName);
    print "\t BKG'S WORKSPACE DONE. Node=",NRnode, 'type=',t
    if opt.verb>0: p6 = printTime(p5,start)

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

    # TODO: This script needs to be included as a py function:
    DCcommand = "python LimitSetting/scripts/DataCardMaker.py -f " + str(newFolder) + str(" -n %d " % NCAT) + "-s " + sigExpStr + " -o " + bkgObsStr;
    if opt.verb==4:
      print DCcommand
    os.system(DCcommand)
    print "\t DATACARD DONE. Node=",NRnode, 'type=',t

    if opt.verb>0: p7 = printTime(p6,start)

    # End of loop over Types
  ## <-- indent
  
  # Here we shall merge datacars of all categories (in this case two)
  cardsToMerge = ''
  for t in signalTypes:
    cardsToMerge += baseFolder+'/'+t+Label+'/datacards/hhbbgg_13TeV_DataCard.txt '

  newDir = baseFolder+'/CombinedCard'+Label
  createDir(newDir)

  combCard = newDir+'/hhbbgg_13TeV_DataCard.txt'
  os.system("combineCards.py "+ cardsToMerge + " > " + combCard+' ')

  # Now we actually need to fix the combined card
  for t in signalTypes:
    strReplace = baseFolder+'/'+t+Label+'/datacards/'
    os.system("sed -i 's|"+strReplace+"./|./|g' "+combCard)
  
  if doCombine:
    runCombine(newDir, doBlinding, Label=Label)
    print "\t COMBINE DONE. Node=",NRnode
  if opt.verb>0: p8 = printTime(p7,start)

  
if __name__ == "__main__":
  print "This is the __main__ part"

  gSystem.Load('libHiggsAnalysisbbggLimits')

  #workingPath = os.getcwd()
  # parentDir = os.path.abspath(os.path.join(workingPath, os.pardir))
  #if opt.verb: print workingPath

  with open(opt.fname, 'r') as fp:
    Params = json.load(fp)

  if opt.verb==2:
    print '\t Input JSON config file:'
    print json.dumps(Params, sort_keys=True,indent=4)


  nodes = Params['signal']['nodes']

  if opt.outDir:
    baseFolder=outDir+"_v"+str(Params['other']["version"])
  else:
    baseFolder="./bbggToolsResults_v"+str(Params['other']["version"])

  createDir(baseFolder)
    
  from multiprocessing import Pool, TimeoutError  
  pool = Pool(processes=opt.ncpu)

  if opt.NRW:
    print 'Running over 5D space points'
    for p in xrange(0,10):
      pool.apply_async(runFullChain, args = (Params, None,p,))
  else:
    print 'Running over nodes'
    for n in nodes:
      #for n in ['2','SM']:
      # print 'Node = ', n
    
      # Run on multiple cores:
      pool.apply_async(runFullChain, args = (Params, n,))

      # Use signle core:
      #runFullChain(Params, NRnode=n)
    
  pool.close()  
  pool.join()

  if opt.verb>0:
    end = time.time()
    print '\t Total Time: %.2f'%(end-begin)
