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

opt = parser.parse_args()
#opt.func()

#if opt.hhelp:
#  parser.print_help()
#  sys.exit(1)

if opt.verb in [0,4]:
  gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
  RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)
  RooMsgService.instance().setSilentMode(True)

workingPath = os.getcwd()
parentDir = os.path.abspath(os.path.join(workingPath, os.pardir))

begin = time.clock()

def createDir(myDir):
  if opt.verb>0 and opt.verb<4: print 'Creating a new directory: ', myDir
  if os.path.exists(myDir):
    if opt.verb>0 and opt.verb<4:
      print "\t OOps, this directory already exists:", myDir
    if opt.overwrite:
      if opt.verb>0 and opt.verb<4:
        print ' But we will continue anyway (--rewrite it!)'
    else:
      sys.exit(1)
  else:
    try: os.makedirs(myDir)
    except OSError:
      if os.path.isdir(myDir): pass
      else: raise

def runCombine(inDir, doBlind, asimov=False):
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
  
  cardName  = inDir+"/hhbbgg_13TeV_DataCard.txt"
  logFile   = inDir+"/log.txt"
  outDir    = inDir

  command1="combine -M Asymptotic -m 125 "+blinded+" "+asimovOpt+" --out "+outDir+" "+cardName+" >>"+logFile
  os.system(command1)


def runFullChain(Params, NRnode=None, NRgridPoint=None):
  #print sys._getframe().f_code.co_name

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
  NCAT = Params['other']["ncat"]
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
    newFolder = baseFolder+ str('/'+t+"_Node_"+NRnode)
    print 'Type = ', t, newFolder

    createDir(newFolder)

    HLFactoryname= str(t+"_Node_"+NRnode)
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
    print "\t SIGNAL ADDED"
    print 'Signal of node '+n+' and type '+ t +' added!'
    if opt.verb>0:
      p1 = time.clock()
      print 'Time since beginning: %.2f sec; since previous point: %.2f sec' %(p1-start,p1-start)


    createDir(newFolder+'/workspaces')
    createDir(newFolder+'/datacards')

    theFitter.SigModelFit( mass)
    print "\t SIGNAL FITTED"
    if opt.verb>0:
      p1ex = time.clock()
      print 'Time since beginning: %.2f sec; since previous point: %.2f sec' %(p1ex-start,p1ex-p1)

    fileBaseName = "hhbbgg.mH"+str(mass)[0:3]+"_13TeV"
    theFitter.MakeSigWS( fileBaseName)
    print "\t SIGNAL'S WORKSPACE DONE"
    if opt.verb>0:
      p2 = time.clock()
      print 'Time since beginning: %.2f sec; since previous point: %.2f sec' %(p2-start,p2-p1ex)

    theFitter.MakePlots( mass)
    print "\t SIGNAL'S PLOT DONE"
    if opt.verb>0:
      p3 = time.clock()
      print 'Time since beginning: %.2f sec; since previous point: %.2f sec' %(p3-start,p3-p2)


    if addHiggs:
      print 'Here will add SM Higgs contributions'
      # theFitter.AddHigData( mass,direc,1)

    ddata = str(dataDir + '/LT_'+dataName+'.root')

    theFitter.AddBkgData(ddata)
    print "\t BKG ADDED"
    if opt.verb>0:
      p4 = time.clock()
      print 'Time since beginning: %.2f sec; since previous point: %.2f sec' %(p4-start,p4-p3)

    if opt.verb==3:
      TheFitter.PrintWorkspace();

    fitresults = theFitter.BkgModelFit( doBands, addHiggs)
    print "\t BKG FITTED"
    if opt.verb>0:
      p5 = time.clock()
      print 'Time since beginning: %.2f sec; since previous point: %.2f sec' %(p5-start,p5-p4)
    if fitresults==None:
      print "PROBLEM with fitresults !!"
      sys.exit(1)

    if opt.verb in [1,2,3]:
      fitresults.Print()

    wsFileBkgName = "hhbbgg.inputbkg_13TeV"
    theFitter.MakeBkgWS( wsFileBkgName);
    print "\t BKG'S WORKSPACE DONE"
    if opt.verb>0:
      p6 = time.clock()
      print 'Time since beginning: %.2f sec; since previous point: %.2f sec' %(p6-start,p6-p5)


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
    print "\t DATACARD DONE"

    if opt.verb>0:
      p7 = time.clock()
      print 'Time since beginning: %.2f sec; since previous point: %.2f sec' %(p7-start,p7-p6)

      
  # Here we shall merge datacars of all categories (in this case two)
  cardsToMerge = ''
  for t in signalTypes:
    cardsToMerge += baseFolder+'/'+t+"_Node_"+NRnode+'/datacards/hhbbgg_13TeV_DataCard.txt '

  newDir = baseFolder+'/CombinedCard_Node_'+NRnode
  createDir(newDir)

  combCard = newDir+'/hhbbgg_13TeV_DataCard.txt'
  os.system("combineCards.py "+ cardsToMerge + " > " + combCard+' ')

  # Now we actually need to fix the combined card
  for t in signalTypes:
    strReplace = baseFolder+'/'+t+"_Node_"+NRnode+'/datacards/'
    os.system("sed -i 's|"+strReplace+"./|./|g' "+combCard)
  
  if doCombine:
    runCombine(newDir, doBlinding)
    print "\t COMBINE DONE"
  if opt.verb>0:
    p8 = time.clock()
    print 'Time since beginning: %.2f sec; since previous point: %.2f sec' %(p8-start,p8-p7)



if __name__ == "__main__":
  print "This is the __main__ part"

  gSystem.Load('libHiggsAnalysisbbggLimits')

  if opt.verb: print workingPath

  with open(opt.fname, 'r') as fp:
    Params = json.load(fp)

  if opt.verb==2:
    print '\t Input JSON config file:'
    print json.dumps(Params, sort_keys=True,indent=4)


  nodes = Params['signal']['nodes']

  for n in nodes:
    print 'Node = ', n
    runFullChain(Params, n)


  if opt.verb>0:
    end = time.clock()
    print '\t Total Time: %.2f'%(end-begin)
