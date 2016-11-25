#!/usr/bin/env python

from ROOT import *
import os,sys,json,time,re
from shutil import copy
gROOT.SetBatch()

__author__ = 'Andrey Pozdnyakov'
__BAD__ = 666

import argparse

def parseNumList(string):
  m = re.match(r'(\d+)(?:-(\d+))?$', string)
  # ^ (or use .split('-'). anyway you like.)
  if not m:
    raise argparse.ArgumentTypeError("'" + string + "' is not a range of number. Expected forms like '0-5' or '2'.")
  start = m.group(1)
  end = m.group(2) or start
  return list(range(int(start,10), int(end,10)+1))

parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument('-f', '--inputFile', dest="fname", type=str, default=None, required=True,
                    help="Json config file")
parser.add_argument('-o', '--outDir', dest="outDir", type=str, default=None,
                    help="Output directory (will be created).")
parser.add_argument('--nodes', dest="nodes", default=None, type=str, nargs='+',
                    choices=['2','3','4','5','6','7','8','9','10','11','12','13','SM','box','all'],
                    help = "Choose the nodes to run")
parser.add_argument('--points', dest="points", default=None, type=parseNumList,
                    help = "Choose the points in the grid to run")
parser.add_argument('--overwrite', dest="overwrite", action="store_true", default=False,
                    help="Overwrite the results into the same directory")
parser.add_argument("-v", dest="verb", type=int, default=0,
                    help="Verbosity level: 0 - Minimal; 1 - Talk to me; 2 - talk more, I like to listen; 3 - not used; 4 - Debug messages only (minimal of other messages)")
parser.add_argument('-j', '--ncpu',dest="ncpu", type=int, default=2,
                    help="Number of cores to run on.")
parser.add_argument('-t', '--timeout',dest="timeout", type=int, default=800,
                    help="Per job timeout (in seconds) for multiprocessing. Jobs will be killed if run longer than this.")

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
    if opt.verb<=4:
      print "\t This directory already exists:", myDir
    if opt.overwrite:
      if opt.verb<=4:
        print "  But we will continue anyway (-- I'll overwrite it!)"
    else:
      if opt.verb<=4:
        print ' And so I exit this place...'
      sys.exit(1)
  else:
    try: os.makedirs(myDir)
    except OSError:
      if os.path.isdir(myDir): pass
      else: raise

def printTime(t1, t2):
  tNew = time.time()
  print 'Time since start of worker: %.2f sec; since previous point: %.2f sec' %(tNew-t2,tNew-t1)
  return tNew

def runCombine(inDir, doBlind, Label = None, asimov=False):
  print 'Running combine tool.  Dir:', inDir, 'Blind:', doBlind
  if opt.verb==4:
    print '[DBG] inDir should be the immediate directory where the card is located'

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

  combExitCode = os.system(command1)

  fName = 'higgsCombine'+Label+'.Asymptotic.mH125.root'

  outDir = inDir
  os.rename(fName, outDir+'/'+fName)

  return combExitCode

#def giveMeName(pName='UnKnown'):
#  # This function is used for naming the processes in the Pool later on.
#  return pName

def runFullChain(Params, NRnode=None, NRgridPoint=-1):
  #print 'Running: ', sys._getframe().f_code.co_name, " Node=",NRnode
  # print sys._getframe().f_code
  PID = os.getpid()

  if opt.verb>1:
    print 'Node=',NRnode, '  gridPoint=',NRgridPoint,' PID=',PID, '\n Options: ', opt

  if NRnode!=None and NRgridPoint!=-1:
    print 'WARning: cannot have both the Node and grid Point. Chose one and try again'
    return __BAD__
  elif NRnode!=None:
    Label = "_Node_"+str(NRnode)
  elif NRgridPoint!=-1:
    Label = "_gridPoint_"+str(NRgridPoint)
  else:
    print 'WARning: must provide one of these: NRnode or NRgridPoint'
    return __BAD__

  # Create PID file to track the job:
  pidfile = "/tmp/PoolWorker"+Label+".pid"
  if os.path.isfile(pidfile):
    print "%s already exists, exiting" % pidfile
    return __BAD__
  file(pidfile, 'w').write(str(PID))


  start = time.time()

  LTDir_type  = os.getenv("CMSSW_BASE")+Params['LTDIR']
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
    return __BAD__

  doCombine       = Params['other']["runCombine"]
  useSigTheoryUnc = Params['other']["useSigTheoryUnc"]
  analysisType = Params['other']["analysisType"]
  HH   = Params['other']["HH"]
  base = Params['other']["base"]
  low  = Params['other']["low"]
  obs  = Params['other']["obs"]
  twotag=Params['other']["twotag"]

  dataName = Params['data']['name']

  sigCat = 0
  if NRnode==None:
    sigCat = -1
  elif NRnode == 'SM':
    sigCat = 0
  elif NRnode == 'box':
    sigCat = 1
  else:
    sigCat = int(NRnode)


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

  # TODO: unify this:
  if 'APZ' in LTDir_type:
    NonResSignalFile = "/LT_output_GluGluToHHTo2B2G_node_"+str(NRnode)+"_13TeV-madgraph.root"
  else:
    NonResSignalFile = "/LT_output_GluGluToHHTo2B2G_node_"+str(NRnode)+"_13TeV-madgraph_0.root"

  if NRgridPoint >= 0:
    NonResSignalFile = "/LT_NR_Nodes_2to13_merged.root"

  if opt.verb==4:
    print '[DBG]:', NonResSignalFile, signalTypes


  for t in signalTypes:
    newFolder = baseFolder+ str('/'+t+Label)
    if opt.verb>1:
      print 'Type = ', t, newFolder

    createDir(newFolder)

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
                          massCuts[9],massCuts[10],massCuts[11], NRgridPoint)

    theFitter.SetVerbosityLevel(opt.verb)
    LTDir = LTDir_type.replace('TYPE', t)
    mass = 125.0

    openStatus = theFitter.AddSigData( mass, str(LTDir+NonResSignalFile))
    if openStatus==-1:
      print 'There is a problem with openStatus'
      return __BAD__
    print "\t SIGNAL ADDED. Node=",NRnode, '  GridPoint=',NRgridPoint, 'type=',t
    if opt.verb>0: p1 = printTime(start, start)

    createDir(newFolder+'/workspaces')
    createDir(newFolder+'/datacards')

    theFitter.SigModelFit( mass)
    print "\t SIGNAL FITTED. Node=",NRnode, '  GridPoint=',NRgridPoint, 'type=',t
    if opt.verb>0: p2 = printTime(p1,start)

    fileBaseName = "hhbbgg.mH"+str(mass)[0:3]+"_13TeV"
    theFitter.MakeSigWS( fileBaseName)
    print "\t SIGNAL'S WORKSPACE DONE. Node=",NRnode, '  GridPoint=',NRgridPoint, 'type=',t
    if opt.verb>0: p3 = printTime(p2,start)

    theFitter.MakePlots( mass)
    print "\t SIGNAL'S PLOT DONE.. Node=",NRnode, '  GridPoint=',NRgridPoint, 'type=',t
    if opt.verb>0: p4 = printTime(p3,start)

    if addHiggs:
      print 'Here will add SM Higgs contributions'
      # theFitter.AddHigData( mass,direc,1)

    ddata = str(LTDir + '/LT_'+dataName+'.root')

    theFitter.AddBkgData(ddata)
    print "\t BKG ADDED. Node=",NRnode, '  GridPoint=',NRgridPoint, 'type=',t
    if opt.verb>0: p4 = printTime(p3,start)

    if opt.verb==3:
      TheFitter.PrintWorkspace();

    fitresults = theFitter.BkgModelFit( doBands, addHiggs)
    print "\t BKG FITTED. Node=",NRnode, '  GridPoint=',NRgridPoint, 'type=',t
    if opt.verb>0: p5 = printTime(p4,start)
    if fitresults==None:
      print "PROBLEM with fitresults !!"
      return __BAD__

    if opt.verb in [1,2,3]:
      fitresults.Print()

    wsFileBkgName = "hhbbgg.inputbkg_13TeV"
    theFitter.MakeBkgWS( wsFileBkgName);
    print "\t BKG'S WORKSPACE DONE. Node=",NRnode, '  GridPoint=',NRgridPoint, 'type=',t
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
    print "\t DATACARD DONE. Node=",NRnode, '  GridPoint=',NRgridPoint, 'type=',t

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
    combStatus = runCombine(newDir, doBlinding, Label=Label, asimov=Params['other']['AdaptivePseudoAsimov'])
    print "\t COMBINE DONE. Node=",NRnode, '  GridPoint=',NRgridPoint
    if combStatus!=0:
      return __BAD__

  if opt.verb>0: p8 = printTime(p7,start)

  os.remove(pidfile)
  return 42

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

  copy(opt.fname, baseFolder)

  # Before we run let's clean-up tmp directory off the files may be left from previous failed jobs:
  import glob
  filelist = glob.glob("/tmp/PoolWorker*.pid")
  for f in filelist:
    os.remove(f)


  # import pebble as pb
  from multiprocessing import Pool, TimeoutError, active_children

  pool = Pool(processes=opt.ncpu)

  res_Nodes = []
  if opt.nodes!=None:
    print 'Running over nodes:',opt.nodes

    if 'all' in opt.nodes:
      myNodes=['SM','box','2','3','4','5','6','7','8','9','10','11','12','13']
    else:
      myNodes = opt.nodes

    for n in myNodes:
      #for n in ['2','SM']:
      # Run on multiple cores:
      res_Nodes.append((n,pool.apply_async(runFullChain, args = (Params, n,))))

      # Use signle core:
      #runFullChain(Params, NRnode=n)

  #for p in [0, 10, 200, 400, 600, 1200, 1500, 1505]:

  res_Points = []
  if opt.points!=None:
    print 'Running over 5D space points:', opt.points
    for p in opt.points:

      res_Points.append((str(p), pool.apply_async(runFullChain, args = (Params, None,p,))))


  pool.close()


  # APZ. The code below tryies to kill the processes which take too long.
  # This implementation is ugly. The better way to do this is to use pebble,
  # But it's only available in Python 3...
  # Useful posts:
  # [1] http://stackoverflow.com/questions/20055498/python-multiprocessing-pool-kill-specific-long-running-or-hung-process
  # [2] http://stackoverflow.com/questions/20991968/asynchronous-multiprocessing-with-a-worker-pool-in-python-how-to-keep-going-aft
  # [3] http://stackoverflow.com/questions/26063877/python-multiprocessing-module-join-processes-with-timeout


  # Using a modified implementation of [1]:

  badJobs = []
  for i, r in enumerate([res_Nodes, res_Points]):
    badJobs.append([])

    if opt.verb==4:
      print 'Type of r:',i, 'Length of r:', len(r)

    while r:
      try:
        j, res = r.pop(0)
        procCheckT = time.time()
        print res.get(opt.timeout), ' Job %s has been finished. Was waiting only for %f Seconds.' % (j, time.time()-procCheckT)

      except Exception as e:
        print str(e)
        if opt.verb==4:
          print "%s is timed out! It's been %f sec that you're running, dear %s" % (j, time.time()-procCheckT, j)
          print "That is too long. Because of that we've gotta kill you. Sorry."

        # We know which process gave us an exception: it is "j" in "i", so let's kill it!
        # First, let's get the PID of that process:
        if i==0:
          pidfile = '/tmp/PoolWorker_Node_'+str(j)+'.pid'
        elif i==1:
          pidfile = '/tmp/PoolWorker_gridPoint_'+str(j)+'.pid'
        PID = None
        if os.path.isfile(pidfile):
          PID = str(open(pidfile).read())

        for p in pool._pool:
          # Here we loop over all running processes and check if PID matches with the one who's overtime:
          # print p, p.pid
          if str(p.pid)==PID:
            if opt.verb==4:
              print 'Found it still running indeed!', p, p.pid, p.is_alive(), p.exitcode

            # We can also double-check how long it's been running with system 'ps' command:"
            # tt = str(subprocess.check_output('ps -p "'+str(p.pid)+'" o etimes=', shell=True)).strip()
            # print 'Run time from OS (may be way off the real time..) = ', tt

            # Now, KILL the m*$@r:
            p.terminate
            pool._pool.remove(p)
            pool._repopulate_pool()

            badJobs[i].append(j)

            if opt.verb==4:
              print 'Here you go,',p.name, ', pid=', p.pid, ', you have been served.'

            os.remove(pidfile)
            break

    if opt.verb==4:
      print 'Broke out of the while loop..'

  pool.terminate()
  pool.join()



  # Just in case, let's remove all pid files from /tmp
  filelist = glob.glob("/tmp/PoolWorker*.pid")
  for f in filelist:
    os.remove(f)

  print 'Bad jobs:', badJobs
  if opt.verb>0:
    end = time.time()
    print '\t Total Time: %.2f'%(end-begin)
