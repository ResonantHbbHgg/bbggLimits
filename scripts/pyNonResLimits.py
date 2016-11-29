#!/usr/bin/env python

from ROOT import *
import os,sys,json,time,re
import logging
from shutil import copy
from pprint import pformat
#from sets import Set
gROOT.SetBatch()

__author__ = 'Andrey Pozdnyakov'
__BAD__ = 666

import argparse

def parseNumList(string):
  # This function is used to pass arguments like:
  # --points 2,4 5-20, 60, 200-400
  # That it, it will parse all those combinations and create lists of points to run over

  # print 'Input string:',string
  chunks = re.split('[,]', string)
  chunks = filter(None, chunks)
  mylist = []
  if len(chunks)==0:
    return None
  for ch in chunks:
    # print '\t This chunk=', ch
    try:
      mylist.append(int(ch))
      # print "It's an int. Append it"
    except ValueError:
      # print 'It is not an int. Try a pattern: number1-number2'
      m = re.match(r'\d*-\d*', ch)
      if m:
        # print 'Matched the chunk:', ch
        start,end = m.group().split('-')
        mylist.extend(range(int(start),int(end)+1))
      else:
        raise argparse.ArgumentTypeError("'" + string + "' is not in acceptable format.")
  # print mylist
  return list(set(mylist))

parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument('-f', '--inputFile', dest="fname", type=str, default=None, required=True,
                    help="Json config file")
parser.add_argument('-o', '--outDir', dest="outDir", type=str, default=None,
                    help="Output directory (will be created).")
parser.add_argument('--nodes', dest="nodes", default=None, type=str, nargs='+',
                    choices=['2','3','4','5','6','7','8','9','10','11','12','13','SM','box','all'],
                    help = "Choose the nodes to run")
parser.add_argument('--points', dest="points", default=None, type=parseNumList, nargs='+',
                    help = "Choose the points in the grid to run")
parser.add_argument('--overwrite', dest="overwrite", action="store_true", default=False,
                    help="Overwrite the results into the same directory")
parser.add_argument("-v", dest="verb", type=int, default=0,
                    help="Verbosity level: 0 - Minimal or no messages; 1 - INFO; 2 - DEBUG; 3 - Go crazy")
parser.add_argument('-j', '--ncpu',dest="ncpu", type=int, default=2,
                    help="Number of cores to run on.")
parser.add_argument('-t', '--timeout',dest="timeout", type=int, default=800,
                    help="Per job timeout (in seconds) for multiprocessing. Jobs will be killed if run longer than this.")

opt = parser.parse_args()
print opt
# opt.func()
#  parser.print_help()

if opt.verb==0:
  logLvl = logging.ERROR
elif opt.verb==1:
  logLvl = logging.INFO
else:
  logLvl = logging.DEBUG

if opt.verb < 3:
  gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
  RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)
  RooMsgService.instance().setSilentMode(True)


begin = time.time()

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

def printTime(t1, t2, log):
  tNew = time.time()
  log.debug('Time since start of worker: %.2f sec; since previous point: %.2f sec' %(tNew-t2,tNew-t1))
  return tNew

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
    combineMethod = 'HybridNew'
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


def runFullChain(Params, NRnode=None, NRgridPoint=-1):
  #print 'Running: ', sys._getframe().f_code.co_name, " Node=",NRnode
  # print sys._getframe().f_code
  PID = os.getpid()

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
  pidfile = "/tmp/PIDs/PoolWorker"+Label+".pid"
  file(pidfile, 'w').write(str(PID))

  try:
    logging.basicConfig(level=logLvl,
                        format='%(asctime)s PID:%(process)d %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename='/tmp/logs/processLog'+str(Label)+'.log',
                        filemode='w')
  except:
    print 'I got excepted!'
    return __BAD__

  procLog = logging.getLogger('Process.Log')
  procLog.info("This log filename="+logging.getLoggerClass().root.handlers[0].baseFilename)
  procLog.info('Process Log started. PID = %d', PID)
  procLog.info('Node=%r  gridPoint=%r PID=%r \n Options: %s',NRnode, NRgridPoint, PID, pformat(opt))

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
    procLog.error("Error NCAT>3!")
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

  NonResSignalFile = "/LT_output_GluGluToHHTo2B2G_node_"+str(NRnode)+"_13TeV-madgraph.root"

  if NRgridPoint >= 0:
    NonResSignalFile = "/LT_NR_Nodes_2to13_merged.root"

  procLog.debug('%s, %s', NonResSignalFile, pformat(signalTypes))


  for t in signalTypes:
    newFolder = baseFolder+ str('/'+t+Label)
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

    openStatus = theFitter.AddSigData( mass, str(LTDir+NonResSignalFile))
    if openStatus==-1:
      procLog.error('There is a problem with openStatus')
      return __BAD__
    procLog.info("\t SIGNAL ADDED. Node=%r, GridPoint=%r, type=%r", NRnode,NRgridPoint,t)
    if opt.verb>0: p1 = printTime(start, start, procLog)

    createDir(newFolder+'/workspaces',procLog)
    createDir(newFolder+'/datacards',procLog)

    theFitter.SigModelFit( mass)
    procLog.info("\t SIGNAL FITTED. Node=%r, GridPoint=%r, type=%r", NRnode,NRgridPoint,t)
    if opt.verb>0: p2 = printTime(p1,start, procLog)

    fileBaseName = "hhbbgg.mH"+str(mass)[0:3]+"_13TeV"
    theFitter.MakeSigWS( fileBaseName)
    procLog.info("\t SIGNAL'S WORKSPACE DONE. Node=%r, GridPoint=%r, type=%r", NRnode,NRgridPoint,t)
    if opt.verb>0: p3 = printTime(p2,start,procLog)

    theFitter.MakePlots( mass)
    procLog.info("\t SIGNAL'S PLOT DONE. Node=%r, GridPoint=%r, type=%r", NRnode,NRgridPoint,t)
    if opt.verb>0: p4 = printTime(p3,start,procLog)

    if addHiggs:
      procLog.debug('Here will add SM Higgs contributions')
      # theFitter.AddHigData( mass,direc,1)

    ddata = str(LTDir + '/LT_'+dataName+'.root')

    theFitter.AddBkgData(ddata)
    procLog.info("\t BKG ADDED. Node=%r, GridPoint=%r, type=%r", NRnode,NRgridPoint,t)
    if opt.verb>0: p4 = printTime(p3,start, procLog)

    if opt.verb>1:
      TheFitter.PrintWorkspace();

    fitresults = theFitter.BkgModelFit( doBands, addHiggs)
    procLog.info("\t BKG FITTED. Node=%r, GridPoint=%r, type=%r", NRnode,NRgridPoint,t)
    if opt.verb>0: p5 = printTime(p4,start,procLog)
    if fitresults==None:
      procLog.error("PROBLEM with fitresults !!")
      return __BAD__

    if opt.verb>1:
      fitresults.Print()

    wsFileBkgName = "hhbbgg.inputbkg_13TeV"
    theFitter.MakeBkgWS( wsFileBkgName);
    procLog.info("\t BKG'S WORKSPACE DONE. Node=%r, GridPoint=%r, type=%r", NRnode,NRgridPoint,t)
    if opt.verb>0: p6 = printTime(p5,start,procLog)

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
    procLog.debug(DCcommand)
    os.system(DCcommand)
    procLog.info("\t DATACARD DONE. Node=%r, GridPoint=%r, type=%r", NRnode,NRgridPoint,t)
    if opt.verb>0: p7 = printTime(p6,start,procLog)

    # End of loop over Types
  ## <-- indent

  # Here we shall merge datacars of all categories (in this case two)
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
    combineOpt = Params['other']['combineOption']
    for method in [1,2,3]:
      # If options 1,2,3 are provided - run the corresponding limits:
      # 1 - asymptotic, 2 - asymptotoc with adaptive azimov option; 3 - hybridnew
      # If combineOpt==4: run all of them at once
      if combineOpt!=4 and method!=combineOpt: continue
      try:
        combStatus = runCombine(newDir, doBlinding, procLog, method, Label=Label)
      except:
        return __BAD__
      procLog.info("\t COMBINE with Option=%r is DONE. Node=%r, GridPoint=%r, type=%r \n \t Status = %r",
                    method, NRnode,NRgridPoint,t, combStatus)
      if combStatus!=0:
        procLog.error('Combine failed...')
        # return __BAD__

  if opt.verb>0: p8 = printTime(p7,start,procLog)

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

  if opt.verb>1:
    print '\t Input JSON config file:'
    print json.dumps(Params, sort_keys=True,indent=4)


  nodes = Params['signal']['nodes']

  if opt.outDir:
    baseFolder=outDir+"_v"+str(Params['other']["version"])
  else:
    baseFolder="./bbggToolsResults_v"+str(Params['other']["version"])

  createDir(baseFolder, over=opt.overwrite)

  copy(opt.fname, baseFolder)

  # import pebble as pb
  from multiprocessing import Pool, TimeoutError

  pool = Pool(processes=opt.ncpu)

  logging.basicConfig(level=logLvl,
                      format='%(asctime)s PID:%(process)d %(name)-12s %(levelname)-8s %(message)s',
                      datefmt='%m-%d %H:%M',
                      filename=baseFolder+'/mainLog_'+time.strftime("%Y%m%d-%H%M%S")+'.log',
                      filemode='w')

  mainLog = logging.getLogger('Main.Log')
  mainLog.info('Main Log started')

  mainLog.info(pformat(opt))

  createDir('/tmp/PIDs/',mainLog,True)
  createDir('/tmp/logs/',mainLog,True)

  res_Nodes = []
  if opt.nodes!=None:

    mainLog.info('Running over nodes:\n'+pformat(opt.nodes))
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


  res_Points = []
  if opt.points!=None:
    listOfPoints = list(set([item for sublist in opt.points for item in sublist]))

    mainLog.debug('Running over 5D space points:\n'+pformat(opt.points))
    for p in listOfPoints:
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

  pCount=0
  totJobs = len(res_Nodes)+len(res_Points)

  badJobs = []
  for i, r in enumerate([res_Nodes, res_Points]):
    badJobs.append([])

    #mainLog.debug('Type of r: %r,  length of r: %r', i, len(r))

    while r:
      sys.stdout.write("\r Progress: %.1f%%\n" % (float(pCount)*100/totJobs))
      sys.stdout.flush()
      pCount+=1

      try:
        j, res = r.pop(0)
        procCheckT = time.time()
        procRes = res.get(opt.timeout)
        mainLog.info('%d Job %s has been finished. Was waiting only for %f Seconds.' % (procRes, j, time.time()-procCheckT))

      except Exception as e:
        mainLog.warning(str(e))
        mainLog.warning("%s is timed out! It's been %f sec that you're running, dear %s" % (j, time.time()-procCheckT, j))
        mainLog.warning("That is too long... Because of that we've gotta kill you. Sorry.")

        # We know which process gave us an exception: it is "j" in "i", so let's kill it!
        # First, let's get the PID of that process:
        if i==0:
          pidfile = '/tmp/PIDs/PoolWorker_Node_'+str(j)+'.pid'
        elif i==1:
          pidfile = '/tmp/PIDs/PoolWorker_gridPoint_'+str(j)+'.pid'
        PID = None
        if os.path.isfile(pidfile):
          PID = str(open(pidfile).read())

        for p in pool._pool:
          # Here we loop over all running processes and check if PID matches with the one who's overtime:
          # print p, p.pid
          if str(p.pid)==PID:
            mainLog.debug('Found it still running indeed! :: %r, %r, %r %r', p, p.pid, p.is_alive(), p.exitcode)

            # We can also double-check how long it's been running with system 'ps' command:"
            # tt = str(subprocess.check_output('ps -p "'+str(p.pid)+'" o etimes=', shell=True)).strip()
            # print 'Run time from OS (may be way off the real time..) = ', tt

            # Now, KILL the m*$@r:
            p.terminate()
            pool._pool.remove(p)
            pool._repopulate_pool()

            badJobs[i].append(j)

            mainLog.debug('Here you go, %s, pid=%d, you have been served.', p.name, p.pid)

            os.remove(pidfile)
            break

    mainLog.debug('Broke out of the while loop.. for '+pformat(r))
  mainLog.debug('Broke out of the enumerate loop...')


  pool.terminate()
  pool.join()

  mainLog.info('Bad jobs:'+ pformat(badJobs))

  end = time.time()
  mainLog.info('\t Total Time: %.2f'%(end-begin))

  mainLog.info('My Main Log Finished')
