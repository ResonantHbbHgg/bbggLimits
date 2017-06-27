#!/usr/bin/env python

from ROOT import *
import os,sys,json,time,re
import logging
from shutil import copy
from pprint import pformat
# import pebble as pb
from multiprocessing import Pool, TimeoutError, current_process
from HiggsAnalysis.bbggLimits.LimitsUtil import *

gROOT.SetBatch()

__author__ = 'Rafael Teixeira de Lima & Andrey Pozdnyakov'
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
parser.add_argument('--mass', dest="mass", default=None, type=str, nargs='+',
                    choices=['250','260','270','280','300','320','340','350','400','450','500','550','600','650','700', '750', '800', '900', 'all'],
                    help = "Choose the resonant mass to run")
parser.add_argument('--points', dest="points", default=None, type=parseNumList, nargs='+',
                    help = "Choose the points in the grid to run")
parser.add_argument('--overwrite', dest="overwrite", action="store_true", default=False,
                    help="Overwrite the results into the same directory")
parser.add_argument("-v", dest="verb", type=int, default=0,
                    help="Verbosity level: 0 - Minimal or no messages; 1 - INFO; 2 - DEBUG; 3 - Go crazy")
parser.add_argument('-j', '--ncpu',dest="ncpu", type=int, default=2,
                    help="Number of cores to run on.")
parser.add_argument('-t', '--timeout',dest="timeout", type=int, default=None,
                    help="Per job timeout (in seconds) for multiprocessing. Jobs will be killed if run longer than this.")
parser.add_argument('--extraLabel', dest='extraLabel', default='',
                    help='Extra label')
parser.add_argument('--analyticalRW', dest='analyticalRW', action='store_true', default=False)
parser.add_argument('--kl', dest='ARW_kl', type=float, default=1.0)
parser.add_argument('--kt', dest='ARW_kt', type=float, default=1.0)
parser.add_argument('--cg', dest='ARW_cg', type=float, default=0.0)
parser.add_argument('--c2', dest='ARW_c2', type=float, default=0.0)
parser.add_argument('--c2g', dest='ARW_c2g', type=float, default=0.0)

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


  if opt.outDir:
    baseFolder=opt.outDir+"_v"+str(Params['other']["version"])
  else:
    baseFolder="./bbggToolsResults_v"+str(Params['other']["version"])

  createDir(baseFolder, over=opt.overwrite)

  copy(opt.fname, baseFolder)

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

  res_Masses = []
  if opt.mass!=None:
    mainLog.info('Running over masses:\n'+pformat(opt.mass))
    if 'all' in opt.mass:
      myNodes=['250','260','270','280','300','320','340','350','400','450','500','550','600','650','700', '750', '800', '900']
    else:
      myNodes = opt.mass
    for n in myNodes:
      res_Masses.append((n,pool.apply_async(runFullChain, args = (opt, Params, n,-1,str(opt.extraLabel)))))

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
      res_Nodes.append((n,pool.apply_async(runFullChain, args = (opt, Params, n,-1,str(opt.extraLabel)))))

      # Use signle core:
      #runFullChain(Params, NRnode=n)

  res_Points = []
  if opt.points!=None:
    listOfPoints = list(set([item for sublist in opt.points for item in sublist]))

    mainLog.debug('Running over 5D space points:\n'+pformat(opt.points))
    for p in listOfPoints:
      res_Points.append((str(p), pool.apply_async(runFullChain, args = (opt, Params, None,p,opt.extraLabel))))

  pool.close()


  # APZ. The code below tryies to kill the processes which take too long.
  # This implementation is ugly. The better way to do this is to use pebble,
  # But it's only available in Python 3...
  # Useful posts:
  # [1] http://stackoverflow.com/questions/20055498/python-multiprocessing-pool-kill-specific-long-running-or-hung-process
  # [2] http://stackoverflow.com/questions/20991968/asynchronous-multiprocessing-with-a-worker-pool-in-python-how-to-keep-going-aft
  # [3] http://stackoverflow.com/questions/26063877/python-multiprocessing-module-join-processes-with-timeout


  # Using a modified implementation of [1]:

  print "Running over:", myNodes

  pCount=0
  totJobs = len(res_Nodes)+len(res_Points)+len(res_Masses)

  badJobs = []
  for i, r in enumerate([res_Masses, res_Nodes, res_Points]):
    badJobs.append([])

    #mainLog.debug('Type of r: %r,  length of r: %r', i, len(r))

    while r:
      print r
      sys.stdout.write("\r Progress: %.1f%%\n" % (float(pCount)*100/totJobs))
      sys.stdout.flush()
      pCount+=1

      try:
        j, res = r.pop(0)
        procCheckT = time.time()
        procRes = res.get(opt.timeout)
        mainLog.info('%d Job %s has been finished. Was waiting only for %f Seconds.' % (procRes, j, time.time()-procCheckT))

      except Exception as e:
        mainLog.warning("We have reached an exception related to %s" % (str(e)))
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
