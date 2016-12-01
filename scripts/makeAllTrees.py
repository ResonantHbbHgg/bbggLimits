#!/usr/bin/env python

from ROOT import *
import os,sys


import argparse
parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument('-x', nargs='+', choices=['res', 'nonres'], required=True, default=None,
                    help = "Choose which samlples to create the trees from.")
parser.add_argument('--NRW', dest="NRW", action="store_true", default=False,
                                        help="add non-resonant weights")
parser.add_argument("-v", "--verbosity",  dest="verb", action="store_true", default=False,
                                        help="Print out more stuff")

opt = parser.parse_args()


if 'nonres' in opt.x:
  #nodes = [[2, 50000]]
  
  nodes = [ ["box", 50000], ["SM", 50000], [2, 49600], [3, 50000], [4, 50000], [5, 50000], [6, 50000],
            [7, 50000], [8, 50000], [9, 49600], [10, 49800], [11, 50000], [12, 50000], [13, 50000] ]
  

  # APZ trees:
  Signals = "/afs/cern.ch/user/a/andrey/work/hh/CMSSW_8_0_8_patch1/src/APZ/fgg-ana/NotATestNov12/output_GluGluToHHTo2B2G_node_THENODE_13TeV-madgraph.root" 
  # bbggTools trees:
  #Signals = "/afs/cern.ch/user/a/andrey/work/hh/CMSSW_8_0_8_patch1/src/flashgg/bbggTools/test/RunJobs/NonResAll/output_GluGluToHHTo2B2G_node_THENODE_13TeV-madgraph_0.root" 
  #Signals = "root://eoscms//eos/cms/store/user/rateixei/HHbbgg/FlatTrees/ICHEP_Regressed4b/output_GluGluToHHTo2B2G_node_THENODE_13TeV-madgraph.root"
  #Signals = "root://eoscms//eos/cms/store/user/rateixei/HHbbgg/FlatTrees/ICHEP_Regressed4b/output_GluGluToHHTo2B2G_node_THENODE_13TeV-madgraph.root"

  Data = "root://eoscms//eos/cms/store/user/rateixei/HHbbgg/FlatTrees/ICHEP_Regressed4b/DoubleEG.root"

  dirPrefix = "LT-Nov18-APZ"

  postFix = " --MX --btagWP 0.8 "
  SFs = " --bVariation 0 --phoVariation 0"

  directory = dirPrefix
  os.system( "mkdir " + directory+"_LowMass" )
  os.system( "mkdir " + directory+"_HighMass" )


  # Sum of events of all nodes 2-13: 
  N0 = 599000
  
  for MM in nodes:
    i = MM[0]
    sigScale = 2.7/float(MM[1])
    print "DOING LowMassCat Signal, node ", i

    if opt.NRW:
      NRW = ' --NRW --NRscale '+ str(2.7/float(N0))
    else:
      NRW = ''
    command = "pyLimitTreeMaker.py -f " + Signals.replace("THENODE", str(i)) + " -o " + directory+"_LowMass" + " --min 0 --max 350 --scale " + str(sigScale) + postFix + SFs + NRW
    print command
    os.system(command)
    #	continue
    print "DOING HighMassCat Signal, node ", i
    command = "pyLimitTreeMaker.py -f " + Signals.replace("THENODE", str(i)) + " -o " + directory+"_HighMass" + " --min 350 --max 35000 --scale " + str(sigScale) + postFix + SFs + NRW
    print command
    os.system(command)

    # End of nodes
    
  # <-- indent

  if opt.NRW:
    # Merge the Nodes 2-13:
    os.system("hadd %s/LT_NR_Nodes_2to13_merged.root %s/LT_output_GluGluToHHTo2B2G_node_[1-9]*.root"%(directory+"_HighMass", directory+"_HighMass"))
    os.system("hadd %s/LT_NR_Nodes_2to13_merged.root %s/LT_output_GluGluToHHTo2B2G_node_[1-9]*.root"%(directory+"_LowMass",  directory+"_LowMass"))
  
  print "DOING LowMassCat Data"
  command = "pyLimitTreeMaker.py -f " + Data + " -o " +   directory+"_LowMass" + " --min 0 --max 350 --scale 1." + postFix
  os.system(command)
  print "DOING HighMassCat Data"
  command = "pyLimitTreeMaker.py -f " + Data + " -o " +   directory+"_HighMass" + " --min 350 --max 35000 --scale 1." + postFix
  os.system(command)

    
else:
  print 'These options are not covered yet...', opt.x
  sys.exit(1)
  
