#!/usr/bin/env python

from ROOT import *
import os,sys
import HiggsAnalysis.bbggLimits.MassWindows as MW

import argparse
parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument('-x', nargs='+', choices=['res', 'nonres'], required=True, default=None,
                    help = "Choose which samlples to create the trees from.")
parser.add_argument('--NRW', dest="NRW", action="store_true", default=False,
                                        help="add non-resonant weights")
parser.add_argument("-v", "--verbosity",  dest="verb", action="store_true", default=False,
                                        help="Print out more stuff")
parser.add_argument("-d", "--dataDir", dest="dataDir", default=None,
                                       help="Input data directory location")
parser.add_argument("-s", "--signalDir", dest="signalDir", default=None,
                                       help="Input signal directory location")
parser.add_argument("--massNR", dest="massNR", default="350",
                                       help="Non resonant M(4body) mass categorization threshold")
parser.add_argument("-l", "--lumi", dest="lumi", default=36.5,
                                       help="Integrated lumi to scale signal")
parser.add_argument("-f", "--folder", dest="folder", default="LT_",
                                       help="Output folder name")
parser.add_argument("--resType", dest="resType", choices=['Radion', 'BulkGraviton'], default="Radion",
                                       help="Radion or BulkGraviton")
parser.add_argument("--highMassRes", dest="isHighMassRes", action="store_true", default=False,
                                       help="Do high mass resonant categorization scheme (res option)")
parser.add_argument("--doPhotonCR", dest="isPhotonCR", action="store_true", default=False,
                                       help="Use photon control region")
parser.add_argument("--doPhotonCRSignalNorm", dest="isPhotonCRSignalNorm", action="store_true", default=False,
                                       help="Pick events from photon control region to match event yield of signal region")
parser.add_argument("--resMass", dest="resMass", default=-100, help="Do one specific mass")


opt = parser.parse_args()

doPhotonControlRegion = ''
if opt.isPhotonCR:
  doPhotonControlRegion = ' --photonCR '
if opt.isPhotonCRSignalNorm:
  doPhotonControlRegion = ' --photonCRNormToSig '

if 'nonres' in opt.x:
  #nodes = [[2, 50000]]
  
  nodes = [ ["box", 50000], ["SM", 50000], [2, 49600], [3, 50000], [4, 50000], [5, 50000], [6, 50000],
            [7, 50000], [8, 50000], [9, 49600], [10, 49800], [11, 50000], [12, 50000], [13, 50000] ]
  

  # APZ trees:
  SignalFiles = "/output_GluGluToHHTo2B2G_node_THENODE_13TeV-madgraph.root"
  if opt.signalDir is None:
    Signals = "/afs/cern.ch/user/a/andrey/work/hh/CMSSW_8_0_8_patch1/src/APZ/fgg-ana/NotATestNov12/" + SignalFiles
  else:
    Signals = opt.signalDir + SignalFiles

  # bbggTools trees:
  #Signals = "/afs/cern.ch/user/a/andrey/work/hh/CMSSW_8_0_8_patch1/src/flashgg/bbggTools/test/RunJobs/NonResAll/output_GluGluToHHTo2B2G_node_THENODE_13TeV-madgraph_0.root" 
  #Signals = "root://eoscms//eos/cms/store/user/rateixei/HHbbgg/FlatTrees/ICHEP_Regressed4b/output_GluGluToHHTo2B2G_node_THENODE_13TeV-madgraph.root"
  #Signals = "root://eoscms//eos/cms/store/user/rateixei/HHbbgg/FlatTrees/ICHEP_Regressed4b/output_GluGluToHHTo2B2G_node_THENODE_13TeV-madgraph.root"

  DataFiles = "/DoubleEG.root"
  if opt.dataDir is None:
    Data = "root://eoscms//eos/cms/store/user/rateixei/HHbbgg/FlatTrees/ICHEP_Regressed4b/DoubleEG.root"
  else:
    Data = opt.dataDir + DataFiles

  dirPrefix = "LT-Jan25-APZ"

  postFix = " --MX --doCatNonRes --btagTight 0.935 --btagMedium 0.8 --btagLoose 0.46 "
  SFs = " --bVariation 0 --phoVariation 0"

  directory = dirPrefix
  os.system( "mkdir " + directory+"_LowMass" )
  os.system( "mkdir " + directory+"_HighMass" )


  # Sum of events of all nodes 2-13: 
  N0 = 599000
  
  for MM in nodes:
    i = MM[0]
    sigScale = float(opt.lumi)/float(MM[1])
    print "DOING LowMassCat Signal, node ", i

    if opt.NRW:
      NRW = ' --NRW --NRscale '+ str(float(opt.lumi)/float(N0))
    else:
      NRW = ''
    command = "pyLimitTreeMaker.py -f " + Signals.replace("THENODE", str(i)) + " -o " + directory+"_LowMass" + " --min 0 --max " + opt.massNR + " --scale " + str(sigScale) + postFix + SFs + NRW
    print command
    os.system(command)
    #	continue
    print "DOING HighMassCat Signal, node ", i
    command = "pyLimitTreeMaker.py -f " + Signals.replace("THENODE", str(i)) + " -o " + directory+"_HighMass" + " --min " + opt.massNR + " --max 35000 --scale " + str(sigScale) + postFix + SFs + NRW
    print command
    os.system(command)

    # End of nodes
    
  # <-- indent

  if opt.NRW:
    # Merge the Nodes 2-13:
    os.system("hadd %s/LT_NR_Nodes_2to13_merged.root %s/LT_output_GluGluToHHTo2B2G_node_[1-9]*.root"%(directory+"_HighMass", directory+"_HighMass"))
    os.system("hadd %s/LT_NR_Nodes_2to13_merged.root %s/LT_output_GluGluToHHTo2B2G_node_[1-9]*.root"%(directory+"_LowMass",  directory+"_LowMass"))
  
  print "DOING LowMassCat Data"
  command = "pyLimitTreeMaker.py -f " + Data + " -o " +   directory+"_LowMass" + " --min 0 --max " + opt.massNR + " --scale 1." + postFix + doPhotonControlRegion
  os.system(command)
  print "DOING HighMassCat Data"
  command = "pyLimitTreeMaker.py -f " + Data + " -o " +   directory+"_HighMass" + " --min " + opt.massNR + " --max 35000 --scale 1." + postFix + doPhotonControlRegion 
  os.system(command)

elif 'res' in opt.x:
  masses = {
  'Radion' : [[250,49800],[260,50000],[270,48400],[280,50000],[300,49200],[320,50000],[340,50000],[350,50000],[400,50000],[450,50000],[500,49200],[550,50000],[600,50000],[650,50000],[700,50000],[750,50000],[800,50000],[900,50000]],
  'BulkGraviton' : [[250,50000], [260,50000], [270,50000], [280,49600], [300,50000], [320,50000], [340,50000], [350,50000], [400,50000], [450,50000], [500,50000], [550,50000], [600,50000], [650,50000], [700,49200], [750,50000], [800,49800], [900,50000], [1000,50000]]
  }

  # APZ trees:
  SignalFiles = "/output_GluGluTo" + opt.resType + "ToHHTo2B2G_M-MASS_narrow_13TeV-madgraph.root"
  if opt.signalDir is None:
    print "You need to specify the input directory"
    sys.exit(1)
  else:
    Signals = opt.signalDir + SignalFiles

  DataFiles = "/DoubleEG.root"
  if opt.dataDir is None:
    Data = "root://eoscms//eos/cms/store/user/rateixei/HHbbgg/FlatTrees/ICHEP_Regressed4b/DoubleEG.root"
  else:
    Data = opt.dataDir + DataFiles

  dirPrefix = opt.folder
  
  catScheme = " --doCatLowMass "
  if opt.isHighMassRes:
    catScheme = " --doCatHighMass "

  postFix = " --MX --tilt " + catScheme+ " --btagTight 0.935 --btagMedium 0.8 --btagLoose 0.46 " 
  SFs = " --bVariation 0 --phoVariation 0"

  directory = dirPrefix + opt.resType
  os.system( "mkdir " + directory )

  for MM in masses[opt.resType]:
    i = MM[0]
    if str(i) != str(opt.resMass) and int(opt.resMass) > 0: continue
    sigScale = float(opt.lumi)/float(MM[1])
    print "DOING LowMassCat Signal, node ", i
    minMass = str(MW.bbgg_MinMass(i))
    maxMass = str(MW.bbgg_MaxMass(i))
    command = "pyLimitTreeMaker.py -f " + Signals.replace("MASS", str(i)) + " -o " + directory + " --min "+minMass+" --max " + maxMass + " --scale " + str(sigScale) + postFix + SFs
    print command
    os.system(command)
    #   continue
    print "DOING HighMassCat Signal, node ", i
    command = "pyLimitTreeMaker.py -f " + Data + " -o " + directory + " --min "+minMass+" --max " + maxMass +" --scale 1.0 " + postFix + doPhotonControlRegion + " && mv " + directory + "/LT_DoubleEG.root "+ directory + "/LT_DoubleEG_M-" + str(i) + ".root"
    print command
    os.system(command)

    
else:
  print 'These options are not covered yet...', opt.x
  sys.exit(1)
  
