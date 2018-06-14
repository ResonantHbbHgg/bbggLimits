#!/usr/bin/env python

from ROOT import *
import os,sys
import HiggsAnalysis.bbggLimits.MassWindows as MW
import HiggsAnalysis.bbggLimits.SMHiggsSamples as SMHiggsSamples

import argparse
parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument('-x', nargs='+', choices=['res', 'nonres'], required=True, default=None,
                    help = "Choose which samlples to create the trees from.")
parser.add_argument('--NRW', dest="NRW", action="store_true", default=False,
                                        help="add non-resonant weights")
parser.add_argument("-v", "--verbosity",  dest="verb", action="store_true", default=False,
                                        help="Print out more stuff")
parser.add_argument("-d", "--dataDir", dest="dataDir", default=None, type=str,
                                       help="Input data directory location")
parser.add_argument("-s", "--signalDir", dest="signalDir", default=None,
                                       help="Input signal directory location")
parser.add_argument("--massNR", dest="massNR", default="400",
                                       help="Non resonant M(4body) mass categorization threshold")
parser.add_argument("-l", "--lumi", dest="lumi", default=35.87,
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
parser.add_argument("--doSMHiggs", dest='doSMHiggs', action="store_true", default=False, help="Make SM single H trees")
parser.add_argument("--ctsCut", dest="ctsCut", default=-10)
parser.add_argument("--resMass", dest="resMass", default=-100, help="Do one specific mass")
parser.add_argument('--doCatMVA', dest="doCatMVA", action="store_true", default=False,
                    help="Do MVA categorization")
parser.add_argument('--doETH', dest="doETH", action="store_true", default=False,
                    help="Do ETH MVA categorization")
parser.add_argument('--MVAHMC0', dest='MVAHMC0', type=float, default=0.982, help="MVAHMC0")
parser.add_argument('--MVAHMC1', dest='MVAHMC1', type=float, default=0.875, help="MVAHMC1")
parser.add_argument('--MVALMC0', dest='MVALMC0', type=float, default=0.982, help="MVALMC0")
parser.add_argument('--MVALMC1', dest='MVALMC1', type=float, default=0.875, help="MVALMC1")
parser.add_argument('--onlySMHH', dest='onlysmhh', action='store_true', default=False)
parser.add_argument('--LMLJBTC', dest='LMLJBTC', type=float, default=-10)
parser.add_argument('--HMLJBTC', dest='HMLJBTC', type=float, default=-10)
parser.add_argument('--LMSJBTC', dest='LMSJBTC', type=float, default=-10)
parser.add_argument('--HMSJBTC', dest='HMSJBTC', type=float, default=-10)
parser.add_argument('--btagDiffVar', dest='btagdiffvar', type=str, default='central')
parser.add_argument('--genDiPhotonFilter', dest='gendiphofilter', action='store_true', default=False)

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
  

  SignalFiles = "/output_GluGluToHHTo2B2G_node_THENODE_13TeV-madgraph.root"
  if opt.signalDir is None:
    Signals = "root://eoscms//eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/Mar82018_ForPubli_RafStyle/Signal/Hadd/" + SignalFiles
  else:
    if 'ForPubli' in opt.signalDir or 'ttHBDT' in opt.signalDir:
      SignalFiles = "/output_GluGluToHHTo2B2G_node_THENODE_13TeV-madgraph_0.root" #(different name: added _0)
      
    Signals = opt.signalDir + SignalFiles

  if opt.doSMHiggs:
    nodes = SMHiggsSamples.SMHiggsNodes
    if opt.signalDir is None:
      Signals = "root://eoscms//eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/Mar82018_ForPubli_RafStyle/Background/Hadd/THENODE"
    else:
      Signals = opt.signalDir + '/THENODE'


  DataFiles = "/DoubleEG.root"
  if opt.dataDir is None:
    Data = "root://eoscms//eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/Mar82018_ForPubli_RafStyle/Data/Hadd/DoubleEG.root"
  else:
    Data = opt.dataDir + DataFiles

  dirPrefix = opt.folder

  massOpt = " --MX --massThreshold " + str(opt.massNR) + " "
  catscheme = " --doCatNonRes --btagTight 0.9535 --btagMedium 0.8484 --btagLoose 0.5426 "
  if opt.doCatMVA:
    catscheme = " --doCatMVA --MVAHMC0 " + str(opt.MVAHMC0) + " --MVAHMC1 " + str(opt.MVAHMC1) + " --MVALMC0 " + str(opt.MVALMC0)+ " --MVALMC1 " + str(opt.MVALMC1)+ " "
    catscheme += ' --LMLJBTC ' + str(opt.LMLJBTC) + ' --HMLJBTC ' + str(opt.HMLJBTC) + ' --LMSJBTC ' + str(opt.LMSJBTC) + ' --HMSJBTC ' + str(opt.HMSJBTC) + ' '
  postFix = massOpt + catscheme + " --cosThetaStarHigh " + str(opt.ctsCut) + " "

  if (opt.gendiphofilter): postFix += " --genDiPhotonFilter "
  if (opt.doETH): postFix += " --isETH "
  SFs = " --bVariation 0 --phoVariation 0 --bDiffVariation " + opt.btagdiffvar + ' '

  directory = dirPrefix
  os.system( "mkdir " + directory+"_LowMass" )
  os.system( "mkdir " + directory+"_HighMass" )


  # Sum of events of all nodes 2-13: 
  N0 = 599000
  
  for MM in nodes:
    i = MM[0]
    if opt.onlysmhh == True and 'SM' not in str(i): continue
    sigScale = float(opt.lumi)/float(MM[1])
    if opt.doSMHiggs:
      sigScale = float(opt.lumi)*float(MM[2])/float(MM[1])

      # Temporary:
      #if 'bb' not in MM: continue
      
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
  
  if opt.dataDir != "0":
    print "DOING LowMassCat Data:"
    command = "pyLimitTreeMaker.py -f " + Data + " -o " +   directory+"_LowMass" + " --min 0 --max " + opt.massNR + " --scale 1." + postFix + doPhotonControlRegion
    print command
    os.system(command)
    print "DOING HighMassCat Data"
    command = "pyLimitTreeMaker.py -f " + Data + " -o " +   directory+"_HighMass" + " --min " + opt.massNR + " --max 35000 --scale 1." + postFix + doPhotonControlRegion 
    print command
    os.system(command)

elif 'res' in opt.x:
  masses = {
  'Radion' : [[250,49800],[260,50000],[270,48400],[280,50000],[300,49200],[320,50000],[340,50000],[350,50000],
              [400,50000],[450,50000],[500,49200],[550,50000],[600,50000],[650,50000],[700,50000],[750,50000],
              [800,50000],[900,50000]],
  'BulkGraviton' : [[250,50000], [260,50000], [270,50000], [280,49600], [300,50000], [320,50000], [340,50000], [350,50000],
                    [400,50000], [450,50000], [500,50000], [550,50000], [600,50000], [650,50000], [700,49200], [750,50000],
                    [800,49800], [900,50000], [1000,50000]]
  }

  # APZ trees:
  SignalFiles = "/output_GluGluTo" + opt.resType + "ToHHTo2B2G_M-MASS_narrow_13TeV-madgraph.root"
  if 'ForPubli' in opt.signalDir:
    SignalFiles = "/output_GluGluTo" + opt.resType + "ToHHTo2B2G_M-MASS_narrow_13TeV-madgraph_0.root"
  if opt.signalDir is None:
    print "You need to specify the input directory"
    sys.exit(1)
  else:
    Signals = opt.signalDir + SignalFiles

  DataFiles = "/DoubleEG.root"
  if opt.dataDir is None:
    Data = "root://eoscms//eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/Mar82018_ForPubli_RafStyle/Data/Hadd/DoubleEG.root"
  else:
    Data = opt.dataDir + DataFiles

  dirPrefix = opt.folder
  
  catScheme = " --doCatLowMass "
  if opt.isHighMassRes:
    catScheme = " --doCatHighMass "
  if opt.doCatMVA:
    catScheme = " --doCatMVA --MVAHMC0 " + str(opt.MVAHMC0) + " --MVAHMC1 " + str(opt.MVAHMC1) + " --MVALMC0 " + str(opt.MVALMC0)+ " --MVALMC1 " + str(opt.MVALMC1)+ " "
    catScheme += ' --LMLJBTC ' + str(opt.LMLJBTC) + ' --HMLJBTC ' + str(opt.HMLJBTC) + ' --LMSJBTC ' + str(opt.LMSJBTC) + ' --HMSJBTC ' + str(opt.HMSJBTC) + ' '

  postFix = " --isRes --MX " + catScheme+ " --btagTight 0.9535 --btagMedium 0.8484 --btagLoose 0.5426 --massThreshold " + str(opt.massNR) + " "
  #  SFs = " --bVariation 0 --phoVariation 0"
  SFs = " --bVariation 0 --phoVariation 0 --bDiffVariation " + opt.btagdiffvar + ' '

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
  
