#!/usr/bin/env python

from ROOT import *
import os,sys

__author__ = 'Andrey Pozdnyakov & Rafael Teixeira de Lima'
import argparse
parser =  argparse.ArgumentParser(description='Limit Tree maker')
inputs = parser.add_mutually_exclusive_group(required=True)
#i/o
inputs.add_argument('-f', '--inputFile', dest="fname", type=str, default=None,
                    help="Filename with Flat Tree")
inputs.add_argument('-i', '--fList', dest="fList", type=str, default=None,
                    help="Text file containing the list of input ROOT files.")
parser.add_argument('-o', '--outDir', dest="outDir", type=str, default=None,
                    required=True, help="Output directory (will be created).")
#mx usage
parser.add_argument('--min', dest="mtotMin", type=float, default=0,
                    help="Mtot minimum cut")
parser.add_argument('--max', dest="mtotMax", type=float, default=10000,
                    help="Mtot maximum cut")
parser.add_argument('--KF', dest="KF", action="store_true", default=False,
                    help="Use Mtot_KF to cut on mass window")
parser.add_argument('--MX', dest="MX", action="store_true", default=False,
                    help="Use MX to cut on mass window")
parser.add_argument('-t', '--tilt', dest="tilt", action="store_true", default=False,
                    help="Select tilted mass window")
#normalization
parser.add_argument('-s', '--scale', dest="scale", type=float, default=1,
                    help="Normalization (lumi*xsec/totEvs), 1 for data")
parser.add_argument('--NRscale', dest="NRscale", type=float, default=1,
                    help="Scale Factor for NonRes weights")
#categorization
parser.add_argument('--doNoCat', dest="doNoCat", action="store_true", default=False,
                    help="Don't cut on categories")
parser.add_argument('--doCatNonRes', dest="doCatNonRes", action="store_true", default=False,
                    help="Non-resonant categorization scheme")
parser.add_argument('--doCatLowMass', dest="doCatLowMass", action="store_true", default=False,
                    help="Low mass resonant categorization scheme")
parser.add_argument('--doCatHighMass', dest="doCatHighMass", action="store_true", default=False,
                    help="High mass resonant categorization scheme")
parser.add_argument('--btagTight', dest="btagTight", type=float, default=0.9,
                    help="Tight b-tagging WP")
parser.add_argument('--btagMedium', dest="btagMedium", type=float, default=0.8,
                    help="Medium b-tagging WP")
parser.add_argument('--btagLoose', dest="btagLoose", type=float, default=0.435,
                    help="Loose b-tagging WP")

#corrections
parser.add_argument('--bVariation', dest="bVariation", type=int, default=-999,
                    help="Apply b-tagging SF factors: 1 or -1")
parser.add_argument('--phoVariation', dest="phoVariation", type=int, default=-999,
                    help="Photon varioation (whatever that means)")

#nonres weights
parser.add_argument('--NRW', dest="NRW", action="store_true", default=False,
                    help="add non-resonant weights")

#other
parser.add_argument("-v", "--verbosity",  dest="verb", action="store_true", default=False,
                    help="Print out more stuff")
parser.add_argument('--photonCR', dest="photonCR", action="store_true", default=False,
                    help="Do photon control region")
parser.add_argument('--photonCRNormToSig', dest="photonCRNormToSig", action="store_true", default=False,
                    help="Do photon control region normalized to signal region")
parser.add_argument('--cosThetaStarLow', dest="cosThetaStarLow", type=float, default=100,
                    help="Lower Cut on cosThetaStar")
parser.add_argument('--cosThetaStarHigh', dest="cosThetaStarHigh", type=float, default=100,
                    help="Upper Cut on cosThetaStar")
parser.add_argument('-H', '--Help', dest="HELP", action="store_true", default=False, help="Display help message")

opt = parser.parse_args()

if opt.HELP:
  parser.print_help()
  sys.exit(1)

if opt.KF and opt.MX:
  print "Sorry, you can't use both --KF and --MX options, choose only one"
  sys.exit(1)

workingPath = os.getcwd()
parentDir = os.path.abspath(os.path.join(workingPath, os.pardir))


def createDir(myDir):
  if opt.verb: print 'Creating a new directory: ', myDir
  if not os.path.exists(myDir):
    try: os.makedirs(myDir)
    except OSError:
      if os.path.isdir(myDir): pass
      else: raise
  else:
    if opt.verb: print "\t OOps, this directory already exists:", myDir


def setAndLoop(fname, options, outFile):

  f = TFile.Open(fname)
  if not f:
    print "The file can't be open... What's up with that?"
    sys.exit(1)
    
  if f.Get("fsDir/bbggSelectionTree"):
    tr = f.Get('fsDir/bbggSelectionTree')
  elif f.Get("bbggSelectionTree"):
    tr = f.Get('bbggSelectionTree')
  else:
    print " Hey! It looks like your =bbggSelectionTree= does not exist in ", f.GetName()
    print "\n Do something about it!"
    sys.exit(1)

  LTM = bbggLTMaker(tr)

#mx usage 
  LTM.SetMax( options.mtotMax )
  LTM.SetMin( options.mtotMin )
  LTM.IsMX( options.MX )
  LTM.IsKinFit( options.KF )
  LTM.SetTilt( options.tilt )
#normalization
  LTM.SetNormalization( options.scale, options.NRscale )
#categorization
  LTM.DoNoCat( options.doNoCat )
  LTM.DoCatNonRes( options.doCatNonRes )
  LTM.DoCatHighMass( options.doCatHighMass )
  LTM.DoCatLowMass( options.doCatLowMass )
  LTM.SetBTagWP_Tight( options.btagTight )
  LTM.SetBTagWP_Medium( options.btagMedium )
  LTM.SetBTagWP_Loose( options.btagLoose )
#corrections
  LTM.DoBVariation( options.bVariation )
  LTM.DoPhoVariation( options.phoVariation )
  LTM.DoNRWeights( options.NRW )
#other
  LTM.IsPhotonCR( options.photonCR )
  LTM.IsPhotonCRNormToSig( options.photonCRNormToSig )
  LTM.SetCosThetaStarLow(options.cosThetaStarLow)
  LTM.SetCosThetaStarHigh(options.cosThetaStarHigh)

  LTM.SetOutFileName( outFile )

  LTM.Loop()

if __name__ == "__main__":
  print "This is the __main__ part"

  if opt.verb: print workingPath
  
  gSystem.Load('libHiggsAnalysisbbggLimits')

  createDir(opt.outDir)

  if opt.fname:
    # Single input file is provided
    inFiles = [opt.fname]
  elif opt.fList:
    # A list of input files is provided
    with open(opt.fList) as li:
      inFiles = li.read().splitlines() 
  else:
    print 'This should never happen.'

    
  for fname in inFiles:
    rootName = fname[fname.rfind('/')+1:]
    if opt.verb: print rootName
    outFile = os.getenv("CMSSW_BASE")+"/src/HiggsAnalysis/bbggLimits/"+opt.outDir+"/LT_"+rootName
    setAndLoop(str(fname), opt, outFile)

