#!/usr/bin/env python

from ROOT import *
import os,sys

__author__ = 'Andrey Pozdnyakov'
import argparse
parser =  argparse.ArgumentParser(description='Limit Tree maker')
inputs = parser.add_mutually_exclusive_group(required=True)
inputs.add_argument('-f', '--inputFile', dest="fname", type=str, default=None,
                    help="Filename with Flat Tree")
inputs.add_argument('-i', '--fList', dest="fList", type=str, default=None,
                    help="Text file containing the list of input ROOT files.")

parser.add_argument('-o', '--outDir', dest="outDir", type=str, default=None,
                    required=True, help="Output directory (will be created).")
parser.add_argument('--min', dest="mtotMin", type=float, default=0,
                    help="Mtot minimum cut")
parser.add_argument('--max', dest="mtotMax", type=float, default=10000,
                    help="Mtot maximum cut")
parser.add_argument('-s', '--scale', dest="scale", type=float, default=1,
                    help="Scale Factor")
parser.add_argument('--NRscale', dest="NRscale", type=float, default=1,
                    help="Scale Factor for NonRes weights")
parser.add_argument('--photonCR', dest="photonCR", action="store_true", default=False,
                    help="Do photon control region")
parser.add_argument('--KF', dest="KF", action="store_true", default=False,
                    help="Use Mtot_KF to cut on mass window")
parser.add_argument('--MX', dest="MX", action="store_true", default=False,
                    help="Use MX to cut on mass window")
parser.add_argument('-t', '--tilt', dest="tilt", action="store_true", default=False,
                    help="Select tilted mass window")
parser.add_argument('--doNoCat', dest="doNoCat", action="store_true", default=False,
                    help="Don't cut on categories")
parser.add_argument('--btagWP', dest="btagWP", type=float, default=0.8,
                    help="Set btagging working point for categories")
parser.add_argument('--doCatMixed', dest="doCatMixed", action="store_true", default=False,
                    help="Do categories with mixed btagging. Cat0: 2>low, Cat1: 1<low+1>high")
parser.add_argument('--btagHigh', dest="btagHigh", type=float, default=0.8,
                    help="for mixed cat, highest value")
parser.add_argument('--btagLow', dest="btagLow", type=float, default=0.435,
                    help="for mixed cat, lowest value")
parser.add_argument('--singleCat', dest="singleCat", action="store_true", default=False,
                    help="only one category")
parser.add_argument('--bVariation', dest="bVariation", type=int, default=-999,
                    help="Apply b-tagging SF factors: 1 or -1")
parser.add_argument('--phoVariation', dest="phoVariation", type=int, default=-999,
                    help="Photon varioation (whatever that means)")
parser.add_argument('--cosThetaStar', dest="cosThetaStar", type=float, default=100,
                    help="Cut on cosThetaStar")

parser.add_argument('--NRW', dest="NRW", action="store_true", default=False,
                    help="add non-resonant weights")
parser.add_argument("-v", "--verbosity",  dest="verb", action="store_true", default=False,
                    help="Print out more stuff")

opt = parser.parse_args()
#opt.func()

#if opt.hhelp:
#  parser.print_help()
#  sys.exit(1)

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
 
  LTM.SetMax( options.mtotMax )
  LTM.SetMin( options.mtotMin )
  LTM.SetNormalization( options.scale, options.NRscale )
  LTM.IsPhotonCR( options.photonCR )
  LTM.IsMX( options.MX )
  LTM.IsKinFit( options.KF )
  LTM.SetTilt( options.tilt )
  LTM.DoNoCat( options.doNoCat )
  LTM.SetBTagWP( options.btagWP )
  LTM.DoCatMixed( options.doCatMixed )
  LTM.SetBTagWP_High( options.btagHigh )
  LTM.SetBTagWP_Low( options.btagLow )
  LTM.DoSingleCat( options.singleCat )
  LTM.DoBVariation( options.bVariation )
  LTM.DoPhoVariation( options.phoVariation )
  LTM.DoNRWeights( options.NRW )

  LTM.SetCosThetaStar(options.cosThetaStar, 0 )
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

