#!/usr/bin/env python

from ROOT import *
import os,sys
import numpy as np
gROOT.SetBatch()

__author__ = 'Andrey Pozdnyakov'

import argparse
parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument('-b','--blind', dest="blind", action="store_true", default=False,
                    help="Do not try to get observed limits.")
parser.add_argument("-c", dest="combineOpt", type=int, default=1,
                    help="Pick which limits to plot. 1-Asymptotic; 2 - Asymptotic with adaptive asimov; 3 - HybridNew")
parser.add_argument('-d', '--dir', dest="inDir", type=str, default=None, required=True,
                    help="Input directory")
parser.add_argument('-x', choices=['res', 'nodes', 'grid', 'lambda'], required=True, default=None,
                    help = "Choose which Limit plot to make.")
parser.add_argument("-v", dest="verb", type=int, default=0,
                    help="Verbosity level: 0 is minimal")

opt = parser.parse_args()


if opt.verb>0:
  print opt
#if opt.verb < 1:
#  gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
#  RooMsgService.instance().setSilentMode(True)


def getValuesFromFile(fname):
  f = TFile.Open(fname)
  if not f:
    print "This file does not exist.. Just skip this point!"
    return None

  if opt.verb>0: 
    f.Print()

  res = []
  tree = f.Get("limit")
  for i,l in enumerate(tree):
    # print i, l, l.limit
    if i==0: res.append(float(l.limit))
    if i==1: res.append(float(l.limit))
    if i==2: res.append(float(l.limit))
    if i==3: res.append(float(l.limit))
    if i==4: res.append(float(l.limit))
    if i==5: res.append(float(l.limit))

  f.Close()
  return res


if __name__ == "__main__":
  print "This is the __main__ part"
  
  #gROOT.ProcessLine(".L ./tdrstyle.C")
  gROOT.LoadMacro("./CMS_lumi.C")
  #setTDRStyle()
  #gROOT.ForceStyle()
  TH1.SetDefaultSumw2(kTRUE)

  xAxis = []
  xErr  = []
  obs   = []
  expMean    = []
  exp1SigHi  = []
  exp1SigLow = []
  exp2SigHi  = []
  exp2SigLow = []
  
  missedPoints = []

  if opt.combineOpt==0:
    fTail = '.Asymptotic.mH125.root'
  elif opt.combineOpt==1:
    fTail = '.Asymptotic.mH125_1.root'
  elif opt.combineOpt==2:
    fTail = '.Asymptotic.mH125_2.root'
  elif opt.combineOpt==3:
    fTail = '.HybridNew.mH125_3.root'
  else:
    print 'Sorry this option is not supported:', opt.combineOpt
    sys.exit(0)
    
  if opt.x == 'res':
    print 'This option is not implemented yet:', opt.x

  elif opt.x=='nodes':
    print 'Making limit plot for nodes'

    myNodes=['SM','box','2','3','4','5','6','7','8','9','10','11','12','13']

    for i,n in enumerate(myNodes):
      l = getValuesFromFile(opt.inDir+'/CombinedCard_Node_'+n+'/higgsCombine_Node_'+n+fTail)
      if l==None:
        missedPoints.append(n)
        continue

      exp2SigLow.append(l[0])
      exp1SigLow.append(l[1])
      expMean.append(l[2])
      exp1SigHi.append(l[3])
      exp2SigHi.append(l[4])
      if not opt.blind:
        obs.append(l[5])

      if opt.verb>0:
        print n,l

      xAxis.append(float(i))
      xErr.append(0.5)

  elif opt.x=='grid':
    print 'Making limit plot for 0-1507 grid points'
    
    for n in xrange(0,1506):
      l = getValuesFromFile(opt.inDir+'/CombinedCard_gridPoint_'+str(n)+'/higgsCombine_gridPoint_'+str(n)+fTail)
      if opt.verb>0:
        print n,l
      if l==None or len(l)<5:
        missedPoints.append(n)
        continue

      exp2SigLow.append(l[0])
      exp1SigLow.append(l[1])
      expMean.append(l[2])
      exp1SigHi.append(l[3])
      exp2SigHi.append(l[4])
      if not opt.blind:
        obs.append(l[5])


      xAxis.append(float(n))
      xErr.append(0.5)


  if len(xAxis) == 0:
    print 'There are no points to plot! Exiting...'
    sys.exit(0)


  print 'Missed points:', missedPoints
  # Create the arrays for graphs
  zeros_Array = np.zeros(len(xAxis),dtype = float)
  xAxis_Array = np.array(xAxis)
  xErr_Array  = zeros_Array

  if opt.x=='grid':
    xErr_Array = np.array(xErr)

  if not opt.blind:
    obs_Array = np.array(obs)
  else:
    obs_Array = zeros_Array
  exp_Array = np.array(expMean)
  exp2SigLowErr_Array = np.array([a-b for a,b in zip(expMean,exp2SigLow)])
  exp1SigLowErr_Array = np.array([a-b for a,b in zip(expMean,exp1SigLow)])
  exp1SigHiErr_Array  = np.array([b-a for a,b in zip(expMean,exp1SigHi)])
  exp2SigHiErr_Array  = np.array([b-a for a,b in zip(expMean,exp2SigHi)])

  mg = TMultiGraph()
  mg.SetTitle('')

  nPoints  = len(xAxis)
  expected = TGraphAsymmErrors(nPoints,xAxis_Array,exp_Array,zeros_Array,zeros_Array,zeros_Array,zeros_Array)
  oneSigma = TGraphAsymmErrors(nPoints,xAxis_Array,exp_Array,xErr_Array,xErr_Array,exp1SigLowErr_Array,exp1SigHiErr_Array)
  twoSigma = TGraphAsymmErrors(nPoints,xAxis_Array,exp_Array,xErr_Array,xErr_Array,exp2SigLowErr_Array,exp2SigHiErr_Array)
  observed = TGraphAsymmErrors(nPoints,xAxis_Array,obs_Array,zeros_Array,zeros_Array,zeros_Array,zeros_Array)


  if opt.x=='nodes':

    twoSigma.SetLineWidth(8)
    twoSigma.SetLineColor(kYellow)

    oneSigma.SetMarkerColor(kBlue+1)
    oneSigma.SetMarkerStyle(21)
    oneSigma.SetLineColor(kGreen+1)
    oneSigma.SetLineWidth(7)

    observed.SetMarkerStyle(20)

    mg.Add(twoSigma,'PZ')
    mg.Add(oneSigma, 'EPZ')
    
    #mg.Add(observed)

    mg.Draw('APZ')


  if opt.x=='grid':

    limDict = {}
    for i in xrange(0,nPoints):
      if opt.verb > 0:
        print i, expMean[i], exp1SigLow[i], exp1SigHi[i], exp2SigLow[i], exp2SigHi[i]

      limDict[str(i)] = {"expected": expMean[i], 
                         "one_sigma": [exp1SigLow[i], exp1SigHi[i]],
                         "two_sigma": [exp2SigLow[i], exp2SigHi[i]] }
      if not opt.blind:
        if opt.verb>0: 
          print obs[i]
        limDict[str(i)]['observed'] = obs[i]
    
    if opt.verb > 0:
      print limDict

    import json
    with open('limits_grid.json', 'w') as fp:
      json.dump(limDict, fp, sort_keys=True, indent=4)
  
    twoSigma.SetLineColor(kYellow)
    twoSigma.SetLineWidth(1)
    twoSigma.SetFillColor(kYellow)

    oneSigma.SetLineColor(kGreen+1)
    oneSigma.SetFillColor(kGreen+1)

    expected.SetLineColor(kBlue+1)
    expected.SetLineWidth(2)

    observed.SetMarkerStyle(21)
    observed.SetMarkerSize(1)


    mg.Add(twoSigma)
    mg.Add(oneSigma)
    mg.Add(expected,'L')

    if not opt.blind:
      mg.Add(observed,'L')

    mg.Draw('A')


  mg.SetMinimum(0)

  mg.GetXaxis().SetTitle('Node Number')
  if opt.x=='nodes':
    mg.GetXaxis().SetLimits(-1, 14)
  if opt.x=='grid':
    mg.GetXaxis().SetLimits(-10, 1520)
  mg.GetYaxis().SetTitle('#sigma(pp #rightarrow HH) #times B(HH #rightarrow bb#gamma#gamma)_{95% CL} (fb)')
  mg.SetMaximum(50)

  gPad.RedrawAxis()


  leg = TLegend(0.60,0.66,0.85,0.89)
  leg.SetTextFont(42)
  leg.SetTextSize(0.04)
  leg.SetFillStyle(0)
  leg.SetBorderSize(0)

  if opt.x=='nodes':
    leg.AddEntry(observed,"Observed", "p")
    leg.AddEntry(oneSigma,"Expected", "p")
    leg.AddEntry(oneSigma,"Expected #pm 1#sigma", "l")
    leg.AddEntry(twoSigma,"Expected #pm 2#sigma", "l")
  if opt.x=='grid':
    leg.AddEntry(observed,"Observed", "p")
    leg.AddEntry(expected,"Expected", "l")
    leg.AddEntry(oneSigma,"Expected #pm 1#sigma", "f")
    leg.AddEntry(twoSigma,"Expected #pm 2#sigma", "f")


  leg.Draw()

  CMS_lumi(c1, 4, 11, "")
  for e in ['.png']:
    c1.SaveAs(opt.inDir+'/limitPlot_'+opt.x+'_'+str(opt.combineOpt)+e)  
