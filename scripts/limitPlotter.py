#!/usr/bin/env python

from ROOT import *
import os,sys
import numpy as np
gROOT.SetBatch()

__author__ = 'Andrey Pozdnyakov'

import argparse
parser =  argparse.ArgumentParser(description='Limit plotting script')
parser.add_argument('-b','--blind', dest="blind", action="store_true", default=False,
                    help="Do not try to get observed limits.")
parser.add_argument("-c", dest="combineOpt", type=int, default=1,
                    help="Pick which limits to plot. 1-Asymptotic; 2 - Asymptotic with adaptive asimov; 3 - HybridNew")
parser.add_argument('-d', '--dir', dest="inDir", type=str, default=None, required=True,
                    help="Input directory")
parser.add_argument('-x', choices=['res','nodes','grid','lambda','yt','c2_1','c2_2','cg_1','cg_2','bench','klkt'],
                    required=True, default=None,
                    help = "Choose which Limit plot to make.")
parser.add_argument("--log", dest="log", action="store_true", default=False,
                    help="Make log scale (in y)")
parser.add_argument("--pdf", dest="pdf", action="store_true", default=False,
                    help="Make PDF plots along with PNGs.")
parser.add_argument("-v", dest="verb", type=int, default=0,
                    help="Verbosity level: 0 is minimal")

opt = parser.parse_args()

import HiggsAnalysis.bbggLimits.ParametersGrid as pg
import HiggsAnalysis.bbggLimits.TdrStyle as tdr
import HiggsAnalysis.bbggLimits.CMS_lumi as CMS_lumi
gridMap = pg.loadMapping_()

br = 0.26 / 100.

def filterMyPoints(x, fixedVals=None):
  if fixedVals==None:
    fixedVals = {'yt':1, 'c2':0, 'cg':0, 'c2g':0}


  for key, value in fixedVals.iteritems():
    if key=='cg' and value=='-c2g':
      # This is a special case. Lets take care of it:
      if x['cg'] != -x['c2g']:
        return False
    elif x[key] != value:
      return False
  return True


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
  if tree==None:
    if opt.verb>0:
      print "The limit tree in the file does not exist.. Just skip this point!"
    return None

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


def LambdaCurve(x, par):
  yt  = par[0]
  c2  = par[1]
  cg  = par[2]
  c2g = par[3]

  return br*pg.getCrossSectionForParameters(x[0], yt, c2, cg, c2g)[0]

def ytCurve(x, par):
  lamb  = par[0]
  c2  = par[1]
  cg  = par[2]
  c2g = par[3]

  return br*pg.getCrossSectionForParameters(lamb, x[0], c2, cg, c2g)[0]

def c2Curve(x, par):
  lamb  = par[0]
  yt  = par[1]
  cg  = par[2]
  c2g = par[3]

  return br*pg.getCrossSectionForParameters(lamb, yt, x[0], cg, c2g)[0]

def cgCurve(x, par):
  lamb  = par[0]
  yt  = par[1]
  c2  = par[2]

  return br*pg.getCrossSectionForParameters(lamb, yt, c2, x[0], -x[0])[0]

def c2gCurve(x, par):
  lamb  = par[0]
  yt  = par[1]
  c2  = par[2]
  cg  = par[3]

  return br*pg.getCrossSectionForParameters(lamb, yt, c2, cg, x[0])[0]

if __name__ == "__main__":
  print "This is the __main__ part"

  lambdaPoints = pg.getPoints(filterMyPoints, gridMap)
  print 'Lambda points:\n', lambdaPoints

  filt = {'lambda':1, 'c2':0, 'cg':0, 'c2g':0}
  ytPoints = pg.getPoints(filterMyPoints, gridMap, filt)
  print 'Yt points:\n', ytPoints

  filt = {'lambda':1, 'yt':1, 'cg':0, 'c2g':0}
  c2Points_1 = pg.getPoints(filterMyPoints, gridMap, filt)
  print 'C2_1 points:\n', c2Points_1

  filt = {'lambda':1, 'yt':1, 'cg':0.2, 'c2g':-0.2}
  c2Points_2 = pg.getPoints(filterMyPoints, gridMap, filt)
  print 'C2_2 points:\n', c2Points_2

  filt = {'lambda':1, 'yt':1, 'c2':0, 'cg':'-c2g'}
  cgPoints_1 = pg.getPoints(filterMyPoints, gridMap, filt)
  print 'Cg_1 points:\n', cgPoints_1

  filt = {'lambda':1, 'yt':1, 'c2':1, 'cg':'-c2g'}
  cgPoints_2 = pg.getPoints(filterMyPoints, gridMap, filt)
  print 'Cg_2 points:\n', cgPoints_2

  filt = {'c2':0, 'cg':0, 'c2g':0}
  klktPoints = pg.getPoints(filterMyPoints, gridMap, filt)
  print 'kt-kl points:\n', klktPoints
  klktScanList = ''

  #gROOT.LoadMacro("./CMS_lumi.C")
  tdr.setTDRStyle()
  tdrStyle.SetTitleSize(0.054, "Y")
  tdrStyle.SetTitleYOffset(1.1)

  TH1.SetDefaultSumw2(kTRUE)

  latex = TLatex()
  latex.SetNDC()
  latex.SetTextAngle(0)
  latex.SetTextColor(kBlack)

  xAxis = []
  xErr  = []
  obs   = []
  expMean    = []
  exp1SigHi  = []
  exp1SigLow = []
  exp2SigHi  = []
  exp2SigLow = []
  theo = []

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
      else:
        obs.append(0)
        
      if opt.verb>0:
        print n,l

      xAxis.append(float(i))
      xErr.append(0.5)


  elif opt.x=='bench':
    print 'Making limit plot for benchmarks'

    for i,n in enumerate(xrange(1507,1519)):
      print i,n
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
      else:
        obs.append(0)

      if opt.verb>0:
        print n,l

      xAxis.append(float(n-1506))
      xErr.append(0.5)

  else:
    print 'Making limit plot for 0-1518 grid points'
    
    count=0
    for n in xrange(0,1519):

      if opt.x=='lambda' and n not in lambdaPoints:
        continue
      if opt.x=='yt' and n not in ytPoints:
        continue
      if opt.x=='c2_1' and n not in c2Points_1:
        continue
      if opt.x=='c2_2' and n not in c2Points_2:
        continue
      if opt.x=='cg_1' and n not in cgPoints_1:
        continue
      if opt.x=='cg_2' and n not in cgPoints_2:
        continue
      if opt.x=='klkt' and n not in klktPoints:
        continue

      if n==324:
        print " This is SM point. It does not exist in the weights."
        print "\t We take it from the Nodes"
        l = getValuesFromFile(opt.inDir+'/CombinedCard_Node_SM/higgsCombine_Node_SM'+fTail)
        if l==None:
          missedPoints.append(n)
          continue
      elif n in [910, 985, 990]:
        continue
      else:
        l = getValuesFromFile(opt.inDir+'/CombinedCard_gridPoint_'+str(n)+'/higgsCombine_gridPoint_'+str(n)+fTail)
        if l==None or len(l)<5:
          missedPoints.append(n)
          continue

      if opt.verb>0:
        print n,l

      exp2SigLow.append(l[0])
      exp1SigLow.append(l[1])
      expMean.append(l[2])
      exp1SigHi.append(l[3])
      exp2SigHi.append(l[4])
      if not opt.blind:
        obs.append(l[5])
      else:
        obs.append(0)

      if opt.x=='grid':
        xAxis.append(float(n))
      if opt.x=='lambda':
        xAxis.append(float(pg.getParametersFromPoint(n,gridMap,True)['lambda']))
        # print n, float(pg.getParametersFromPoint(n,True)['lambda'])
      if opt.x=='yt':
        xAxis.append(float(pg.getParametersFromPoint(n,gridMap,True)['yt']))
      if opt.x in ['c2_1','c2_2']:
        xAxis.append(float(pg.getParametersFromPoint(n,gridMap,True)['c2']))
      if opt.x in ['cg_1','cg_2']:
        xAxis.append(float(pg.getParametersFromPoint(n,gridMap,True)['cg']))
      if opt.x in 'klkt':
        klktScanList += ' '.join([str(count), str(pg.getParametersFromPoint(n,gridMap,True)['lambda']),
                                  str(pg.getParametersFromPoint(n,gridMap,True)['yt']),
                                  str(l[2]), str(l[5]), str(l[3]), str(l[1]), str(l[4]), str(l[0]), '\n'])
        count+=1
        # Add symmetric values:
        klktScanList += ' '.join([str(count), str(-pg.getParametersFromPoint(n,gridMap,True)['lambda']),
                                  str(-pg.getParametersFromPoint(n,gridMap,True)['yt']),
                                  str(l[2]), str(l[5]), str(l[3]), str(l[1]), str(l[4]), str(l[0]), '\n'])
        count+=1
    
      try:
        theo.append(pg.getCrossSectionForPoint(n, gridMap)[0]*br)
      except:
        print 'Exception on getCrossSectionForPoint(). Point=', n
        theo.append(0)

      xErr.append(0.5)

  if opt.x=='klkt':
    # Adding fake values at the edges to make plotting script work:
    for kl in [-20, 20]:
      for kt in np.linspace(-2.5, 2.5, 11):
        klktScanList += ' '.join([str(count), str(kl), str(kt), '3.0 4.0 5.0 3.0 6.0 2.0 \n'])
        count+=1
    for kt in [-2.5, 2.5]:
      for kl in np.linspace(-20, 20, 11):
        klktScanList += ' '.join([str(count), str(kl), str(kt), '3.0 4.0 5.0 3.0 6.0 2.0 \n'])
        count+=1
    outList = open("KlKtList.txt", "w+")
    outList.write(klktScanList)
    outList.close()
    if opt.blind:
      os.system("python scripts/MakeKLKTplot.py --limitsFile KlKtList.txt --outFile KlKtList.root")
    else:
      os.system("python scripts/MakeKLKTplot.py --limitsFile KlKtList.txt --outFile KlKtList.root --unblind")
    
  if len(xAxis) == 0:
    print 'There are no points to plot! Exiting...'
    sys.exit(0)


  print 'Missed points:', missedPoints
  # Create the arrays for graphs
  zeros_Array = np.zeros(len(xAxis),dtype = float)
  xAxis_Array = np.array(xAxis)
  xErr_Array  = zeros_Array

  if len(theo)==len(xAxis):
    theo_Array = np.array(theo)
  else:
    theo_Array = np.zeros(len(xAxis),dtype = float)
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

  print expMean

  mg = TMultiGraph()
  mg.SetTitle('')

  nPoints  = len(xAxis)
  expected = TGraphAsymmErrors(nPoints,xAxis_Array,exp_Array,zeros_Array,zeros_Array,zeros_Array,zeros_Array)
  oneSigma = TGraphAsymmErrors(nPoints,xAxis_Array,exp_Array,xErr_Array,xErr_Array,exp1SigLowErr_Array,exp1SigHiErr_Array)
  twoSigma = TGraphAsymmErrors(nPoints,xAxis_Array,exp_Array,xErr_Array,xErr_Array,exp2SigLowErr_Array,exp2SigHiErr_Array)
  observed = TGraphAsymmErrors(nPoints,xAxis_Array,obs_Array,zeros_Array,zeros_Array,zeros_Array,zeros_Array)

  theory = TGraphAsymmErrors(nPoints,xAxis_Array,theo_Array,zeros_Array,zeros_Array,zeros_Array,zeros_Array)

  if opt.x=='grid':
    # Make a JSON file
    limDict = {}
    for i in xrange(0,nPoints):
      if opt.verb > 0:
        print i, xAxis[i], expMean[i], exp1SigLow[i], exp1SigHi[i], exp2SigLow[i], exp2SigHi[i]

      p = int(xAxis[i])
      limDict[p] = {"expected": expMean[i],
                    "one_sigma": [exp1SigLow[i], exp1SigHi[i]],
                    "two_sigma": [exp2SigLow[i], exp2SigHi[i]] }
      if not opt.blind:
        if opt.verb>0:
          print obs[i]
          limDict[p]['observed'] = obs[i]

    if opt.verb > 0:
      print limDict

    import json
    from json import encoder
    encoder.FLOAT_REPR = lambda o: format(o, '.4f')

    with open('limits_grid.json', 'w') as fp:
      json.dump(limDict, fp, sort_keys=True, indent=4)

  if opt.x=='grid':
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

    # mg.Add(theory,'L')

    if not opt.blind:
      mg.Add(observed,'L')

    mg.Draw('A')
    mg.GetXaxis().SetTitle('Node Number')

  else:
    twoSigma.SetLineWidth(8)
    twoSigma.SetLineColor(kYellow)
    twoSigma.SetMarkerStyle(1)

    oneSigma.SetMarkerColor(kBlue+1)
    oneSigma.SetMarkerStyle(21)
    oneSigma.SetLineColor(kGreen+1)
    oneSigma.SetLineWidth(7)

    observed.SetMarkerStyle(20)

    mg.Add(twoSigma,'PZ')
    mg.Add(oneSigma, 'EPZ')

    if not opt.blind:
      mg.Add(observed)

    mg.Draw('APZ')

    if opt.x == 'nodes':
      mg.GetXaxis().SetTitle('Node Number')
    if opt.x == 'bench':
      mg.GetXaxis().SetTitle('Benchmark Number')

    if opt.x in ['lambda', 'yt','c2_1','c2_2','cg_1','cg_2']:
      theory.SetMarkerStyle(22)
      theory.SetMarkerSize(1.2)
      theory.SetMarkerColor(kRed+2)
      #mg.Add(theory,'PC')

      if opt.x == 'lambda':
        mg.GetXaxis().SetTitle('#kappa_{#lambda}')
        thFunc = TF1('lambdaFunc', LambdaCurve, -16, 16, 4)
        thFunc.SetParameters(1,0,0,0)

      if opt.x == 'yt':
        mg.GetXaxis().SetTitle('#kappa_{t}')
        thFunc = TF1('ytFunc', ytCurve, 0.3, 2.6, 4)
        thFunc.SetParameters(1,0,0,0)

      if opt.x in ['c2_1', 'c2_2']:
        mg.GetXaxis().SetTitle('c_{2}')

        thFunc = TF1('c2Func', c2Curve, -2.5, 3, 4)
        if opt.x=='c2_1':
          thFunc.SetParameters(1,1,0,0)
        if opt.x=='c2_2':
          thFunc.SetParameters(1,1,0.2,-0.2)

      if opt.x in ['cg_1', 'cg_2']:
        mg.GetXaxis().SetTitle('c_{g}')

        thFunc = TF1('cgFunc', cgCurve, -1.1, 1.1, 3)
        if opt.x=='cg_1':
          thFunc.SetParameters(1,1,0)
        if opt.x=='cg_2':
          thFunc.SetParameters(1,1,1)

      thFunc.SetLineWidth(2)
      thFunc.SetLineColor(kRed+2)
      thFunc.Draw('L same')

  mg.SetMinimum(0)

  if opt.x=='nodes':
    mg.GetXaxis().SetLimits(-1, 14)
  if opt.x=='bench':
    mg.GetXaxis().SetLimits(0, 13)
  if opt.x=='grid':
    mg.GetXaxis().SetLimits(-10, 1520)
  if opt.x=='lambda':
    mg.GetXaxis().SetLimits(-16, 16)
  if opt.x=='yt':
    mg.GetXaxis().SetLimits(0, 3)
  if opt.x in ['c2_1','c2_2']:
    mg.GetXaxis().SetLimits(-4, 4)
  if opt.x in ['cg_1','cg_2']:
    mg.GetXaxis().SetLimits(-1.3, 1.3)
  mg.GetYaxis().SetTitle('#sigma(pp #rightarrow HH) #times B(HH #rightarrow bb#gamma#gamma)_{95% CL} (fb)')

  mg.SetMaximum(11)

  gPad.RedrawAxis()


  leg = TLegend(0.60,0.68,0.85,0.91)
  leg.SetTextFont(42)
  leg.SetTextSize(0.04)
  leg.SetFillStyle(0)
  leg.SetBorderSize(0)

  if opt.x != 'grid':
    leg.AddEntry(observed,"Observed", "p")
    leg.AddEntry(oneSigma,"Expected", "p")
    leg.AddEntry(oneSigma,"Expected #pm 1#sigma", "l")
    leg.AddEntry(twoSigma,"Expected #pm 2#sigma", "l")
    if opt.x in ['lambda', 'yt', 'c2_1', 'c2_2', 'cg_1', 'cg_2']:
      #leg.AddEntry(theory,"Theory", "p")
      leg.AddEntry(thFunc,"Theory prediction", "l")
      latex.SetTextFont(42)
      latex.SetTextSize(0.03)
      if opt.x=='lambda':
        latex.DrawLatex(0.2,0.8, '#kappa_{t}=#kappa_{t}^{SM}, c_{2}=c_{2}^{SM}, c_{g}=c_{g}^{SM}, c_{2g}=c_{2g}^{SM}')
      if opt.x=='yt':
        latex.DrawLatex(0.2,0.8, '#kappa_{#lambda}=#kappa_{#lambda}^{SM}, c_{2}=c_{2}^{SM}, c_{g}=c_{g}^{SM}, c_{2g}=c_{2g}^{SM}')
      if opt.x=='c2_1':
        latex.DrawLatex(0.2,0.8, '#kappa_{#lambda}=#kappa_{#lambda}^{SM}, #kappa_{t}=#kappa_{t}^{SM}, c_{g}=c_{g}^{SM}, c_{2g}=c_{2g}^{SM}')
      if opt.x=='c2_2':
        latex.DrawLatex(0.2,0.8, '#kappa_{#lambda}=#kappa_{#lambda}^{SM}, #kappa_{t}=#kappa_{t}^{SM}, c_{g}=-c_{2g}=0.2')
      if opt.x=='cg_1':
        latex.DrawLatex(0.2,0.8, '#kappa_{#lambda}=#kappa_{#lambda}^{SM}, #kappa_{t}=#kappa_{t}^{SM}, c_{2}=c_{2}^{SM}, c_{g}=-c_{2g}')
      if opt.x=='cg_2':
        latex.DrawLatex(0.2,0.8, '#kappa_{#lambda}=#kappa_{#lambda}^{SM}, #kappa_{t}=#kappa_{t}^{SM}, c_{2}=1, c_{g}=-c_{2g}')

  if opt.x=='grid':
    leg.AddEntry(observed,"Observed", "p")
    leg.AddEntry(expected,"Expected", "l")
    leg.AddEntry(oneSigma,"Expected #pm 1#sigma", "f")
    leg.AddEntry(twoSigma,"Expected #pm 2#sigma", "f")

  leg.Draw()


  if opt.log:
    mg.SetMinimum(0.03)
    mg.SetMaximum(170)
    gPad.SetLogy()


  CMS_lumi.CMS_lumi(c1, 4, 11)

  ext = ['.png']
  if opt.pdf: ext.append('.pdf')
  for e in ext:
    c1.SaveAs(opt.inDir+'/limitPlot_'+opt.x+'_'+str(opt.combineOpt)+e)
