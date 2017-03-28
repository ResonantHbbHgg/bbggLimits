#!/usr/bin/env python

from ROOT import *
import sys, getopt, os
import argparse
import glob
import HiggsAnalysis.bbggLimits.tdrStyle as tdr

gROOT.SetBatch(kTRUE)

parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument('-f', '--inputFolders', dest="folders", default=None, type=str, nargs='+', required=True,
                    help="input folders")
parser.add_argument('-n', '--inputNames', dest="names", default=None, type=str, nargs='+', required=True,
                    help="input folders cat names")
parser.add_argument('-l', '--lumi', dest='lumi', default='36.5', type=str, help='Integrated luminosoty')
parser.add_argument('--label', dest='label', default='', type=str, help='Label')
parser.add_argument('--log', dest='log', action='store_true', default=False)
parser.add_argument('--isAsymptotic', dest='asymp', action='store_true', default=False)
parser.add_argument('--normSM', dest='normSM', action='store_true', default=False)
parser.add_argument('--unblind', dest='unblind', action='store_true', default=False)

opt = parser.parse_args()

'''
-rw-r--r--. 1 rateixei zh 6548 Jan 30 02:26 higgsCombine_Node_SM_CatBased400CTS_qt_central.HybridNew.mH125.quant0.500.root
-rw-r--r--. 1 rateixei zh 6529 Jan 30 03:33 higgsCombine_Node_SM_CatBased400CTS_qt_m1s.HybridNew.mH125.quant0.160.root
-rw-r--r--. 1 rateixei zh 6532 Jan 30 01:36 higgsCombine_Node_SM_CatBased400CTS_qt_p1s.HybridNew.mH125.quant0.840.root
-rw-r--r--. 1 rateixei zh 6537 Jan 30 01:45 higgsCombine_Node_SM_CatBased400CTS_qt_p2s.HybridNew.mH125.quant0.975.root
'''

def main(argv):
  print opt.folders, opt.names
  tdr.setTDRStyle()

  quantiles = ['0.025', '0.160', '0.500', '0.840', '0.975', '-1']
  qt_names = ['m2s', 'm1s', 'central', 'p1s', 'p2s', 'obs']
  gr_1s = TGraphAsymmErrors()
  gr_1s.SetFillColor(kGreen+1)
#  gr_1s.SetLineColor(kBlack)
  gr_2s = TGraphAsymmErrors()
  gr_2s.SetFillColor(kOrange)
#  gr_2s.SetLineColor(kBlack)
  gr_ce = TGraphErrors()
  gr_ce.SetLineColor(kBlue)
  gr_ce.SetMarkerColor(kBlue)
  gr_ce.SetMarkerStyle(24)
  gr_ce.SetMarkerSize(0)
  gr_observed = TGraphErrors()
  gr_observed.SetLineColor(kRed+3)
  gr_observed.SetLineWidth(3)
  gr_observed.SetLineStyle(kDashed)
  gr_observed.SetMarkerColor(kRed+3)
  gr_observed.SetMarkerStyle(20)
  gr_observed.SetMarkerSize(1)
  thisMax = 0

  normalization = 1.
  if opt.normSM: normalization = (33.45*0.0026)

  for iff,ff in enumerate(opt.folders):
    qts = {}
    for iqt, qt in enumerate(quantiles):
      qts[qt] = 0
      fs = glob.glob(ff+"/*"+qt+".root")
      if opt.asymp:
        fs = glob.glob(ff+"*Asymptotic*.root")
      print fs, ff
      if len(fs) > 0:
        tfile = TFile(fs[0], "READ")
        tree = tfile.Get("limit")
        tree.Draw("limit", "quantileExpected>"+str(float(qt)-0.001) + ' && quantileExpected < ' +str(float(qt)+0.001), "goff")
        qts[qt] = tree.GetV1()[0]/(normalization)
#        print qt, qts[qt]
    gr_1s.SetPoint(iff, iff+0.5, qts['0.500'])
    gr_2s.SetPoint(iff, iff+0.5, qts['0.500'])
    gr_ce.SetPoint(iff, iff+0.5, qts['0.500'])
    gr_ce.SetPointError(iff, 0.5, 0)
    e1s = qts['0.500'] - qts['0.160']
    if qts['0.160'] == 0: e1s = 0
    e2s =  qts['0.500'] - qts['0.025']
    if qts['0.025'] == 0: e2s = 0
    gr_1s.SetPointError(iff, 0.5,0.5, e1s, max(qts['0.840'] - qts['0.500'], 0) )
    gr_2s.SetPointError(iff, 0.5,0.5, e2s, max(qts['0.975'] - qts['0.500'], 0) )
    if qts['0.975'] > thisMax: thisMax = qts['0.975']
    gr_observed.SetPoint(iff, iff+0.5, qts['-1'])
    gr_observed.SetPointError(iff, 0.5, 0.)

  c0 = TCanvas('a', 'a', 800, 600)
  c0.SetGrid()
  if opt.log: c0.SetLogy()
  gr_2s.Draw("AE2Z")
  if opt.log:
    gr_2s.GetYaxis().SetRangeUser(0.5/(normalization), thisMax*8)
    gr_2s.GetYaxis().SetLimits(0.5/(normalization), thisMax*8)
  else:
    gr_2s.GetYaxis().SetRangeUser(0., thisMax*1.2)
    gr_2s.GetYaxis().SetLimits(0., thisMax*1.2)
  gr_2s.GetXaxis().CenterLabels()
  if opt.normSM: gr_2s.GetYaxis().SetTitle("#sigma(pp#rightarrowHH) x #it{B}(HH#rightarrowb#bar{b}#gamma#gamma)/#sigma_{SM} x #it{B}_{SM}")
  else: gr_2s.GetYaxis().SetTitle("#sigma(pp#rightarrowHH) x #it{B}(HH#rightarrowb#bar{b}#gamma#gamma) [fb]")
  nbins = gr_2s.GetXaxis().GetNbins()
  rat = float(float(nbins+1.5)/float(len(opt.folders)))
  for iff,ff in enumerate(opt.folders):
#    tbin = gr_2s.GetXaxis().FindBin(iff*2+0.5)
    tbin = int(iff*rat + rat/2.)
    gr_2s.GetXaxis().SetBinLabel(tbin, opt.names[iff])
    gr_2s.GetXaxis().SetNdivisions(10)
    gr_2s.GetXaxis().LabelsOption("hc")
    gr_2s.GetXaxis().SetLabelSize(0.055)
  gr_2s.GetXaxis().SetLimits(0, len(opt.folders))
#  gr_2s.GetXaxis().SetRangeUser(0, len(opt.folders))
  gr_2s.Draw("AE2Z")
  c0.Update()
  gr_1s.Draw("E2Z same")
  gr_ce.Draw("EZ same")
  if opt.unblind: gr_observed.Draw("EPZ same")

  leg = TLegend(0.5, 0.51, 0.89, 0.95)
  headerTitle = "pp#rightarrowHH#rightarrowb#bar{b}#gamma#gamma (SM)"
  leg.SetHeader(headerTitle)
  header =leg.GetListOfPrimitives().First()
  header.SetTextSize(.045)
  leg.SetFillStyle(0)
  leg.SetLineWidth(0)
  leg.SetBorderSize(0)
  leg.AddEntry(gr_observed, 'Observed 95% C.L. upper limit', 'lp')
  leg.AddEntry(gr_ce, 'Expected 95% C.L. upper limit', 'l')
  leg.AddEntry(gr_1s, 'Expected #pm 1 std. dev.', 'f')
  leg.AddEntry(gr_2s, 'Expected #pm 2 std. dev.', 'f')
  leg.SetTextSize(0.035)

  tdr.cmsPrel(float(opt.lumi)*1000,  "13",  0, True,  0, 1.25)
  leg.Draw("same")
  c0.SaveAs("test"+opt.label+"_SM.pdf")

if __name__ == "__main__":
  main(sys.argv[1:])

