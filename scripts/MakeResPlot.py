#!/usr/bin/env python

from ROOT import *
import sys, getopt, os
import argparse
import glob
#import HiggsAnalysis.bbggLimits.TdrStyle as tdr
#import HiggsAnalysis.bbggLimits.CMS_lumi as CMS_lumi
from HiggsAnalysis.bbggLimits.ResonantCrossSections import *
from HiggsAnalysis.bbggLimits.NiceColors import *
from HiggsAnalysis.bbggLimits.MyCMSStyle import *

gROOT.SetBatch(kTRUE)

parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument('-f', '--inputFolder', dest="folder", default=None, type=str, required=True,
                    help="input folders")
parser.add_argument('-n', '--inputName', dest="name", type=str, default="MyResonantAnalysis",
                    help="input folders cat names")
parser.add_argument('-l', '--lumi', dest='lumi', default='36.5', type=str, help='Integrated luminosoty')
parser.add_argument('--label', dest='label', default='', type=str, help='Label')
parser.add_argument('--log', dest='log', action='store_true', default=False)
parser.add_argument('--isAsymptotic', dest='asymp', action='store_true', default=False)
parser.add_argument('--isGrav', dest='isGrav', action='store_true', default=False)
parser.add_argument('--hmFolder', dest='hmfolder', default=None)
parser.add_argument('--max', dest='max',  default=None, type=float)
parser.add_argument('--min', dest='min',  default=0.001, type=float)
parser.add_argument('--observed', dest='observed', action='store_true', default=False)

opt = parser.parse_args()

'''

-rw-r--r--. 1 rateixei zh 6548 Jan 30 02:26 higgsCombine_Node_SM_CatBased400CTS_qt_central.HybridNew.mH125.quant0.500.root
-rw-r--r--. 1 rateixei zh 6529 Jan 30 03:33 higgsCombine_Node_SM_CatBased400CTS_qt_m1s.HybridNew.mH125.quant0.160.root
-rw-r--r--. 1 rateixei zh 6532 Jan 30 01:36 higgsCombine_Node_SM_CatBased400CTS_qt_p1s.HybridNew.mH125.quant0.840.root
-rw-r--r--. 1 rateixei zh 6537 Jan 30 01:45 higgsCombine_Node_SM_CatBased400CTS_qt_p2s.HybridNew.mH125.quant0.975.root
'''

#masses = [250, 260, 270, 280, 300, 320, 340, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900]
#massesLM = [250, 260, 270, 280, 300, 320, 340, 350, 400, 450, 500, 550]
#massesHM = [550, 600, 650, 700, 750, 800, 900]

legL = TLegend(0.20, 0.65, 0.48, 0.80)
leg = TLegend(0.49, 0.60, 0.87, 0.85)
#leg.SetNColumns(2)
legL.SetFillColorAlpha(kWhite, 0.9)
leg.SetFillColorAlpha(kWhite, 0.9)
def main(argv):
  print opt.folder, opt.name
#  tdr.setTDRStyle()
  #change the CMS_lumi variables (see CMS_lumi.py)
#  CMS_lumi.lumi_13TeV = opt.lumi+" fb^{-1}"
#  CMS_lumi.writeExtraText = 1
#  CMS_lumi.extraText = "Preliminary"
#  CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)


  quantiles = ['0.025', '0.160', '0.500', '0.840', '0.975']
  if opt.observed: quantiles.append('-1')
  qt_names = ['m2s', 'm1s', 'central', 'p1s', 'p2s']
  if opt.observed: qt_names.append('obs')
  gr_1s = TGraphAsymmErrors()
  gr_1s.SetFillColor(cNiceGreen2)
  gr_1s.SetLineColor(cNiceGreen2)
  gr_2s = TGraphAsymmErrors()
  gr_2s.SetFillColor(kOrange)
  gr_2s.SetLineColor(kOrange)
#  gr_2s.SetLineColor(cNiceYellow2)
  gr_ce = TGraphErrors()
  gr_ce.SetLineColor(cNiceBlueDark)
  gr_ce.SetLineWidth(2)
  gr_ce.SetLineStyle(kDashed)
  gr_ce.SetMarkerStyle(24)
  gr_ce.SetMarkerColor(cNiceBlueDark)
  gr_observed = TGraphErrors()
  gr_observed.SetLineColor(kBlack)
  gr_observed.SetLineWidth(2)
  gr_observed.SetMarkerColor(kBlack)
  gr_observed.SetMarkerStyle(20)
  hmgr_1s = TGraphAsymmErrors()
  hmgr_1s.SetFillColor(cNiceGreen2)
  hmgr_1s.SetLineColor(cNiceGreen2)
  hmgr_2s = TGraphAsymmErrors()
  hmgr_2s.SetFillColor(kOrange)
  hmgr_2s.SetLineColor(kOrange)
#  hmgr_2s.SetFillColor(cNiceYellow2)
#  hmgr_2s.SetLineColor(cNiceYellow2)
  hmgr_ce = TGraphErrors()
  hmgr_ce.SetLineColor(cNiceBlueDark)
  hmgr_ce.SetLineWidth(2)
  hmgr_ce.SetLineStyle(kDashed)
  hmgr_ce.SetMarkerStyle(24)
  hmgr_ce.SetMarkerColor(cNiceBlueDark)
  hmgr_observed = TGraphErrors()
  hmgr_observed.SetLineColor(kBlack)
  hmgr_observed.SetLineWidth(2)
  hmgr_observed.SetMarkerColor(kBlack)
  hmgr_observed.SetMarkerStyle(20)
  thisMax = 0
  masses = [260, 270, 280, 300, 320, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900]
  massesLM = [260, 270, 280, 300, 320, 350, 400, 450, 500, 550, 600]
  massesHM = [600, 650, 700, 750, 800, 900]

  if opt.hmfolder: masses = massesLM
  for iff,m in enumerate(masses):
#lims_Res_NewBTagWPRadion_v66/Radion_Node_250_Rad/datacards/higgsCombineRadion_Node_250_Rad_qt_
    ff = opt.folder.replace("MASS", str(m))
    qts = {}
    for iqt, qt in enumerate(quantiles):
      qts[qt] = 0
      fs = glob.glob(ff+"*"+qt+".root")
      if opt.asymp:
        fs = glob.glob(ff+"*Asymptotic*.root")
      print fs, ff
      if len(fs) > 0:
        tfile = TFile(fs[0], "READ")
        tree = tfile.Get("limit")
        tree.Draw("limit", "quantileExpected>"+str(float(qt)-0.001) + ' && quantileExpected < ' +str(float(qt)+0.001), "goff")
        qts[qt] = tree.GetV1()[0]
        print qt, qts[qt]
    gr_1s.SetPoint(iff, m, qts['0.500'])
    gr_2s.SetPoint(iff, m, qts['0.500'])
    gr_ce.SetPoint(iff, m, qts['0.500'])
    gr_ce.SetPointError(iff, 0, 0)
    e1s = qts['0.500'] - qts['0.160']
    if qts['0.160'] == 0: e1s = 0
    e2s =  qts['0.500'] - qts['0.025']
    if qts['0.025'] == 0: e2s = 0
    gr_1s.SetPointError(iff, 0,0, e1s, max(qts['0.840'] - qts['0.500'], 0) )
    gr_2s.SetPointError(iff, 0,0, e2s, max(qts['0.975'] - qts['0.500'], 0) )
    if opt.observed: gr_observed.SetPoint(iff, m, qts['-1'])
    if qts['0.975'] > thisMax: thisMax = qts['0.975']

  if opt.hmfolder:
    for iff,m in enumerate(massesHM):
      ff = opt.hmfolder.replace("MASS", str(m))
      qts = {}
      for iqt, qt in enumerate(quantiles):
        qts[qt] = 0
        fs = glob.glob(ff+"*"+qt+".root")
        if opt.asymp:
          fs = glob.glob(ff+"*Asymptotic*.root")
        print fs, ff
        if len(fs) > 0:
          tfile = TFile(fs[0], "READ")
          tree = tfile.Get("limit")
          tree.Draw("limit", "quantileExpected>"+str(float(qt)-0.001) + ' && quantileExpected < ' +str(float(qt)+0.001), "goff")
          qts[qt] = tree.GetV1()[0]
          print qt, qts[qt]
      hmgr_1s.SetPoint(iff, m, qts['0.500'])
      hmgr_2s.SetPoint(iff, m, qts['0.500'])
      hmgr_ce.SetPoint(iff, m, qts['0.500'])
      hmgr_ce.SetPointError(iff, 0, 0)
      e1s = qts['0.500'] - qts['0.160']
      if qts['0.160'] == 0: e1s = 0
      e2s =  qts['0.500'] - qts['0.025']
      if qts['0.025'] == 0: e2s = 0
      hmgr_1s.SetPointError(iff, 0,0, e1s, max(qts['0.840'] - qts['0.500'], 0) )
      hmgr_2s.SetPointError(iff, 0,0, e2s, max(qts['0.975'] - qts['0.500'], 0) )
      if qts['0.975'] > thisMax: thisMax = qts['0.975']
      if opt.observed: hmgr_observed.SetPoint(iff, m, qts['-1'])


  SetGeneralStyle()
  c0 = TCanvas('a', 'a', 800, 600)
#  c0.SetGrid()
#  SetFrameStyle(c0)
  SetPadStyle(c0)
  if opt.log: c0.SetLogy()
  gr_2s.Draw("A3Z")
  gr_2s.GetYaxis().SetRangeUser(0, thisMax*1.2)
  gr_2s.GetYaxis().SetTitle("#sigma(pp#rightarrowX) #times #it{B}(X#rightarrowHH#rightarrow#gamma#gammab#bar{b}) [fb]")
  gr_2s.GetXaxis().SetTitle("Resonance Mass [GeV]")
  myMin = 0
  myMax = thisMax*1.2
  if opt.max :
    gr_2s.GetYaxis().SetRangeUser(float(opt.min), float(opt.max))
    myMin = float(opt.min)
    myMax = float(opt.max)
  gr_2s.GetXaxis().SetLimits(245, 905)
  SetAxisTextSizes(gr_2s)
#  gr_2s.Draw("A3Z")
#  CMS_lumi.CMS_lumi(c0, 0,11)
#  c0.cd()
#  c0.Update()
#  c0.RedrawAxis()
#  frame = c0.GetFrame()
#  frame.Draw()
#  c0.Update()
  gr_1s.Draw("3Z same")
  gr_ce.Draw("LZP same")
  if opt.observed: gr_observed.Draw("PL same")
  if opt.hmfolder:
    hmgr_2s.Draw("3Z same")
    hmgr_1s.Draw("3Z same")
    hmgr_ce.Draw("LZP same")
    gr_ce.Draw("LZP same")
    line = TLine(600, myMin, 600, 6)
    line.SetLineStyle(kDashed)
    line.SetLineWidth(2)
    line.Draw("same")
    if opt.observed:
      hmgr_observed.Draw("PL same")
      gr_observed.Draw("PL same")
#  tdrStyle.SetTitleSize(0.003, "XYZ")
#  tdrStyle.SetTitleYOffset(.8)
  c0.Update()
  c0.cd()
  c0.Update()
  c0.RedrawAxis()
#  c0.SetGrid()

  headerTitleL = "#font[61]{pp#rightarrowX#rightarrowHH#rightarrow#gamma#gammab#bar{b} (spin-0)}\n"
  if opt.isGrav: headerTitleL = "#font[61]{pp#rightarrowX#rightarrowHH#rightarrow#gamma#gammab#bar{b} (spin-2)}\n"

  ltx = TLatex(260, 150, headerTitleL)
  ltx.SetTextSize(.035)
  ltx.DrawLatex(280, myMax*0.6, headerTitleL)

  legL.SetLineWidth(0)
  legL.SetBorderSize(0)

  
  headerTitle = "95% CL upper limits"
  leg.SetHeader(headerTitle)
  header = leg.GetListOfPrimitives().First()
  header.SetTextSize(.035)
#  leg.SetFillStyle(0)


  leg.SetLineWidth(0)
  leg.SetBorderSize(0)
  leg.AddEntry(gr_observed, 'Observed', 'lp')
  leg.AddEntry(gr_ce, 'Expected', 'lp')
  leg.AddEntry(gr_1s, 'Expected #pm 1 std. deviation', 'f')
  leg.AddEntry(gr_2s, 'Expected #pm 2 std. deviation', 'f')
  leg.SetTextSize(0.035)
  if opt.isGrav:
    gr_grav_0p5kpl.SetLineColor(cNicePurple)
    gr_grav_1p0kpl.SetLineColor(cNiceRed)
    gr_grav_0p5kpl.Draw("same LZ")
    gr_grav_1p0kpl.Draw("same LZ")
    headerTitle = "Bulk KK-graviton"
    legL.SetHeader(headerTitle)
    legL.SetTextSize(.035)
    legL.AddEntry(gr_grav_0p5kpl, grav_0p5kpl_leg, 'l')
    legL.AddEntry(gr_grav_1p0kpl, grav_1p0kpl_leg, 'l')
  else:
    gr_rad_2tev.SetLineColor(cNicePurple)
    gr_rad_3tev.SetLineColor(cNiceRed)
    gr_rad_2tev.Draw("same LZ")
    gr_rad_3tev.Draw("same LZ")
    headerTitle = "Bulk radion"
    legL.SetHeader(headerTitle)
    legL.SetTextSize(.035)
    legL.AddEntry(gr_rad_2tev, rad_2tev_leg, 'l')
    legL.AddEntry(gr_rad_3tev, rad_3tev_leg, 'l')
 
  legL.SetTextSize(0.035)

  leg.Draw()
  legL.Draw()

#  tdr.cmsPrel(float(opt.lumi)*1000,  "13",  0, True,  0, 1.25)
  DrawCMSLabels(c0, '35.9')
  c0.SaveAs(opt.name+".pdf")

if __name__ == "__main__":
  main(sys.argv[1:])

