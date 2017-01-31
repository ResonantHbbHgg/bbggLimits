#!/usr/bin/env python

from ROOT import *
import sys, getopt, os
import argparse
import glob

parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument('-f', '--inputFolder', dest="folder", default=None, type=str, required=True,
                    help="input folders")
parser.add_argument('-n', '--inputName', dest="name", default=None, type=str, required=True,
                    help="input folders cat names")
parser.add_argument('-l', '--lumi', dest='lumi', default='36.5', type=str, help='Integrated luminosoty')
parser.add_argument('--label', dest='label', default='', type=str, help='Label')
parser.add_argument('--log', dest='log', action='store_true', default=False)

opt = parser.parse_args()

'''
-rw-r--r--. 1 rateixei zh 6548 Jan 30 02:26 higgsCombine_Node_SM_CatBased400CTS_qt_central.HybridNew.mH125.quant0.500.root
-rw-r--r--. 1 rateixei zh 6529 Jan 30 03:33 higgsCombine_Node_SM_CatBased400CTS_qt_m1s.HybridNew.mH125.quant0.160.root
-rw-r--r--. 1 rateixei zh 6532 Jan 30 01:36 higgsCombine_Node_SM_CatBased400CTS_qt_p1s.HybridNew.mH125.quant0.840.root
-rw-r--r--. 1 rateixei zh 6537 Jan 30 01:45 higgsCombine_Node_SM_CatBased400CTS_qt_p2s.HybridNew.mH125.quant0.975.root
'''

masses = [250, 260, 270, 280, 300, 320, 340, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900]

def main(argv):
  print opt.folder, opt.name

  quantiles = ['0.025', '0.160', '0.500', '0.840', '0.975']
  qt_names = ['m2s', 'm1s', 'central', 'p1s', 'p2s']
  gr_1s = TGraphAsymmErrors()
  gr_1s.SetFillColor(kGreen)
  gr_1s.SetLineColor(kGreen)
  gr_2s = TGraphAsymmErrors()
  gr_2s.SetFillColor(kYellow)
  gr_2s.SetLineColor(kYellow)
  gr_ce = TGraphErrors()
  gr_ce.SetLineColor(kBlack)
  gr_ce.SetMarkerColor(kBlack)
  thisMax = 0
  for iff,m in enumerate(masses):
#lims_Res_NewBTagWPRadion_v66/Radion_Node_250_Rad/datacards/higgsCombineRadion_Node_250_Rad_qt_
    ff = opt.folder.replace("MASS", str(m))
    qts = {}
    for iqt, qt in enumerate(quantiles):
      qts[qt] = 0
      fs = glob.glob(ff+"*"+qt+".root")
      print fs, ff
      if len(fs) > 0:
        tfile = TFile(fs[0], "READ")
        tree = tfile.Get("limit")
        tree.Draw("limit", "", "goff")
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
    if qts['0.975'] > thisMax: thisMax = qts['0.975']

  c0 = TCanvas('a', 'a', 800, 600)
  c0.SetGrid()
  if opt.log: c0.SetLogy()
  gr_2s.Draw("A3Z")
  gr_2s.GetYaxis().SetRangeUser(0, thisMax*1.2)
#  gr_2s.GetXaxis().CenterLabels()
#  for iff,ff in enumerate(opt.folders):
#    tbin = gr_2s.GetXaxis().FindBin(iff+0.5)
#    gr_2s.GetXaxis().SetBinLabel(tbin, opt.names[iff])
#    gr_2s.GetXaxis().SetNdivisions(10)
#    gr_2s.GetXaxis().LabelsOption("hc")
#    gr_2s.GetXaxis().SetLabelSize(0.055)
#  gr_2s.GetXaxis().SetLimits(0, len(opt.folders))
#  gr_2s.GetXaxis().SetRangeUser(0, len(opt.folders))
  gr_2s.Draw("A3Z")
  c0.Update()
  gr_1s.Draw("3Z same")
  gr_ce.Draw("LZ same")
  c0.SaveAs("test.pdf")

if __name__ == "__main__":
  main(sys.argv[1:])

