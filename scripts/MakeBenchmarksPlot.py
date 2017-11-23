#!/usr/bin/env python

from ROOT import *
from array import array
import argparse, os, sys
from HiggsAnalysis.bbggLimits.MyCMSStyle import *
from HiggsAnalysis.bbggLimits.NiceColors import *
from HiggsAnalysis.bbggLimits.DefineScans import *

gROOT.SetBatch()
gStyle.SetOptStat(0)

parser =  argparse.ArgumentParser(description='Benchmark plot maker')
parser.add_argument("limdir")
parser.add_argument('-u',"--unblind", dest="unblind", action='store_true', default=False)
opt = parser.parse_args()


myLineHeight = 0.02
myLineWidth = 0.05

quantiles = ['0.025', '0.160', '0.500', '0.840', '0.975', '-1']
lims = {}
plots = {}
for qt in quantiles:
  lims[qt] = []
  plots[qt] = TGraphAsymmErrors()
  plots[qt].SetName('plot_'+qt.replace('.', 'p').replace("-", "m"))


for ii in xrange(0, len(klJHEP)):
  kl = klJHEP[ii]
  kt = ktJHEP[ii]
  c2 = c2JHEP[ii]
  cg = cgJHEP[ii]
  c2g = c2gJHEP[ii]

  # Use this one for Rafael's results: 
  #nodename = 'Node_SM'+'_'.join(['kl'+str(kl), 'kt' + str(kt), 'cg'+ str(cg), 'c2' + str(c2), 'c2g' + str(c2g)]).replace('.', 'p').replace('-', 'm')

  nodename = "_".join(['ARW','kl',str(float(kl)),'kt',str(float(kt)),'cg',str(float(cg)),'c2',str(float(c2)),'c2g',str(float(c2g))]).replace('.', 'p').replace('-', 'm')
  name = opt.limdir + '/CombinedCard_'+nodename  + '/higgsCombine_' + nodename + '.Asymptotic.mH125_1.root'

  print name
  
  rfile = TFile(name, 'READ')
  if rfile.IsZombie(): continue
  tree = rfile.Get("limit")
  if tree == None or tree.GetEntriesFast() < 1: continue
  for qt in quantiles:
    tree.Draw("limit", "quantileExpected>"+str(float(qt)-0.001) + ' && quantileExpected < ' +str(float(qt)+0.001), "goff")
    lims[qt].append(tree.GetV1()[0])
  rfile.Close()
  print 'Expected limit at r=0.5: ', lims['0.500']
  # print plots

h1 = TH1F('h1', '', 14, 0.5, 14.5)
h1.GetXaxis().SetTitle("Shape benchmark hypothesis")
h1.GetYaxis().SetTitle("#sigma(pp#rightarrowHH) #times #it{B}(HH#rightarrow#gamma#gammab#bar{b}) [fb]")
h1.SetMaximum(10)

for ib in range(0, len(klJHEP)):

  ibin = ib
  if ib == 0: ibin = 13
  if ib == 13: ibin = 14

  for qt in quantiles:
    # print lims['-1']
    if '-1' in qt:
      plots[qt].SetPoint(ib, ibin, lims['-1'][ib])
    else:
      plots[qt].SetPoint(ib, ibin, lims['0.500'][ib])
  
  plots['0.160'].SetPointError(ib, 0,0, lims['0.500'][ib] - lims['0.160'][ib], lims['0.840'][ib] - lims['0.500'][ib] )
  plots['0.025'].SetPointError(ib, 0,0, lims['0.500'][ib] - lims['0.025'][ib], lims['0.975'][ib] - lims['0.500'][ib] )
  plots['-1'].SetPointError(ib, myLineWidth,myLineWidth, myLineHeight, myLineHeight )
  plots['0.500'].SetPointError(ib, myLineWidth,myLineWidth, myLineHeight, myLineHeight )
  thisbin = h1.FindBin(float(ibin - 0.01))
  h1.GetXaxis().SetBinLabel(thisbin, str(ibin))

smbin = h1.FindBin(12.99)
topbin = h1.FindBin(13.99)
h1.GetXaxis().SetBinLabel(smbin, 'SM')
h1.GetXaxis().SetBinLabel(topbin, '#kappa_{#lambda}= 0')
SetAxisTextSizes(h1)
h1.GetXaxis().SetLabelSize(0.052)

#s2_col = kYellow
s2_col = kOrange
#s1_col = kGreen
s1_col = TColor.GetColor(NiceGreen2)
#th_col = kRed
th_col = TColor.GetColor(NiceRed)
th2_col = TColor.GetColor(NiceOrange)
ob_col = kBlack
#ob_col = TColor.GetColor(NiceBlueDark)
cen_col = cNiceBlueDark

SetGeneralStyle()
c0 = TCanvas("c", "c", 800, 600)
#c0.SetGrid()
h1.Draw()
#plots['0.025'].SetMaximum(11)
plots['0.025'].Draw("PZ same")
plots['0.025'].SetLineColor(s2_col)
plots['0.025'].SetFillColor(s2_col)
plots['0.025'].SetLineWidth(10)
#plots['0.025'].SetTitle("")
#plots['0.025'].GetXaxis().SetLimits(-1, len(klJHEP))
#plots['0.025'].GetYaxis().SetRangeUser(0, 10)
#plots['0.025'].GetXaxis().SetTitle("Shape Benchmark Points")
#plots['0.025'].GetYaxis().SetTitle("#sigma(pp#rightarrowHH)#times#it{B}(HH#rightarrowb#bar{b}#gamma#gamma) [fb]")
c0.Update()
plots['0.160'].Draw("EPZ same")
#plots['0.160'].SetMarkerStyle(21)
plots['0.160'].SetMarkerColor(s1_col)
plots['0.160'].SetLineWidth(10)
plots['0.160'].SetLineColor(s1_col)
plots['0.160'].SetFillColor(s1_col)
c0.Update()
SetPadStyle(c0)
c0.Update()
#SetAxisTextSizes(plots['0.025'])
c0.Update()
plots['0.500'].SetLineWidth(2)
plots['0.500'].SetLineColor(cen_col)
plots['0.500'].SetFillStyle(1)
plots['0.500'].SetFillColor(cen_col)
plots['0.500'].SetMarkerColor(cen_col)
plots['0.500'].SetMarkerStyle(24)
plots['0.500'].SetMarkerSize(1.1)
plots['-1'].SetLineWidth(2)
plots['-1'].SetLineColor(kBlack)
plots['-1'].SetFillStyle(1)
plots['-1'].SetFillColor(kBlack)
plots['-1'].SetMarkerColor(kBlack)
plots['-1'].SetMarkerStyle(20)
plots['-1'].SetMarkerSize(1.1)

plots['0.500'].Draw("P same")
if opt.unblind: plots['-1'].Draw("P same")


ltx = TLatex()
ltx.SetNDC()
ltx.SetTextSize(.035)
ltx.DrawLatex(0.14,0.83,'#font[61]{pp#rightarrowHH#rightarrow#gamma#gammab#bar{b}}')


leg = TLegend(0.50, 0.6, 0.80, 0.87)
leg.SetFillColorAlpha(kWhite, 0.8)
leg.SetBorderSize(0)
headerTitle = "95% CL upper limits"
leg.SetHeader(headerTitle)
header = leg.GetListOfPrimitives().First()
leg.SetTextSize(.035)
leg.AddEntry(plots['-1'], 'Observed', 'p')
leg.AddEntry(plots['0.500'], 'Expected', 'p')
leg.AddEntry(plots['0.160'], 'Expected #pm 1 std. deviation', 'l')
leg.AddEntry(plots['0.025'], 'Expected #pm 2 std. deviation', 'l')
leg.Draw("same")


DrawCMSLabels(c0, '35.9')

c0.SaveAs(opt.limdir+"/BenchmarkPlot.pdf")
c0.SaveAs(opt.limdir+"/BenchmarkPlot.png")


plots['0.500'].Print("all")

sys.exit()

'''
plots['0.160']_hh = TGraphAsymmErrors(len(nodesList), array('d', nodesList), array('d', [i/2.6 for i in centralVal]), array('d', zeros), array('d', zeros), array('d', [i/2.6 for i in s1_down]), array('d', [i/2.6 for i in s1_up]))
plots['0.025']_hh = TGraphAsymmErrors(len(nodesList), array('d', nodesList), array('d', [i/2.6 for i in centralVal]), array('d', zeros), array('d', zeros), array('d', [i/2.6 for i in s2_down]), array('d', [i/2.6 for i in s2_up]))

c0 = TCanvas("c", "c", 800, 600)
plots['0.025']_hh.SetMaximum(Maximum*2/2.6)
plots['0.025']_hh.Draw("APZ")
plots['0.025']_hh.SetLineColor(11111)
plots['0.025']_hh.SetLineWidth(8)
plots['0.025']_hh.SetTitle("")
plots['0.025']_hh.GetXaxis().SetLimits(-1, 14)
plots['0.025']_hh.GetXaxis().SetTitle("Benchmark Points")
plots['0.025']_hh.GetYaxis().SetTitle("#sigma(pp#rightarrowHH)/Br(HH#rightarrowbb#gamma#gamma) [pb]")
c0.Update()
plots['0.160']_hh.Draw("EPZ")
plots['0.160']_hh.SetMarkerStyle(21)
plots['0.160']_hh.SetMarkerColor(kBlue+1)
plots['0.160']_hh.SetLineWidth(7)
plots['0.160']_hh.SetLineColor(kGreen+1)

tlatex.SetTextAngle(0)
tlatex.SetTextColor(kBlack)
tlatex.SetTextFont(63)
tlatex.SetTextAlign(11)
tlatex.SetTextSize(25)
tlatex.DrawLatex(0.11, 0.91, "CMS")
tlatex.SetTextFont(53)
tlatex.DrawLatex(0.18, 0.91, "Preliminary")
tlatex.SetTextFont(43)
tlatex.DrawLatex(0.6, 0.91, "#sqrt{s} = 13 TeV, L = " + str(opt.lumi) + " fb^{-1}")
leg.Draw("same")
c0.SaveAs(folder+"/NonResPlot_hh.pdf")
c0.SaveAs(folder+"/NonResPlot_hh.png")
'''
