from ROOT import *
from array import array
import argparse, os, sys
from HiggsAnalysis.bbggLimits.MyCMSStyle import *
from HiggsAnalysis.bbggLimits.NiceColors import *
from HiggsAnalysis.bbggLimits.DefineScans import *

gROOT.SetBatch()

parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument("-f", "--folder", dest="f", type=str)
parser.add_argument("-o", "--outfile", dest="out", type=str)
parser.add_argument('-l', '--lumi', dest='lumi', type=str, default='35.9')
parser.add_argument("--unblind", dest="unblind", action='store_true', default=False)
opt = parser.parse_args()

outfile = TFile(opt.out, 'RECREATE')

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
# for ikt in scan['kt']:
#  for cg in scan['cg']:
#   for c2 in scan['c2']:
#    for c2g in scan['c2g']:

 kl = klJHEP[ii]
 kt = ktJHEP[ii]
 c2 = c2JHEP[ii]
 cg = cgJHEP[ii]
 c2g = c2gJHEP[ii]

 nodename = 'kl'+str(kl).replace('.', 'p').replace('-', 'm') + '_kt' + str(kt).replace('.', 'p').replace('-', 'm') + '_cg'+ str(cg).replace('.', 'p').replace('-', 'm') + '_c2' + str(c2).replace('.', 'p').replace('-', 'm') + '_c2g' + str(c2g).replace('.', 'p').replace('-', 'm')
 name = opt.f + '/CombinedCard_Node_SM'+nodename  + '/higgsCombineCombinedCard_Node_SM' + nodename + '.Asymptotic.mH125.root'
 print name
 rfile = TFile(name, 'READ')
 if rfile.IsZombie(): continue
 tree = rfile.Get("limit")
 if tree == None or tree.GetEntriesFast() < 1: continue
 for qt in quantiles:
  tree.Draw("limit", "quantileExpected>"+str(float(qt)-0.001) + ' && quantileExpected < ' +str(float(qt)+0.001), "goff")
  lims[qt].append(tree.GetV1()[0])
 rfile.Close()

print plots

for ib in range(0, len(klJHEP)):
  for qt in quantiles:
    print lims['-1']
    if '-1' in qt:
      plots[qt].SetPoint(ib, ib, lims['-1'][ib])
    else:
      plots[qt].SetPoint(ib, ib, lims['0.500'][ib])
  plots['0.160'].SetPointError(ib, 0,0, lims['0.500'][ib] - lims['0.160'][ib], lims['0.840'][ib] - lims['0.500'][ib] )
  plots['0.025'].SetPointError(ib, 0,0, lims['0.500'][ib] - lims['0.025'][ib], lims['0.975'][ib] - lims['0.500'][ib] )
  plots['-1'].SetPointError(ib, myLineWidth,myLineWidth, myLineHeight, myLineHeight )
  plots['0.500'].SetPointError(ib, myLineWidth,myLineWidth, myLineHeight, myLineHeight )


#s2_col = kYellow
s2_col = TColor.GetColor(NiceYellow2)
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
c0.SetGrid()
plots['0.025'].SetMaximum(11)
plots['0.025'].Draw("APZ")
plots['0.025'].SetLineColor(s2_col)
plots['0.025'].SetLineWidth(8)
plots['0.025'].SetTitle("")
plots['0.025'].GetXaxis().SetLimits(-1, len(klJHEP))
plots['0.025'].GetYaxis().SetRangeUser(0, 10)
plots['0.025'].GetXaxis().SetTitle("Shape Benchmark Points")
plots['0.025'].GetYaxis().SetTitle("#sigma(pp#rightarrowHH)#times#it{B}(HH#rightarrowb#bar{b}#gamma#gamma) [fb]")
c0.Update()
plots['0.160'].Draw("EPZ")
plots['0.160'].SetMarkerStyle(21)
plots['0.160'].SetMarkerColor(s1_col)
plots['0.160'].SetLineWidth(7)
plots['0.160'].SetLineColor(s1_col)
c0.Update()
SetPadStyle(c0)
c0.Update()
SetAxisTextSizes(plots['0.025'])
c0.Update()
plots['0.500'].SetLineWidth(2)
plots['0.500'].SetLineColor(cen_col)
plots['0.500'].SetFillStyle(1)
plots['0.500'].SetFillColor(cen_col)
plots['-1'].SetLineWidth(2)
plots['-1'].SetLineColor(kBlack)
plots['-1'].SetFillStyle(1)
plots['-1'].SetFillColor(kBlack)

plots['0.500'].Draw("2 same")
if opt.unblind: plots['-1'].Draw("2 same")

leg = TLegend(0.12, 0.45, 0.45, 0.89)
leg.SetFillColorAlpha(kWhite, 0.8)
leg.SetBorderSize(0)
leg.SetHeader("#font[61]{pp#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}")#{#kappa_{t} = 1, c_{g} = c_{2g} = c_{2} = 0}")
leg.AddEntry(plots['-1'], 'Observed 95% C.L. limit', 'l')
leg.AddEntry(plots['0.500'], 'Expected 95% C.L. limit', 'l')
leg.AddEntry(plots['0.160'], 'Expected #pm 1#sigma', 'l')
leg.AddEntry(plots['0.025'], 'Expected #pm 2#sigma', 'l')
leg.Draw("same")


DrawCMSLabels(c0, '35.9')

c0.SaveAs(opt.out+"NonResPlot.pdf")
c0.SaveAs(opt.out+"NonResPlot.png")

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
