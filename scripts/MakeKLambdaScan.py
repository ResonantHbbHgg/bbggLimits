from ROOT import *
import argparse, os
from HiggsAnalysis.bbggLimits.NiceColors import *
import HiggsAnalysis.bbggLimits.CMS_lumi as CMS_lumi
from HiggsAnalysis.bbggLimits.MyCMSStyle import *
from HiggsAnalysis.bbggLimits.DefineScans import *

gROOT.SetBatch()

A13tev = [2.09078, 10.1517, 0.282307, 0.101205, 1.33191, -8.51168, -1.37309, 2.82636, 1.45767, -4.91761, -0.675197, 1.86189, 0.321422, -0.836276, -0.568156]

def functionGF(kl,kt,c2,cg,c2g,A): return A[0]*kt**4 + A[1]*c2**2 + (A[2]*kt**2 + A[3]*cg**2)*kl**2  + A[4]*c2g**2 + ( A[5]*c2 + A[6]*kt*kl )*kt**2  + (A[7]*kt*kl + A[8]*cg*kl )*c2 + A[9]*c2*c2g  + (A[10]*cg*kl + A[11]*c2g)*kt**2+ (A[12]*kl*cg + A[13]*c2g )*kt*kl + A[14]*cg*c2g*kl


parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument("-f", "--folder", dest="f", type=str)
parser.add_argument("--unblind", dest="unblind", action='store_true', default=False)
parser.add_argument('-o', '--outFile', dest='outf', type=str)
opt = parser.parse_args()

outfile = TFile(opt.outf, 'RECREATE')

quantiles = ['0.025', '0.160', '0.500', '0.840', '0.975', '-1']
lims = {}
plots = {}
for qt in quantiles:
  lims[qt] = []
  plots[qt] = TGraphAsymmErrors()
  plots[qt].SetName('plot_'+qt.replace('.', 'p').replace("-", "m"))

'''
        tfile = TFile(fs[0], "READ")
        tree = tfile.Get("limit")
        tree.Draw("limit", "quantileExpected>"+str(float(qt)-0.001) + ' && quantileExpected < ' +str(float(qt)+0.001), "goff")
        qts[qt] = tree.GetV1()[0]/(normalization)
'''

myKl = []
notworked = open('klscan_notworked.txt', 'w+')
for kl in scan_kl['kl']:
#  fname = opt.f + '/HighMass_Node_SMkl' + str(kl).replace('.', 'p').replace('-', 'm') + '_kt1p0_cg0p0_c20p0_c2g0p0/datacards/higgsCombineHighMass_Node_SMkl' + str(kl).replace('.', 'p').replace('-', 'm') + '_kt1p0_cg0p0_c20p0_c2g0p0.Asymptotic.mH125.root'
  fname = opt.f + '/CombinedCard_Node_SMkl' + str(kl).replace('.', 'p').replace('-', 'm') + '_kt1p0_cg0p0_c20p0_c2g0p0/higgsCombineCombinedCard_Node_SMkl' + str(kl).replace('.', 'p').replace('-', 'm') + '_kt1p0_cg0p0_c20p0_c2g0p0.Asymptotic.mH125.root'
  tfile = TFile(fname, "READ")
  if tfile.IsZombie() == 1:
    notworked.write(str(kl) + " 1.0 0.0 0.0 0.0 \n")
    continue
  tree = tfile.Get('limit')
  if tree == None:
    notworked.write(str(kl) + " 1.0 0.0 0.0 0.0 \n")
    continue
  myKl.append(kl)
  for qt in quantiles:
    tree.Draw("limit", "quantileExpected>"+str(float(qt)-0.001) + ' && quantileExpected < ' +str(float(qt)+0.001), "goff")
    lims[qt].append(tree.GetV1()[0])
  if lims['0.500'] == 0 or lims['-1'] == 0:
    notworked.write(str(kl) + " 1.0 0.0 0.0 0.0 \n")
    continue
  tfile.Close()

nonresXSEC = TGraphAsymmErrors()
nonresXSEC.SetName("nonresXsec")
nonresXSEC_2 = TGraphAsymmErrors()
nonresXSEC_2.SetName("nonresXsec_2")
for ikl,kl in enumerate(myKl):
  xsec = 33.45*0.0026*functionGF(kl, 1.0, 0.0, 0.0, 0.0, A13tev)
  xsec_2 = 33.45*0.0026*functionGF(kl*2.0, 2.0, 0.0, 0.0, 0.0, A13tev)
  nonresXSEC.SetPoint(ikl, kl, xsec)
  nonresXSEC.SetPointError(ikl, 0,0, xsec*0.1, xsec*0.1)
  nonresXSEC_2.SetPoint(ikl, kl, xsec_2)
  nonresXSEC_2.SetPointError(ikl, 0,0, xsec_2*0.1, xsec_2*0.1)
  for qt in quantiles:
    if '-1' in qt:
      plots[qt].SetPoint(ikl, kl, lims['-1'][ikl])
    else:
      plots[qt].SetPoint(ikl, kl, lims['0.500'][ikl])
  plots['0.160'].SetPointError(ikl, 0,0, lims['0.500'][ikl] - lims['0.160'][ikl], lims['0.840'][ikl] - lims['0.500'][ikl] ) 
  plots['0.025'].SetPointError(ikl, 0,0, lims['0.500'][ikl] - lims['0.025'][ikl], lims['0.975'][ikl] - lims['0.500'][ikl] ) 


SetGeneralStyle()
c = TCanvas("c", "c", 800, 600)
c.SetGrid()
SetPadStyle(c)
#s2_col = kYellow
s2_col = TColor.GetColor(NiceYellow2)
#s1_col = kGreen
s1_col = TColor.GetColor(NiceGreen2)
#th_col = kRed
th_col = TColor.GetColor(NiceRed)
th2_col = TColor.GetColor(NiceBlue)
ob_col = kBlack
#ob_col = TColor.GetColor(NiceBlueDark)

plots['0.025'].SetFillColor(s2_col)
plots['0.025'].SetLineColor(s2_col)
plots['0.160'].SetFillColor(s1_col)
plots['0.160'].SetLineColor(s1_col)
plots['0.500'].SetLineColor(cNiceBlueDark)
plots['0.500'].SetLineWidth(2)
plots['0.500'].SetLineStyle(kDashed)
plots['-1'].SetLineColor(ob_col)
plots['-1'].SetLineWidth(3)
nonresXSEC.SetLineColor(th_col)
nonresXSEC.SetFillColorAlpha(th_col, 0.5)
nonresXSEC_2.SetLineColor(cNicePurple)
nonresXSEC_2.SetFillColorAlpha(cNicePurple, 0.5)
plots['0.025'].Draw("A3Z")

plots['0.025'].GetXaxis().SetTitle('#kappa_{#lambda}/#kappa_{t}')
#plots['0.025'].GetYaxis().SetTitle('95% C.L. limit on #sigma(pp#rightarrowHH)#times#font[52]{B}(HH#rightarrowb#bar{b}#gamma#gamma) [fb]')
plots['0.025'].GetYaxis().SetTitle('#sigma(pp#rightarrowHH)#times#font[52]{B}(HH#rightarrowb#bar{b}#gamma#gamma) [fb]')
plots['0.025'].SetMaximum(14)
plots['0.025'].GetXaxis().SetRangeUser(-20, 20)
SetAxisTextSizes(plots['0.025'])
c.Update()

plots['0.160'].Draw("3Z same")
plots['0.500'].Draw("LZ same")
if(opt.unblind): plots['-1'].Draw("CZ same")
nonresXSEC.Draw("3ZL same")
nonresXSEC_2.Draw("3ZL same")

leg = TLegend(0.12, 0.45, 0.45, 0.89)
leg.SetFillColorAlpha(kWhite, 0.8)
leg.SetBorderSize(0)
leg.SetHeader("#font[61]{pp#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}")#{#kappa_{t} = 1, c_{g} = c_{2g} = c_{2} = 0}")
leg.AddEntry(plots['-1'], 'Observed 95% C.L. limit', 'l')
leg.AddEntry(plots['0.500'], 'Expected 95% C.L. limit', 'l')
leg.AddEntry(plots['0.160'], 'Expected #pm 1#sigma', 'f')
leg.AddEntry(plots['0.025'], 'Expected #pm 2#sigma', 'f')
leg.AddEntry(nonresXSEC, 'Theory Prediction (#kappa_{t} = 1)', 'f')
leg.AddEntry(nonresXSEC_2, 'Theory Prediction (#kappa_{t} = 2)', 'f')
leg.Draw("same")

#CMS_lumi.lumi_13TeV = "35.87 fb^{-1}"
#CMS_lumi.writeExtraText = 1
#CMS_lumi.extraText = "Preliminary"
#CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
#CMS_lumi.CMS_lumi(c, 0,11)

latex = TLatex()
latex.SetNDC()
latex.SetTextSize(25)
latex.SetTextAlign(33)
latex.SetTextFont(43)
#latex.DrawLatex(0.89,.89,"#kappa_{t} = 1, c_{g} = c_{2g} = c_{2} = 0")
latex.DrawLatex(0.89,.89,"c_{g} = c_{2g} = c_{2} = 0")

c.Update()
c.cd()
c.Update()
c.RedrawAxis()

DrawCMSLabels(c, '35.9')

c.SaveAs(opt.outf.replace(".root", "")+".pdf")

outfile.cd()
for qt in quantiles:
  plots[qt].Write()

nonresXSEC.Write()
nonresXSEC_2.Write()
outfile.Close()
