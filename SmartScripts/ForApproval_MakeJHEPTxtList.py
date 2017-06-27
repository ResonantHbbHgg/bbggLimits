from ROOT import *
import argparse, os
from HiggsAnalysis.bbggLimits.DefineScans import *

folder = ["LIMS_LT_350_HMHPC_970_HMMPC_600_LMHPC_985_LMMPC_600_v66"]
outf = "JHEP_Text_LIMS_LT_350_HMHPC_970_HMMPC_600_LMHPC_985_LMMPC_600_v66"


quantiles = ['0.025', '0.160', '0.500', '0.840', '0.975', '-1']
lims = {}
for qt in quantiles:
  lims[qt] = []

myKl = []

worked = open(outf+".txt", "w+")
not_worked = open(outf + "_not_worked.txt", "w+")

counter = 0
#npts = 0

myScan = scan_kl

for ff in folder:
  for ii in range(0, len(klJHEP)):
    kl = klJHEP[ii]
    kt = ktJHEP[ii]
    cg = cgJHEP[ii]
    c2 = c2JHEP[ii]
    c2g =c2gJHEP[ii]
#    if npts > 50: break 
#    npts += 1
    thisName = 'kl' + str(kl).replace('.', 'p').replace('-', 'm') + '_kt' + str(kt).replace('.', 'p').replace('-', 'm') + '_cg'+str(cg).replace('.', 'p').replace('-', 'm')+'_c2'+str(c2).replace('.', 'p').replace('-', 'm')+'_c2g'+str(c2g).replace('.', 'p').replace('-', 'm')
    fname_HM = ff + '/HighMass_Node_SM' + thisName + '/datacards/higgsCombineHighMass_Node_SM' + thisName + '.Asymptotic.mH125.root'
    fname_LM = ff + '/LowMass_Node_SM' + thisName + '/datacards/higgsCombineLowMass_Node_SM' + thisName + '.Asymptotic.mH125.root'
    fname = ff + '/CombinedCard_Node_SM' + thisName + '/higgsCombineCombinedCard_Node_SM'+ thisName + '.Asymptotic.mH125.root'
    tfile = TFile(fname, "READ")
    tfile_HM = TFile(fname_HM, "READ")
    tfile_LM = TFile(fname_LM, "READ")
    thesePoints = str(kl) + ' ' + str(kt) + ' ' + str(cg) + ' ' + str(c2) + ' ' + str(c2g)
    if tfile.IsZombie() == 1 or tfile_HM.IsZombie() == 1 or tfile_LM.IsZombie() == 1:
      not_worked.write(thesePoints+'\t isZombie \n')
      continue
    tree = tfile.Get('limit')
    tree_LM = tfile.Get('limit')
    tree_HM = tfile.Get('limit')
    if tree == None or tree_LM == None or tree_HM == None:
      not_worked.write(thesePoints+'\t tree is None \n')
      continue
    if tree.GetEntries() < 2 or tree_LM.GetEntries() < 2 or tree_HM.GetEntries() < 2:
      not_worked.write(thesePoints+'\t tree not filled correctly \n')
      continue
    myKl.append(kl)
    for qt in quantiles:
      tree.Draw("limit", "quantileExpected>"+str(float(qt)-0.001) + ' && quantileExpected < ' +str(float(qt)+0.001), "goff")
      lims[qt].append(tree.GetV1()[0])
    tfile.Close()
    if lims['0.500'][counter] == 0 or lims['-1'][counter] == 0:
      print '#######################################################', lims['0.500'][counter], lims['-1'][counter]
      not_worked.write(thesePoints+'\t Expected or Observed limit is 0 \n')
      continue
# ipt kl kt exp obs exp_p1s exp_m1s exp_p2s exp_m2s
    str_worked = str(counter) + ' '
    str_worked += str(kl) + ' '
    str_worked += str(kt) + ' '
    str_excl = str_worked
    str_worked += str(lims['0.500'][counter]) + ' '
    str_worked += str(lims['-1'][counter]) + ' ' 
    str_worked += str(lims['0.840'][counter]) + ' ' 
    str_worked += str(lims['0.160'][counter]) + ' '
    str_worked += str(lims['0.975'][counter]) + ' ' 
    str_worked += str(lims['0.025'][counter]) + '\n'
    worked.write(str_worked)

    counter += 1

worked.close()
not_worked.close()
'''
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


c = TCanvas("c", "c", 800, 600)
c.SetGrid()
#s2_col = kYellow
s2_col = TColor.GetColor(NiceYellow)
#s1_col = kGreen
s1_col = TColor.GetColor(NiceGreen)
#th_col = kRed
th_col = TColor.GetColor(NiceRed)
th2_col = TColor.GetColor(NiceOrange)
#ob_col = kBlue
ob_col = TColor.GetColor(NiceBlueDark)

plots['0.025'].SetFillColor(s2_col)
plots['0.025'].SetLineColor(s2_col)
plots['0.160'].SetFillColor(s1_col)
plots['0.160'].SetLineColor(s1_col)
plots['0.500'].SetLineColor(kBlack)
plots['0.500'].SetLineWidth(1)
plots['0.500'].SetLineStyle(kDashed)
plots['-1'].SetLineColor(ob_col)
plots['-1'].SetLineWidth(3)
nonresXSEC.SetLineColor(th_col)
nonresXSEC.SetFillColorAlpha(th_col, 0.5)
nonresXSEC_2.SetLineColor(th_col)
nonresXSEC_2.SetFillColorAlpha(th2_col, 0.5)
plots['0.025'].Draw("A3Z")

plots['0.025'].GetXaxis().SetTitle('#kappa_{#lambda}/#kappa_{t}')
plots['0.025'].GetYaxis().SetTitle('95% C.L. limit on #sigma(pp#rightarrowHH)#times#font[52]{B}(HH#rightarrowb#bar{b}#gamma#gamma) [fb]')
plots['0.025'].SetMaximum(12)

plots['0.160'].Draw("3Z same")
plots['0.500'].Draw("LZ same")
#plots['-1'].Draw("CZ same")
nonresXSEC.Draw("3ZL same")
nonresXSEC_2.Draw("3ZL same")

leg = TLegend(0.12, 0.45, 0.45, 0.89)
leg.SetFillColorAlpha(kWhite, 0.8)
leg.SetBorderSize(0)
leg.SetHeader("#font[61]{pp#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}")#{#kappa_{t} = 1, c_{g} = c_{2g} = c_{2} = 0}")
leg.AddEntry(plots['-1'], 'Observed', 'l')
leg.AddEntry(plots['0.500'], 'Expected', 'l')
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
latex.DrawLatex(0.89,.89,"#kappa_{t} = 1, c_{g} = c_{2g} = c_{2} = 0")

c.SaveAs("klscan.pdf")

outfile.cd()
for qt in quantiles:
  plots[qt].Write()

nonresXSEC.Write()
nonresXSEC_2.Write()
outfile.Close()
'''
