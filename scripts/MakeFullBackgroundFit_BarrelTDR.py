from ROOT import *
from HiggsAnalysis.bbggLimits.NiceColors import *
import argparse, sys
from HiggsAnalysis.bbggLimits.MyCMSStyle import *

#python scripts/MakeFullBackgroundFit.py -i TDR_PLOT_14TeV/MaxLikelihoodFitResult.root -o FullBkgPlot_HM --signalNormalization 1 1 --signalFactor 1 1 --text "m_{hh} > 350 GeV" --unblind 

def MakeFullBackgroundPdf(bkg_pdf, bkg_norm, hig_pdfs, hig_norms):
  if len(hig_pdfs) == 0: 
    print 'No HIG background found'
    return bkg_pdf
  argPdfs = RooArgList()
  argNorms = RooArgList()
  argPdfs.add(bkg_pdf)
  argNorms.add(bkg_norm)
  if len(hig_pdfs) != len(hig_norms):
    print "list of higgs pdfs has different size wrt normalizations!"
    return None
  for hh in range(0, len(hig_pdfs)):
#    hig_pdfs[hh].Print()
    argPdfs.add(hig_pdfs[hh])
    argNorms.add(hig_norms[hh])
#  argNorms.Print()
#  argPdfs.Print()
  if argNorms.getSize() != argPdfs.getSize():
    print 'ArgNorms and ArgPdfs have different sizes!'
    return None
  totPdf = RooAddPdf('totBkg', 'Nonresonant + single H background', argPdfs, argNorms)
  return totPdf

gSystem.Load("libHiggsAnalysisCombinedLimit.so")
gROOT.SetBatch(kTRUE)

parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument('-i', '--inputFile', dest="dname", type=str, default=None, required=True)
parser.add_argument('-o', '--outFile', dest="outf", type=str, default=".")
parser.add_argument('-L', '--lumi', dest='lumi', type=str, default='3.0')
parser.add_argument('-Sn', '--signalNormalization', dest='snorm', type=float, default=[10.0,10.0], nargs='+')
parser.add_argument('-Sf', '--signalFactor', dest='fsignal', type=float, default=[10.0,10.0], nargs='+')
parser.add_argument('-H', '--addHiggs', dest='addh', action='store_true', default=False)
parser.add_argument('-Hl', '--higgsList', dest='hlist', type=str, default=None, nargs='+',
                           choices=['ggh', 'vbf', 'tth', 'vh', 'bbh'] )
parser.add_argument('-t', '--text', dest='text', type=str, default='')
parser.add_argument('-u', '--unblind', dest='unblind', default=False, action='store_true')

opt = parser.parse_args()

tfile = TFile(opt.dname, "READ")
ofile = TFile(opt.outf+".root", "RECREATE")

w_all = tfile.Get("MaxLikelihoodFitResult")

w_all.Print()

cats = [2,3]
#if 'High' in opt.dname or 'High' in opt.text: cats = [2,3]

Higgses = opt.hlist
if Higgses == None: Higgses = ['ggh', 'vbf', 'tth', 'vh', 'bbh']

dims = ['mjj', 'mgg']
bins = [16, 50]
binlow = [80, 100]
binhigh = [160, 150]
xtitle = ['M(jj) [GeV]', 'M(#gamma#gamma) [GeV]']
ytitle = ['Events/(5 GeV)', 'Events/(1 GeV)']
for cc in cats:
 for iobs,obs in enumerate(dims):

  intc = cc
  if intc > 1: intc = cc - 2

  var = w_all.var(obs)
  data_cat = w_all.obj('CMS_channel')

  print '=============================== ',data_cat.Print()

  catbin = 'cat'
  data_cat.setRange("catcut",catbin+str(cc))
#  var.Print()

  sig_pdf_name = obs+'Sig_cat'+str(intc)+'_CMS_sig_cat'+str(cc-2)
  print sig_pdf_name
  sig_pdf = w_all.pdf(sig_pdf_name)
  sig_pdf.Print()
  signormName = "n_exp_bincat"
  sig_norm = RooRealVar('sig_norm', 'sig norm', w_all.obj(signormName+str(cc)+'_proc_Sig').getVal())
  print  'sigNorm = ',signormName+str(cc)+'_proc_Sig ',sig_norm.getVal()

  bkg_pdf_name = obs+'BkgTmpBer1_cat'+str(intc)#+'_CMS_Bkg_cat'+str(cc)
  print "==============", bkg_pdf_name
  bkg_pdf = w_all.pdf(bkg_pdf_name)
#  bkg_pdf.Print()
  normName = 'n_exp_final_bincat'
  print normName+str(cc)+'_proc_Bkg'
  if cc == 0 or cc == 1: normName = 'n_exp_final_bincat'
  bkg_norm = RooRealVar('bkg_norm', 'nonres bkg norm', w_all.obj(normName+str(cc)+'_proc_Bkg').getVal())

  data2d = w_all.data("model_sData")
  print "before cut ", data2d.Print()
  
#  data = data2d.reduce(RooArgSet(RooFit.CutRange('catcut')))
  data = data2d.reduce(RooFit.CutRange('catcut'))
  data.Print()
#  sys.exit()

  hig_pdfs = []
  hig_norms = []
  totHiggs = 0
  if 1==1:
    for hh in Higgses:
      hig_pdf_name = obs+'Hig_'+hh+'_cat'+str(intc)+'_CMS_hig_'+hh+'_cat'+str(cc)
      print hig_pdf_name
      hig_pdf = w_all.pdf(hig_pdf_name)
      hig_pdfs.append(hig_pdf)
      normNameHIG = 'n_exp_bincat'
#      if cc == 2 or cc == 3: normNameHIG = 'n_exp_bincat'
      print normNameHIG+str(cc)+'_proc_'+hh
      norm = RooRealVar(hh+'_norm', hh+' bkg norm', w_all.obj(normNameHIG+str(cc)+'_proc_sh_'+hh).getVal() )
      hig_norms.append(norm)
      print hig_pdf_name, norm.getVal()
      totHiggs += w_all.obj(normNameHIG+str(cc)+'_proc_sh_'+hh).getVal()


  argPdfs = RooArgList()
  argNorms = RooArgList()
  for hh in range(0, len(hig_pdfs)):
#    hig_pdfs[hh].Print()
    argPdfs.add(hig_pdfs[hh])
    argNorms.add(hig_norms[hh])

  totHiggsPdf = RooAddPdf('totHiggsBkg', 'single H background', argPdfs, argNorms)

  hig_pdfs.append(sig_pdf)
  hig_norms.append(sig_norm)

  totBkg = MakeFullBackgroundPdf(bkg_pdf, bkg_norm, hig_pdfs, hig_norms)
  totBkgNorm = totHiggs + bkg_norm.getVal()



#  print totBkg
  print 'Total nonres:', bkg_norm.getVal(), 'total higgs:', totHiggs, totBkgNorm, 'total Sig:'
#  sys.exit()

  binning = bins[iobs]

  frame = var.frame(RooFit.Title(" "),RooFit.Bins(binning), RooFit.Range(binlow[iobs],binhigh[iobs]))
  dataind = 0
#  if not opt.unblind:
#    dataind = 2
#    blindedRegions = {}
#    blindedRegions['mgg'] = [100, 120, 130, 180]
#    blindedRegions['mjj'] = [70, 80, 140, 190]
#    var.removeRange("unblindReg_1")
#    var.removeRange("unblindReg_2")
#    var.setRange("unblindReg_1",blindedRegions[var.GetName()][0],blindedRegions[var.GetName()][1])
#    var.setRange("unblindReg_2",blindedRegions[var.GetName()][2],blindedRegions[var.GetName()][3])
#    data.plotOn(frame,RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0), RooFit.CutRange("unblindReg_1"))
#    data.plotOn(frame,RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0), RooFit.CutRange("unblindReg_2"))
#    data.plotOn(frame,RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0), RooFit.Invisible())
#    var.removeRange("unblindReg_1")
#    var.removeRange("unblindReg_2")
#  else:
#    data.plotOn(frame,RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0), RooFit.Range(binlow[iobs], binhigh[iobs]))


  bkg_pdf.plotOn(frame,RooFit.LineColor(cNiceGreenDark), RooFit.LineStyle(kDashed), RooFit.Precision(1E-5), RooFit.Normalization(bkg_norm.getVal(), RooAbsReal.NumEvent))
  totBkg.plotOn(frame,RooFit.LineColor(cNiceMidnight),RooFit.Precision(1E-5), RooFit.Normalization(totBkgNorm, RooAbsReal.NumEvent))
  if 'mjj' in obs:
    sig_pdf.plotOn(frame,RooFit.LineColor(cNiceRed), RooFit.Precision(1E-5), RooFit.Normalization(sig_norm.getVal()*10,RooAbsReal.NumEvent))
    totHiggsPdf.plotOn(frame,RooFit.LineColor(cNiceBlue),RooFit.Precision(1E-5),  RooFit.Normalization(totHiggs*10, RooAbsReal.NumEvent))
  else:
    sig_pdf.plotOn(frame,RooFit.LineColor(cNiceRed), RooFit.Precision(1E-5), RooFit.Normalization(sig_norm.getVal(),RooAbsReal.NumEvent))
    totHiggsPdf.plotOn(frame,RooFit.LineColor(cNiceBlue),RooFit.Precision(1E-5), RooFit.Normalization(totHiggs, RooAbsReal.NumEvent))

  RooRandom.randomGenerator().SetSeed(111)
  binnedData = totBkg.generate( RooArgSet(var), totBkgNorm)
  binnedData.plotOn(frame,RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0), RooFit.Range(binlow[iobs], binhigh[iobs]))
##  binnedData = data.binnedClone() 
  



#  if not opt.unblind:
#    dataind = 2
#    blindedRegions = {}
#    blindedRegions['mgg'] = [100, 120, 130, 180]
#    blindedRegions['mjj'] = [70, 80, 140, 190]
#    var.removeRange("unblindReg_1")
#    var.removeRange("unblindReg_2")
#    var.setRange("unblindReg_1",blindedRegions[var.GetName()][0],blindedRegions[var.GetName()][1])
#    var.setRange("unblindReg_2",blindedRegions[var.GetName()][2],blindedRegions[var.GetName()][3])
#    data.plotOn(frame,RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0), RooFit.CutRange("unblindReg_1"))
#    data.plotOn(frame,RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0), RooFit.CutRange("unblindReg_2"))
#    data.plotOn(frame,RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0), RooFit.Invisible())
#    var.removeRange("unblindReg_1")
#    var.removeRange("unblindReg_2")
#  else:
#    data.plotOn(frame,RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0))

  datahist = frame.getObject(4)
  bkghist = frame.getObject(0)
  totbkgh = frame.getObject(1)
  sigh = frame.getObject(2)
  singlh = frame.getObject(3)

  leg = TLegend(0.5, 0.60, 0.89, 0.89)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetTextFont(43)
  leg.SetTextSize(18)
#  leg.SetNColumns(3)
  leg.AddEntry(datahist, 'Toy data', 'pe')
  leg.AddEntry(totbkgh, 'Total line shape', 'l')
  leg.AddEntry(bkghist, 'Nonresonant background', 'l')

  if 'mjj' in obs:

    leg.AddEntry(singlh, 'Single Higgs (x10)', 'l')
#  sigText = 'SM HH Signal (x20)'
#  if int(intc) == 1:
    sigText = 'SM HH Signal (x10)'
    leg.AddEntry(sigh, sigText, 'l')

  else:

    leg.AddEntry(singlh, 'Single Higgs', 'l')
    sigText = 'SM HH Signal'
    leg.AddEntry(sigh, sigText, 'l')


  SetGeneralStyle()
  c = TCanvas("c", "c", 800, 600)
  frame.Draw()
  frame.GetXaxis().SetTitle(xtitle[iobs])
  frame.GetYaxis().SetTitle(ytitle[iobs])
  frame.SetMaximum(frame.GetMaximum()*1.6)
  frame.SetMinimum(0.0001)
  leg.Draw('same')
  c.Update()
  SetPadStyle(c)
  c.Update()
  SetAxisTextSizes(frame)
  c.Update()

  topy = 0.91
  stepy = 0.08
  tlatex = TLatex()
  tlatex.SetNDC()
  tlatex.SetTextAngle(0)
  tlatex.SetTextColor(kBlack)
  tlatex.SetTextFont(63)
  tlatex.SetTextAlign(11)
  tlatex.SetTextSize(25)
#  tlatex.DrawLatex(0.11, topy, "CMS")
  tlatex.SetTextFont(53)
#  tlatex.DrawLatex(0.18, topy, "Preliminary")
  tlatex.SetTextFont(43)
  tlatex.SetTextSize(20)
  tlatex.SetTextAlign(31)
#  tlatex.DrawLatex(0.9, topy, "L = " + str(opt.lumi) + " fb^{-1} (13 TeV)")
  tlatex.SetTextAlign(11)
  tlatex.SetTextSize(25)
  Cat = "High purity category"
  if int(intc) == 1:
    Cat = "Medium purity category"
  if "|" in opt.text:
    an = opt.text.split("|")
#               tlatex.SetTextFont(63)
    tlatex.DrawLatex(0.14, topy-stepy*1, an[0])
#               tlatex.SetTextFont(43)
    tlatex.DrawLatex(0.14, topy-stepy*2, an[1])
    tlatex.DrawLatex(0.14, topy-stepy*3, Cat)
  else:
#               tlatex.SetTextFont(63)
    tlatex.DrawLatex(0.14, topy-stepy*1, opt.text)
#               tlatex.SetTextFont(43)
    tlatex.DrawLatex(0.14, topy-stepy*2, Cat)

  DrawCMSLabels(c, '3000')
  c.SaveAs(opt.outf+str(cc) + obs+".pdf")
  c.SaveAs(opt.outf+str(cc) + obs+".png")

ofile.Close()
