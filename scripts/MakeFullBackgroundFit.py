from ROOT import *
from HiggsAnalysis.bbggLimits.NiceColors import *
import argparse, sys
from HiggsAnalysis.bbggLimits.MyCMSStyle import *

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
    argPdfs.add(hig_pdfs[hh])
    argNorms.add(hig_norms[hh])

  if argNorms.getSize() != argPdfs.getSize():
    print 'ArgNorms and ArgPdfs have different sizes!'
    return None
  totPdf = RooAddPdf('totBkg', 'Nonresonant + single H background', argPdfs, argNorms)
  return totPdf

gSystem.Load("libHiggsAnalysisCombinedLimit.so")
gROOT.SetBatch(kTRUE)

parser =  argparse.ArgumentParser(description='Background fit plot maker')
parser.add_argument('-i', '--inputFile', dest="dname", type=str, default=None, required=True)
parser.add_argument('-o', '--outFile', dest="outf", type=str, default=".")
parser.add_argument('-L', '--lumi', dest='lumi', type=str, default='35.9')
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

cats = [0,1]
if 'High' in opt.outf or 'High' in opt.text: cats = [2,3]

Higgses = opt.hlist
if Higgses == None: Higgses = ['ggh', 'vbf', 'tth', 'vh', 'bbh']

do1D = 0
dims = ['mjj','mgg']
bins = [24, 80]
xtitle = ['m_{jj} [GeV]', 'm_{#gamma#gamma} [GeV]']
ytitle = ['Events/(5 GeV)', 'Events/(1 GeV)']
#yLimits = {'mgg': [14, 90, 14, 60], 'mjj': [18, 220, 30, 150]}
yLimits = {'mgg': [90, 1100, 240, 1200], 'mjj': [180, 3200, 600, 3500]}

for cc in cats:
 for iobs,obs in enumerate(dims):
  if iobs==0 and do1D: continue 
  intc = cc
  if intc > 1: intc = cc - 2

  var = w_all.var(obs)
  data_cat = w_all.obj('CMS_channel')
  catbin = 'ch1_cat'
  if cc == 0 or cc == 1: catbin = 'ch2_cat'
  data_cat.setRange("catcut",catbin+str(cc))
  #  var.Print()

  if do1D:
    sig_pdf_name = 'shapeSig_Sig_'+catbin+str(cc)
  else:
    sig_pdf_name = obs+'Sig_cat'+str(intc)+'_CMS_sig_cat'+str(cc)

  print cc, sig_pdf_name
  sig_pdf = w_all.pdf(sig_pdf_name)
  sig_pdf.Print()

  if do1D:
    bkg_pdf_name = 'shapeBkg_Bkg_'+catbin+str(cc)
  else:
    bkg_pdf_name = obs+'BkgTmpBer1_cat'+str(intc)+'_CMS_Bkg_cat'+str(cc)
  bkg_pdf = w_all.pdf(bkg_pdf_name)
  #  bkg_pdf.Print()
  normName = 'n_exp_final_binch1_cat'
  if cc == 0 or cc == 1: normName = 'n_exp_final_binch2_cat'
  bkg_norm = RooRealVar('bkg_norm', 'nonres bkg norm', w_all.obj(normName+str(cc)+'_proc_Bkg').getVal())

  data2d = w_all.data("data_obs")
  data2d.Print()
  
  data = data2d.reduce(RooFit.CutRange('catcut'))
  data.Print()

  hig_pdfs = []
  hig_norms = []
  totHiggs = 0
  for hh in Higgses:
    if do1D:
      hig_pdf_name = 'shapeBkg_'+hh+'_'+catbin+str(cc)
    else:
      hig_pdf_name = obs+'Hig_'+hh+'_cat'+str(intc)+'_CMS_hig_'+hh+'_cat'+str(cc)
    hig_pdf = w_all.pdf(hig_pdf_name)
    hig_pdfs.append(hig_pdf)
    normNameHIG = 'n_exp_binch1_cat'
    if cc == 0 or cc == 1: normNameHIG = 'n_exp_binch2_cat'
    norm = RooRealVar(hh+'_norm', hh+' bkg norm', w_all.obj(normNameHIG+str(cc)+'_proc_'+hh).getVal() )
    hig_norms.append(norm)
    print hig_pdf_name, norm.getVal()
    totHiggs += w_all.obj(normNameHIG+str(cc)+'_proc_'+hh).getVal()

  totBkg = MakeFullBackgroundPdf(bkg_pdf, bkg_norm, hig_pdfs, hig_norms)
  totBkgNorm = totHiggs + bkg_norm.getVal()

  print 'Total nonres:', bkg_norm.getVal(), 'total higgs:', totHiggs, totBkgNorm


  binning = bins[iobs]

  frame = var.frame(RooFit.Title(" "),RooFit.Bins(binning))
  dataind = 0
  if not opt.unblind:
    dataind = 2
    blindedRegions = {}
    blindedRegions['mgg'] = [100, 120, 130, 180]
    blindedRegions['mjj'] = [70, 80, 140, 190]
    
    var.setRange("unblindReg_1",blindedRegions[var.GetName()][0],blindedRegions[var.GetName()][1])
    var.setRange("unblindReg_2",blindedRegions[var.GetName()][2],blindedRegions[var.GetName()][3])
    data.plotOn(frame,RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0), RooFit.CutRange("unblindReg_1"))
    data.plotOn(frame,RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0), RooFit.CutRange("unblindReg_2"))
    data.plotOn(frame,RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0), RooFit.Invisible())
    var.removeRange("unblindReg_1")
    var.removeRange("unblindReg_2")
  else:
    data.plotOn(frame,RooFit.DataError(RooAbsData.Poisson),RooFit.XErrorSize(0))

  bkg_pdf.plotOn(frame,RooFit.LineColor(cNiceGreenDark), RooFit.LineStyle(kDashed), RooFit.Precision(1E-5), RooFit.Normalization(bkg_norm.getVal(), RooAbsReal.NumEvent))
  totBkg.plotOn(frame,RooFit.LineColor(cNiceBlueDark),RooFit.Precision(1E-5), RooFit.Normalization(totBkgNorm, RooAbsReal.NumEvent))
  sig_pdf.plotOn(frame,RooFit.LineColor(cNiceRed), RooFit.Precision(1E-5), RooFit.Normalization(opt.snorm[intc]*opt.fsignal[intc],RooAbsReal.NumEvent))
  #  data.plotOn(frame,RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0))

  if not opt.unblind:
    dataind = 2
    blindedRegions = {}
    blindedRegions['mgg'] = [100, 120, 130, 180]
    blindedRegions['mjj'] = [70, 80, 140, 190]

    var.setRange("unblindReg_1",blindedRegions[var.GetName()][0],blindedRegions[var.GetName()][1])
    var.setRange("unblindReg_2",blindedRegions[var.GetName()][2],blindedRegions[var.GetName()][3])
    data.plotOn(frame,RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0), RooFit.CutRange("unblindReg_1"))
    data.plotOn(frame,RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0), RooFit.CutRange("unblindReg_2"))
    data.plotOn(frame,RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0), RooFit.Invisible())
    var.removeRange("unblindReg_1")
    var.removeRange("unblindReg_2")
  else:
    data.plotOn(frame,RooFit.DataError(RooAbsData.Poisson),RooFit.XErrorSize(0))

  datahist = frame.getObject(0)
  bkghist = frame.getObject(dataind+1)
  totbkgh = frame.getObject(dataind+2)
  sigh = frame.getObject(dataind+3)

  leg = TLegend(0.5, 0.55, 0.89, 0.89)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetTextFont(43)
  leg.SetTextSize(20)

  leg.AddEntry(datahist, 'Data', 'pe')
  leg.AddEntry(totbkgh, 'Full background model', 'l')
  leg.AddEntry(bkghist, 'Nonresonant background', 'l')

  sigText = 'SM HH signal (x'+str(int(opt.fsignal[intc]))+')'
  leg.AddEntry(sigh, sigText, 'l')

  SetGeneralStyle()
  c = TCanvas("c", "c", 800, 600)
  frame.Draw()
  frame.GetXaxis().SetTitle(xtitle[iobs])
  frame.GetYaxis().SetTitle(ytitle[iobs])
  frame.SetMaximum(yLimits[obs][cc])
  frame.SetMinimum(0.000)
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

  tlatex.SetTextFont(53)

  tlatex.SetTextFont(43)
  tlatex.SetTextSize(20)
  tlatex.SetTextAlign(31)

  tlatex.SetTextAlign(11)
  tlatex.SetTextSize(25)
  Cat = "High-purity category"
  if int(intc) == 1:
    Cat = "Medium-purity category"
  if "|" in opt.text:
    an = opt.text.split("|")

    tlatex.DrawLatex(0.14, topy-stepy*1, an[0])
    tlatex.DrawLatex(0.14, topy-stepy*2, an[1])
    tlatex.DrawLatex(0.14, topy-stepy*3, Cat)
  else:

    tlatex.DrawLatex(0.14, topy-stepy*1, opt.text)
    tlatex.DrawLatex(0.14, topy-stepy*2, Cat)

  DrawCMSLabels(c, '35.9')
  c.SaveAs(opt.outf+str(cc) + obs+".pdf")
  c.SaveAs(opt.outf+str(cc) + obs+".png")

ofile.Close()
