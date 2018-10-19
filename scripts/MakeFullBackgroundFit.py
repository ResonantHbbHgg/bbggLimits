from ROOT import *
from HiggsAnalysis.bbggLimits.NiceColors import *
import argparse, sys
from HiggsAnalysis.bbggLimits.MyCMSStyle import *
import numpy as n

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


def MakeFullModelPdf(sig_pdf, sig_norm, bkg_pdf, bkg_norm, hig_pdfs, hig_norms):
  if len(hig_pdfs) == 0:
    print 'No HIG background found'
    return bkg_pdf
  argPdfs = RooArgList()
  argNorms = RooArgList()
  argPdfs.add(bkg_pdf)
  argPdfs.add(sig_pdf)
  argNorms.add(bkg_norm)
  argNorms.add(sig_norm)
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
  totModPdf = RooAddPdf('totModPdf', 'Signal + Nonresonant background + single H background', argPdfs, argNorms)
  return totModPdf

gSystem.Load("libHiggsAnalysisCombinedLimit.so")
gROOT.SetBatch(kTRUE)

parser =  argparse.ArgumentParser(description='Background fit plot maker')
parser.add_argument('-i', '--inputFile', dest="dname", type=str, default=None, required=True)
parser.add_argument('-o', '--outFile', dest="outf", type=str, default=".")
parser.add_argument('-L', '--lumi', dest='lumi', type=str, default='3000')
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

dims = ['mjj', 'mgg']
bins = [25, 25]
minval = [93.2, 105]
maxval = [190, 145]

xtitle = ['m_{jj} [GeV]', 'm_{#gamma#gamma} [GeV]']
ytitle = ['Events/(5 GeV)', 'Events/(1 GeV)']
yLimits = {'mgg': [60, 700, 60, 400], 'mjj': [110, 1000, 110, 700]}
#yLimits = {'mgg': [0.03, 0.1, 0.03, 0.1], 'mjj': [0.03, 0.1, 0.03, 0.1]}
#yLimits = {'mgg': [6, 500, 6, 700], 'mjj': [200, 1200, 200, 2000]}

for cc in cats:
 for iobs,obs in enumerate(dims):

  intc = cc
  if intc > 1: intc = cc - 2

  var = w_all.var(obs)
  data_cat = w_all.obj('CMS_channel')
  catbin = 'ch1_cat'
  if cc == 0 or cc == 1: catbin = 'ch2_cat'
  data_cat.setRange("catcut",catbin+str(cc))
#  var.Print()
  var.Print()

  sig_pdf_name = obs+'Sig_cat'+str(intc)+'_CMS_sig_cat'+str(cc)
  sig_pdf = w_all.pdf(sig_pdf_name)
#  print sig_pdf_name
  print sig_pdf_name
  sig_pdf.Print()
  normName = 'n_exp_binch1_cat'
  if cc == 0 or cc == 1: normName = 'n_exp_binch2_cat'
  ####sig_norm = w_all.obj(normName+str(cc)+'_proc_Sig').getVal()
  sig_norm = RooRealVar('sig_norm','signal norm', w_all.obj(normName+str(cc)+'_proc_Sig').getVal())
  #print "@@@ sig_norm: %.1f" % sig_norm
  print "@@@ sig_norm: %.1f" % sig_norm.getVal()

  bkg_pdf_name = obs+'BkgTmpExp1_cat'+str(intc)+'_CMS_Bkg_cat'+str(cc)
  bkg_pdf = w_all.pdf(bkg_pdf_name)
#  bkg_pdf.Print()
  bkg_pdf.Print()
  normName = 'n_exp_final_binch1_cat'
  if cc == 0 or cc == 1: normName = 'n_exp_final_binch2_cat'
  bkg_norm = RooRealVar('bkg_norm', 'nonres bkg norm', w_all.obj(normName+str(cc)+'_proc_Bkg').getVal())


  data2d = w_all.data("model_sData")
#  data2d = w_all.data("data_obs")
  data2d.Print()
  
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

  model = MakeFullModelPdf(sig_pdf, sig_norm, bkg_pdf, bkg_norm, hig_pdfs, hig_norms)

#  print totBkg
  print 'Total nonres:', bkg_norm.getVal(), 'total higgs:', totHiggs, totBkgNorm
#  sys.exit()

  binning = bins[iobs]

  frame = var.frame(RooFit.Title(" "),RooFit.Bins(binning),RooFit.Range(minval[iobs],maxval[iobs]))
  dataind = 0
  if not opt.unblind:
    dataind = 2
    blindedRegions = {}
    blindedRegions['mgg'] = [100, 120, 130, 180]
    blindedRegions['mjj'] = [70, 80, 140, 190]
#    var.removeRange("unblindReg_1")
#    var.removeRange("unblindReg_2")
    var.setRange("unblindReg_1",blindedRegions[var.GetName()][0],blindedRegions[var.GetName()][1])
    var.setRange("unblindReg_2",blindedRegions[var.GetName()][2],blindedRegions[var.GetName()][3])
    data.plotOn(frame,RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0), RooFit.CutRange("unblindReg_1"))
    data.plotOn(frame,RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0), RooFit.CutRange("unblindReg_2"))
    data.plotOn(frame,RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0), RooFit.Invisible())
    var.removeRange("unblindReg_1")
    var.removeRange("unblindReg_2")
  else:
#    data.plotOn(frame,RooFit.DataError(RooAbsData.Poisson),RooFit.XErrorSize(0))
    data.plotOn(frame,RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0))
    
  
  bkg_pdf.plotOn(frame,RooFit.LineColor(cNiceGreenDark), RooFit.LineStyle(kDashed), RooFit.Precision(1E-5), RooFit.Normalization(bkg_norm.getVal(), RooAbsReal.NumEvent))
  totBkg.plotOn(frame,RooFit.LineColor(cNiceBlueDark),RooFit.Precision(1E-5), RooFit.Normalization(totBkgNorm, RooAbsReal.NumEvent))

###  sig_pdf.plotOn(frame,RooFit.LineColor(cNiceRed), RooFit.Precision(1E-5), RooFit.Normalization(opt.snorm[intc]*opt.fsignal[intc],RooAbsReal.NumEvent))

  ##model.plotOn(frame,RooFit.LineColor(cNiceRed), RooFit.Precision(1E-5), RooFit.Normalization(sig_norm.getVal()*opt.fsignal[intc],totBkgNorm,RooAbsReal.NumEvent))
  model.plotOn(frame,RooFit.LineColor(cNiceRed), RooFit.Precision(1E-5))

  ##sig_pdf.plotOn(frame,RooFit.LineColor(cNiceRed), RooFit.LineStyle(kDashed), RooFit.Precision(1E-5), RooFit.Normalization((sig_norm.getVal())*opt.fsignal[intc],RooAbsReal.NumEvent))

  datahist = frame.getObject(0)

  rnd = TRandom(2)

  x=n.zeros(1, dtype=float)
  xE=n.zeros(1, dtype=float)
  y=n.zeros(1, dtype=float)
  yE=n.zeros(1, dtype=float)
  rndvalue=n.zeros(1, dtype=float)


  bkghist = frame.getObject(dataind+1)
  totbkgh = frame.getObject(dataind+2)
  modelhist = frame.getObject(dataind+3)
  ##sigh = frame.getObject(dataind+4)


  leg = TLegend(0.5, 0.55, 0.89, 0.89)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetTextFont(43)
  leg.SetTextSize(20)
#  leg.SetNColumns(3)
  leg.AddEntry(datahist, 'Pseudo-data', 'pe')
  leg.AddEntry(bkghist, 'Nonresonant background', 'l')
  leg.AddEntry(totbkgh, 'Full background', 'l')
  leg.AddEntry(modelhist, 'Signal + Full background', 'l')
  ####leg.AddEntry(sigh, 'SM HH signal', 'l')
#  sigText = 'SM HH Signal (x20)'
#  if int(intc) == 1:
  #####sigText = 'SM HH signal (x'+str(int(opt.fsignal[intc]))+')'
  #####leg.AddEntry(sigh, sigText, 'l')

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
#  tlatex.DrawLatex(0.11, topy, "CMS")
  tlatex.SetTextFont(53)
#  tlatex.DrawLatex(0.18, topy, "Preliminary")
  tlatex.SetTextFont(43)
  tlatex.SetTextSize(20)
  tlatex.SetTextAlign(31)
#  tlatex.DrawLatex(0.9, topy, "L = " + str(opt.lumi) + " fb^{-1} (13 TeV)")
  tlatex.SetTextAlign(11)
  tlatex.SetTextSize(25)
  Cat = "High-purity category"
  if int(intc) == 1:
    Cat = "Medium-purity category"
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

  N = datahist.GetN()
  alpha = 1 - 0.6827

  for i in range(0,N):
    datahist.GetPoint(i, x, y)
    rndvalue[0] = rnd.Poisson(y[0])

    if rndvalue[0] == 0 :
      L = 0
    else : 
      L = Math.gamma_quantile(alpha/2,rndvalue[0],1.)

    U =  Math.gamma_quantile_c(alpha/2,rndvalue[0]+1,1) 
   
    datahist.SetPointEYlow(i, rndvalue[0]-L)
    datahist.SetPointEYhigh(i, U-rndvalue[0])

#    Eup = sqrt()
#    Edo = 
    print "x = ",x," y = ",y, " rndvalue", rndvalue
    datahist.SetPoint(i,x,rndvalue)
#    datahist.SetPointEYhigh(i,x,Eup)
#    datahist.SetPointEYlow(i,x,Edo)


##  N = sigh.GetN()
##  totsig = 0
##  for i in range(1,N-1):
##    sigh.GetPoint(i,x,y)
##    xup = x[0] 
##    yup = y[0]
##    sigh.GetPoint(i-1,x,y)
##    totsig = totsig + yup*(xup-x[0])
#    print "x = ", x[0], " val = ", yup, "width = ", xup-x[0]

##  integral = sigh.Integral()

    
##  print "totsig = ", totsig, " integral = ", integral
##  print "----------------------------------------"


  c.SaveAs(opt.outf+str(cc) + obs+"_poiss.pdf")
  c.SaveAs(opt.outf+str(cc) + obs+"_poiss.png")

ofile.Close()
