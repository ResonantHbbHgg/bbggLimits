from ROOT import *
import os,sys, math
#import RooStats

import argparse
parser =  argparse.ArgumentParser(description='Residual checkerr')
parser.add_argument('-t',  choices=['signal', 'background'], dest='type', required=True, default=None, help = "Background or signal.")
parser.add_argument('-w', '--workspace', dest="workspace", type=str, default=False, required=True, help="Workspace file location")
parser.add_argument('-c', '--category', dest='cat', type=int, default=0, required=True, help="category")
parser.add_argument('--alphaMjj', dest='alphaMjj', type=float, default=10, help='linear correlation factor')
parser.add_argument('--alphaMgg', dest='alphaMgg', type=float, default=1, help='linear correlation factor')
parser.add_argument('--nToys', dest='nToys', type=int, default=5, help='number of toy datasets')
opt = parser.parse_args()

gSystem.Load("libHiggsAnalysisCombinedLimit.so")

print 'Calculating', opt.type, 'residuals'

wfile = TFile(opt.workspace)
w = wfile.Get("w_all")

mgg = w.var("mgg")
mjj = w.var("mjj")

mgg_pdf_p0 = w.var('CMS_hhbbgg_13TeV_mgg_bkg_slope1_cat'+str(opt.cat)).getVal()
mgg_pdf_p1 = w.var('CMS_hhbbgg_13TeV_mgg_bkg_slope2_cat'+str(opt.cat)).getVal()
mgg_pdf_p2 = w.var('CMS_hhbbgg_13TeV_mgg_bkg_slope3_cat'+str(opt.cat)).getVal()

mjj_pdf_p0 = w.var('CMS_hhbbgg_13TeV_mjj_bkg_slope1_cat'+str(opt.cat)).getVal()
mjj_pdf_p1 = w.var('CMS_hhbbgg_13TeV_mjj_bkg_slope2_cat'+str(opt.cat)).getVal()
mjj_pdf_p2 = w.var('CMS_hhbbgg_13TeV_mjj_bkg_slope3_cat'+str(opt.cat)).getVal()

nW = RooWorkspace('nW')

nW.factory('mjj[80,200]')
nW.factory('mgg[100,180]')

#nW.factory("RooRealVar:mjj_c('mjj + alphaMjj*mgg', {mjj,mgg,alphaMjj["+str(opt.alphaMjj)+", -10, 10]})")
#nW.factory("RooRealVar:mgg_c('mgg + alphaMgg*mjj', {mjj,mgg,alphaMgg["+str(opt.alphaMgg)+", -10, 10]})")

nW.factory("RooFormulaVar:m_b2_p0_mgg('b2_p0_mgg*b2_p0_mgg + alphaMjj*mjj', {b2_p0_mgg["+str(mgg_pdf_p0)+", -50, 50], mjj, alphaMjj["+str(opt.alphaMjj)+", -500, 500]})")
nW.factory("RooFormulaVar:m_b2_p1_mgg('b2_p1_mgg*b2_p1_mgg', {b2_p1_mgg["+str(mgg_pdf_p1)+", -50, 50]})")
nW.factory("RooFormulaVar:m_b2_p2_mgg('b2_p2_mgg*b2_p2_mgg', {b2_p2_mgg["+str(mgg_pdf_p2)+", -50, 50]})")
nW.factory("RooBernstein:Bern2_mgg(mgg, {m_b2_p0_mgg, m_b2_p1_mgg, m_b2_p2_mgg})")

#nW.factory("RooFormulaVar:m_b2_p0_mjj('b2_p0_mjj*b2_p0_mjj + alphaMjj*mgg', {b2_p0_mjj["+str(mjj_pdf_p0)+"], mgg, alphaMjj["+str(opt.alphaMgg)+", -10, 10]})")
nW.factory("RooFormulaVar:m_b2_p0_mjj('b2_p0_mjj*b2_p0_mjj', {b2_p0_mjj["+str(mjj_pdf_p0)+", -50, 50]})")
nW.factory("RooFormulaVar:m_b2_p1_mjj('b2_p1_mjj*b2_p1_mjj', {b2_p1_mjj["+str(mjj_pdf_p1)+", -50, 50]})")
nW.factory("RooFormulaVar:m_b2_p2_mjj('b2_p2_mjj*b2_p2_mjj', {b2_p2_mjj["+str(mjj_pdf_p2)+", -50, 50]})")
#nW.factory("RooBernstein:Bern2_mjj(mjj_c, {m_b2_p0_mjj, m_b2_p1_mjj, m_b2_p2_mjj})")
nW.factory("RooBernstein:Bern2_mjj(mjj, {m_b2_p0_mjj, m_b2_p1_mjj, m_b2_p2_mjj})")

nW.factory("PROD:model( Bern2_mgg, Bern2_mjj)")

nW.defineSet("poi", "alphaMjj")
nW.defineSet("obs", "mjj,mgg")

modelCfg = RooStats.ModelConfig("2d")
modelCfg.SetWorkspace(nW)
modelCfg.SetPdf(nW.pdf('model'))
modelCfg.SetParametersOfInterest(nW.set('poi'))
modelCfg.SetObservables(nW.set('obs'))

grs = []
grs_y = []

for alp in range(0, opt.nToys):
  gr = TGraphErrors()
  gr_y = []
  thisAlpha = opt.alphaMjj
  gr.SetName("alphameasurements_"+str(thisAlpha)+"_toy"+str(alp))
  print '####################################################################################'
  print 'Doing alpha', opt.alphaMjj
  print '####################################################################################'
  
  for nt in range(1, 15):
    print '####################################################################################'
    print 'Generating toy #',nt
    print '####################################################################################'
    observables = RooArgSet(nW.set('obs'))
    nEvs = 10*nt
    expected = int(nEvs)
    myModel = nW.pdf('model')
    nW.var('alphaMjj').setVal(float(thisAlpha))
    nW.var('b2_p0_mgg').setVal(mgg_pdf_p0)
    nW.var('b2_p1_mgg').setVal(mgg_pdf_p1)
    nW.var('b2_p2_mgg').setVal(mgg_pdf_p2)
    nW.var('b2_p0_mjj').setVal(mjj_pdf_p0)
    nW.var('b2_p1_mjj').setVal(mjj_pdf_p1)
    nW.var('b2_p2_mjj').setVal(mjj_pdf_p2)
    toy = myModel.generate(observables, expected)
    print '####################################################################################'
    print 'Fitting toy #',nt
    print '####################################################################################'
    nW.var('alphaMjj').setVal(float(0.00000001))
    myModel.fitTo(toy)
    alphaVal = nW.var('alphaMjj').getVal()
    alphaErr = nW.var('alphaMjj').getError()
    gr.SetPoint(nt-1, nEvs, alphaVal)
    if alphaErr < 1:
      gr.SetPointError(nt-1, 0, alphaErr)
    else:
      gr.SetPointError(nt-1, 0, 0)
    gr_y.append(alphaVal)
  grs.append(gr)
  grs_y.append(gr_y)
#  break
  
gr_av = TGraph()
gr_av.SetName("alpha_avg")
for i in range(1,15):
  tot = 0
  for j in grs_y:
    tot += j[i-1]
  avg = float(tot)/float(len(grs_y))
  gr_av.SetPoint(i-1, i*10, avg)
  

#  nW.var('alphaMjj').setRange(alphaVal - alphaErr*10, alphaVal + alphaErr*10)
#  nll = RooNLLVar('nll', 'nll', myModel, toy)
#  pll = RooProfileLL('pll', 'pll', nll, RooArgSet(nW.var('alphaMjj')))
#  mframe = nW.var('alphaMjj').frame(alphaVal - alphaErr*10, alphaVal + alphaErr*10)
#  pll.plotOn(mframe

fout = TFile("out_corr"+str(opt.alphaMjj)+"_2.root", "RECREATE")
fout.cd()
c = TCanvas('c', 'c', 800, 600)
for i,g in enumerate(grs):
  if i == 0:
    g.Draw("APLE*")
    c.Update()
    g.GetXaxis().SetTitle("#Events in toy MC")
    g.GetYaxis().SetTitle("#alpha [1/GeV]")
    g.GetYaxis().SetRangeUser(-1,1)
    g.GetYaxis().SetLimits(-1,1)
  else:
    g.Draw("PLE*")
  g.SetLineColor(i+1)
  g.Write()
#  mframe.Draw()
c.SaveAs("pll_m25"+str(opt.alphaMjj)+"_2.pdf")
gr_av.Write()
fout.Close()
