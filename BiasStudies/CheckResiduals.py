from ROOT import *
import os,sys, math
#0.00194105866831 0.00194105866831 0.00101209105924
def PoissonianError(Exp, Obs):
  alpha = 1 - 0.6827
  if Obs < 1E-5: Obs = 0
  if Exp < 1E-5: Exp = 0
  if Obs == 0: L = 0
  else: L = Math.gamma_quantile(alpha/2., Obs, 1.)
  U = Math.gamma_quantile_c(alpha/2., Obs+1, 1.)
  diff = Obs - Exp
  if diff > 0: err = Obs-L
  else: err = U-Obs
  print err, Obs, Exp, L, U, Obs-Exp, (Obs-Exp)/err
  return err
  

import argparse
parser =  argparse.ArgumentParser(description='Residual checkerr')
parser.add_argument('-t',  choices=['signal', 'background'], dest='type', required=True, default=None, help = "Background or signal.")
parser.add_argument('-w', '--workspace', dest="workspace", type=str, default=False, required=True, help="Workspace file location")
parser.add_argument('-c', '--category', dest='cat', type=int, default=0, required=True, help="category")
opt = parser.parse_args()

gSystem.Load("libHiggsAnalysisCombinedLimit.so")

print 'Calculating', opt.type, 'residuals'

mjjMin = 80
mjjMax = 200
mggMin = 100
mggMax = 180

wfile = TFile(opt.workspace)
w = wfile.Get("w_all")

mgg = w.var("mgg")
mjj = w.var("mjj")

data = 0
pdf_2d = 0
pdf_mgg = 0
pdf_mjj = 0
if opt.type == 'signal':
  data = w.data('Sig_cat'+str(opt.cat))
  pdf_2d = w.pdf('CMS_sig_cat'+str(opt.cat))
  pdf_mgg = w.pdf('mggSig_cat'+str(opt.cat))
  pdf_mjj = w.pdf('mjjSig_cat'+str(opt.cat))
else:
  data = w.data('data_obs_cat'+str(opt.cat))
  pdf_2d = w.pdf('BkgPdf_cat'+str(opt.cat))
  pdf_mgg = w.pdf('mggBkgTmpBer1_cat'+str(opt.cat))
  pdf_mjj = w.pdf('mjjBkgTmpBer1_cat'+str(opt.cat))

frame_mgg = mgg.frame(RooFit.Title(" "),RooFit.Bins(80))
frame_mjj = mjj.frame(RooFit.Title(" "),RooFit.Bins(40))

data.plotOn(frame_mgg,RooFit.DataError(RooAbsData.SumW2))
pdf_mgg.plotOn(frame_mgg)

data.plotOn(frame_mjj,RooFit.DataError(RooAbsData.SumW2))
pdf_mjj.plotOn(frame_mjj)

h_data_mgg = frame_mgg.getObject(0)
h_pdf_mgg = frame_mgg.getObject(1)

h_data_mjj = frame_mjj.getObject(0)
h_pdf_mjj = frame_mjj.getObject(1)

h_data_2d = data.createHistogram(mgg,mjj, 80, 40)
h_pdf_2d = pdf_2d.createHistogram("hh_model",mgg,RooFit.Binning(80),RooFit.YVar(mjj,RooFit.Binning(40)))  #pdf_2d.createHistogram(mgg,mjj, 80, 40)
h_pdf_2d_norm = h_pdf_2d.Integral()
h_pdf_2d.Scale(h_data_2d.Integral()/h_pdf_2d_norm)
h_pdf_2d_mgg = h_pdf_2d.ProjectionX('h_pdf_2d_mgg')
h_pdf_2d_mjj = h_pdf_2d.ProjectionY('h_pdf_2d_mjj')
h_data_2d_mgg = h_data_2d.ProjectionX('h_data_2d_mgg')
h_data_2d_mjj = h_data_2d.ProjectionY('h_data_2d_mjj')

print "Done making histograms, now calculating residuals..."

h_mjj_res = TH1F('h_mjj_res', '', 100, -1, 1)
h_mgg_res = TH1F('h_mgg_res', '', 100, -1, 1)
h_mgg_mjj_res = TH1F('h_mgg_mjj_res', '', 100, 0, 0.5)
h_2d_res = TH1F('h_2d_res', '', 100, -1, 1)
h_2d2_res = TH1F('h_2d2_res', '', 100, 0, 0.5)
h_res_diff = TH1F('h_res_diff', '', 100, -0.1, 0.1)

h_res_2d = TH2F('h_res_2d', '', 20, 115, 135, 40, 80, 200)

print h_data_2d_mgg.GetNbinsX(), h_data_2d_mjj.GetNbinsX()
print h_data_2d.GetNbinsX(), h_data_2d.GetNbinsY()
print h_pdf_2d.GetNbinsX(), h_pdf_2d.GetNbinsY()

for ig in range(1, 80):

  if h_pdf_2d.GetXaxis().GetBinCenter(ig) > 135 and opt.type=='signal': continue
  if h_pdf_2d.GetXaxis().GetBinCenter(ig) < 115 and opt.type=='signal': continue

#  pErr_2d_mgg = PoissonianError(h_data_2d_mgg.GetBinContent(ig), h_pdf_2d_mgg.GetBinContent(ig))
  pErr_2d_mgg = PoissonianError(h_pdf_2d_mgg.GetBinContent(ig), h_data_2d_mgg.GetBinContent(ig))
  v_2d_mgg = (h_data_2d_mgg.GetBinContent(ig) - h_pdf_2d_mgg.GetBinContent(ig))
#  g_res = (h_data_2d_mgg.GetBinContent(ig) - h_pdf_2d_mgg.GetBinContent(ig))/math.sqrt(h_pdf_2d_mgg.GetBinContent(ig))
  g_res = v_2d_mgg/pErr_2d_mgg #(h_data_2d_mgg.GetBinContent(ig) - h_pdf_2d_mgg.GetBinContent(ig))/pErr_2d_mgg
  h_mgg_res.Fill(g_res)
  for ij in range(1, 40):
#    pErr_2d_mjj = PoissonianError(h_data_2d_mjj.GetBinContent(ij), h_pdf_2d_mjj.GetBinContent(ij))
    pErr_2d_mjj = PoissonianError(h_pdf_2d_mjj.GetBinContent(ij), h_data_2d_mjj.GetBinContent(ij))
    v_2d_mjj = (h_data_2d_mjj.GetBinContent(ij) - h_pdf_2d_mjj.GetBinContent(ij))
    j_res = (h_data_2d_mjj.GetBinContent(ij) - h_pdf_2d_mjj.GetBinContent(ij))/math.sqrt(h_pdf_2d_mjj.GetBinContent(ij))
    j_res = v_2d_mjj/pErr_2d_mjj   #(h_data_2d_mjj.GetBinContent(ij) - h_pdf_2d_mjj.GetBinContent(ij))/pErr_2d_mjj
    h_mjj_res.Fill(j_res)
    h_mgg_mjj_res.Fill(abs(g_res*j_res))

#    pErr_2d = PoissonianError(h_data_2d.GetBinContent(ig,ij), h_pdf_2d.GetBinContent(ig,ij))
    pErr_2d = PoissonianError(h_pdf_2d.GetBinContent(ig,ij), h_data_2d.GetBinContent(ig,ij))
    v_2d = (h_data_2d.GetBinContent(ig,ij) - h_pdf_2d.GetBinContent(ig,ij))
    res_2d = (h_data_2d.GetBinContent(ig,ij) - h_pdf_2d.GetBinContent(ig,ij))/math.sqrt(h_pdf_2d.GetBinContent(ig,ij))
    res_2d = v_2d/pErr_2d #(h_data_2d.GetBinContent(ig,ij) - h_pdf_2d.GetBinContent(ig,ij))/pErr_2d
    h_2d_res.Fill(res_2d)
    h_2d2_res.Fill(res_2d*res_2d)

    h_res_diff.Fill(res_2d*res_2d - abs(g_res*j_res))
    h_res_diff.Fill( ((v_2d*v_2d) - (v_2d_mgg*v_2d_mjj))/( pErr_2d*pErr_2d + pErr_2d_mgg*pErr_2d_mjj ) )
    h_res_2d.Fill(h_pdf_2d.GetXaxis().GetBinCenter(ig), h_pdf_2d.GetYaxis().GetBinCenter(ij), res_2d )
    

outFile = TFile("Residuals_"+str(opt.type)+"_cat"+str(opt.cat)+"_2.root", "RECREATE")
outFile.cd()
h_data_2d.Write()
h_data_mgg.Write()
h_data_mjj.Write()
h_pdf_2d.Write()
h_pdf_mgg.Write()
h_pdf_mjj.Write()
h_pdf_2d_mgg.Write()
h_pdf_2d_mjj.Write()
h_data_2d_mgg.Write()
h_data_2d_mjj.Write()
h_mjj_res.Write()
h_mgg_res.Write()
h_mgg_mjj_res.Write()
h_2d_res.Write()
h_res_diff.Write()
h_2d2_res.Write()
h_res_2d.Write()
outFile.Close()
