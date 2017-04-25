from ROOT import *
import os,sys, math
#0.00194105866831 0.00194105866831 0.00101209105924
def PoissonianError(Obs, Exp):
  alpha = 1 - 0.6827
  if Obs < 1E-5: Obs = 0
  if Exp < 1E-5: Exp = 0
  if Obs == 0: L = 0
  else: L = Math.gamma_quantile(alpha/2., Obs, 1.)
  U = Math.gamma_quantile_c(alpha/2., Obs+1, 1.)
  diff = Obs - Exp
  if diff > 0: err = Obs-L
  else: err = U-Obs
#  print err, Obs, Exp, L, U, Obs-Exp, (Obs-Exp)/err
  return err
  
gROOT.SetBatch(True)

import argparse
parser =  argparse.ArgumentParser(description='Residual checkerr')
parser.add_argument('-o', '--observed', dest='observed', type=str, required=True, default=None, help = "Background or signal.")
parser.add_argument('-e', '--expected', dest="expected", type=str, default=None, required=True, help="Workspace file location")
parser.add_argument('-n', '--normalization', dest='normalization', type=float, default=10)
parser.add_argument('--output', '--outputLocation', dest='output', type=str, default='')
parser.add_argument('-l', '--label', dest='label', type=str, required=True)
opt = parser.parse_args()

obs_ar = opt.observed.split(':')
exp_ar = opt.expected.split(':')

if len(obs_ar) < 2 or len(exp_ar) < 2:
  print "-o file:histo -e file:histo"
  sys.exit()

obs_f = TFile(obs_ar[0])
exp_f = TFile(exp_ar[0])

obs_h = obs_f.Get(obs_ar[1])
exp_h = exp_f.Get(exp_ar[1])

nbinsX = obs_h.GetXaxis().GetNbins()
nbinsY = obs_h.GetYaxis().GetNbins()

obs_h.Scale(opt.normalization)
exp_h.Scale(opt.normalization)

residual = TH2F('residual', ';mgg [GeV];mjj [GeV]', nbinsX, 100, 180, nbinsY, 80,200)
residual_u = TH2F('residual_u', ';mgg [GeV];mjj [GeV]', nbinsX, 100, 180, nbinsY, 80,200)

for X in range(1, nbinsX+1):
 for Y in range(1, nbinsY+1):
  print X, Y
  obs_b = obs_h.GetBinContent(X, Y)
  exp_b = exp_h.GetBinContent(X, Y)
  diff = (obs_b - exp_b)
  err = PoissonianError(obs_b, exp_b)
  residual_u.SetBinContent(X,Y, diff/err)
  myres = diff/err
  if myres > 0.15: myres = 0.1499999
  if myres < -0.15: myres = -0.149999
  residual.SetBinContent(X,Y, myres)
  

outFile = TFile("ResidualsTH2F_"+str(opt.label)+".root", "RECREATE")
outFile.cd()
obs_h.Write()
exp_h.Write()
residual.Write()
residual_u.Write()
outFile.Close()

gStyle.SetOptStat(0)

c = TCanvas('c', 'c', 800, 600)
obs_h.Draw("COLZ")
c.SaveAs(opt.output+'res_th2F_obs_'+str(opt.label)+'.pdf')
c.SaveAs(opt.output+'res_th2F_obs_'+str(opt.label)+'.png')

c1 = TCanvas('c1', 'c1', 800, 600)
exp_h.Draw("COLZ")
c1.SaveAs(opt.output+'res_th2F_exp_'+str(opt.label)+'.pdf')
c1.SaveAs(opt.output+'res_th2F_exp_'+str(opt.label)+'.png')

c2 = TCanvas('c2', 'c2', 800, 600)
residual.GetZaxis().SetRangeUser(-0.15,0.15)
residual.Draw("COLZ")
c2.SaveAs(opt.output+'res_th2F_res_'+str(opt.label)+'.pdf')
c2.SaveAs(opt.output+'res_th2F_res_'+str(opt.label)+'.png')
