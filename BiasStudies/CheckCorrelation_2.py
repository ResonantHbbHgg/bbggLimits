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

##Get 2D pdf from signal fit
pdf2d = w.pdf('BkgPdf_cat'+str(opt.cat))

##Create new workspace and pdf' = pdf2d + alpha*mjj*mgg
nW = RooWorkspace('nW')
nW.factory('mjj[80,200]')
nW.factory('mgg[100,180]')
getattr(nW,'import')(pdf2d, RooCmdArg())
nW.factory("EXPR::j1g1( 'mjj*mgg ', {mjj,mgg} )")
nW.factory("EXPR::j2g1( 'mjj*mjj*mgg ', {mjj,mgg} )")
nW.factory("EXPR::j1g2( 'mjj*mgg*mgg ', {mjj,mgg} )")
nW.factory("SUM:model(corr_beta[1, -10, 1000]*j2g1,BkgPdf_cat"+str(opt.cat)+")")
nW.factory('model_norm[1]')
nW.Print()

fout = TFile("corr_workspace_beta.root", "RECREATE")
nW.Write()
fout.Close()
