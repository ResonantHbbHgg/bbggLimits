import sys
import os, math
from ROOT import *
gROOT.SetBatch(True)
#import RooStats

import argparse
parser =  argparse.ArgumentParser(description='Print 2d pdfs')
parser.add_argument('-f', '--file', dest='File', required=True, default=None)
parser.add_argument('-w', '--workspace', dest='workspace', required=True, default=None)
parser.add_argument('-p', '--pdf', dest='pdf', required=True, default=None)
parser.add_argument('-l', '--label', dest='label', default='')
parser.add_argument('-n', '--name', dest='name', default='')

opt = parser.parse_args()

gSystem.Load("libHiggsAnalysisCombinedLimit.so")

tfile = TFile(opt.File)
w = tfile.Get(opt.workspace)
pdf_2d = w.pdf(opt.pdf)
data_obs = w.data('data_obs')
mgg = w.var("mgg")
mjj = w.var("mjj")

### print pdf parameters
params = pdf_2d.getParameters(data_obs)
params.Print()
iter = params.createIterator()
var = iter.Next()
while var :
    print var.GetName(), ' - ', var.getVal()
    var = iter.Next()

h_pdf_2d = pdf_2d.createHistogram("hh_model",mgg,RooFit.Binning(80),RooFit.YVar(mjj,RooFit.Binning(40)))  #pdf_2d.createHistogram(mgg,mjj, 80, 40)
h_pdf_2d_norm = h_pdf_2d.Integral()
h_pdf_2d.Scale(1./h_pdf_2d_norm)

titles = opt.label.split(';')
h_pdf_2d.SetTitle(titles[0])
if len(titles) > 1: h_pdf_2d.GetXaxis().SetTitle(titles[1])
if len(titles) > 2: h_pdf_2d.GetYaxis().SetTitle(titles[2])
h_pdf_2d.SetName('pdf2d_'+opt.name)
c = TCanvas('c', 'c', 800, 600)
h_pdf_2d.Draw("COLZ")
c.SaveAs('pdf2d_'+opt.name+'.pdf')
c.SaveAs('pdf2d_'+opt.name+'.png')

tout = TFile('root_'+opt.name+'.root', 'RECREATE')
tout.cd()
h_pdf_2d.Write()
tout.Close()
