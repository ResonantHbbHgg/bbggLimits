#! /usr/bin/env python
import os,sys, argparse, time
import tempfile, getpass, subprocess
import numpy as np
import matplotlib.pyplot as plt

username = getpass.getuser()
from ROOT import *
gROOT.SetBatch()


# Get the shape of the background from a real card
rooWsFile = TFile('hhbbgg.inputbkg_13TeV.root')
anaWs = rooWsFile.Get('w_all')
anaWs.Print()

mgg = anaWs.var('mgg')
mgg.Print()

anaBkgPdf = anaWs.pdf("mggBkgTmpBer1_cat0_CMS_Bkg_cat2")
fakeData = anaBkgPdf.generate(RooArgSet(mgg), 120)
  
anaBkgPdf.Print()

c = TCanvas("c","c",0,0,900,600)
c.cd()


testFrame = mgg.frame()
fakeData.plotOn(testFrame, RooFit.Binning(80), RooFit.Name('data1'))
anaBkgPdf.plotOn(testFrame, RooFit.Name("Bkg"), RooFit.LineColor(kBlue), RooFit.LineWidth(2), RooFit.FillColor(kCyan-6))
testFrame.Draw()

c.SaveAs('fig_ana.png')
