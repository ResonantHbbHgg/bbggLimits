#!/usr/bin/env python

from ROOT import *
from HiggsAnalysis.bbggLimits.SigPlotter import *
gROOT.SetBatch()

import argparse
parser =  argparse.ArgumentParser(description='Signal shape plot maker')
parser.add_argument("limdir")

opt = parser.parse_args()

pointSM = 'ARW_kl_1p0_kt_1p0_cg_0p0_c2_0p0_c2g_0p0'
D_LM = opt.limdir+"/LowMass_"+pointSM
D_HM = opt.limdir+"/HighMass_"+pointSM

lumi = '35.9'
doBands = 1
name = "#font[61]{pp#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}|SM-like HH"
Label = "SMHHLM"
DSCB = True
bins = [24,160]
xmax = {'mgg':135, 'mjj': 190}
xmin = {'mgg':118, 'mjj': 70}

if __name__ == "__main__":

    gSystem.Load("libHiggsAnalysisCombinedLimit.so")

    for c in [0,1,2,3]:
        ccat = c
        if c in [0,1]:
            D_MCAT = D_LM
            analysis=name+'|Low-mass region'
        else:
            D_MCAT = D_HM
            analysis=name+'|High-mass region'
            if c == 2:
                ccat = 0
            if c == 3:
                ccat = 1

        wfile = D_MCAT+'/workspaces/hhbbgg.mH125_13TeV.inputsig.root'
        wroot = TFile(wfile, "READ")
        workspace = wroot.Get("w_all")

        print 'Doing %r, cat=%r, ccat=%r' % (D_MCAT, c, ccat)
        #wroot.Print("all")
        #workspace.Print('all')

        for i,ob in enumerate(['mjj','mgg']):
            data2D = workspace.data("Sig_cat"+str(c))
            data2D.Print()
            pdf = workspace.pdf(ob+"Sig_cat"+str(ccat)+"_CMS_sig_cat"+str(c))

            var = workspace.var(ob)
            data = data2D.reduce(RooArgSet(var))

            label_ax = "m_{jj} [GeV]"
            if 'mgg' in ob:
                label_ax = "m_{#gamma#gamma} [GeV]"
            MakeSigPlot(data, pdf, var, label_ax, lumi, str(c), analysis, doBands, Label+"_signal_fit_"+ob+"_cat"+str(c), bins[i], xmin[ob], xmax[ob], 1, DSCB, outPath=opt.limdir)
