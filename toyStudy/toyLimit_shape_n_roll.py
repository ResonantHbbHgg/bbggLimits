#! /usr/bin/env python
import os,sys,time
import tempfile, getpass, subprocess
import numpy as np
import matplotlib.pyplot as plt

username = getpass.getuser()
from ROOT import *
gROOT.SetBatch()


import argparse
parser =  argparse.ArgumentParser(description='Shape limits toy game')
parser.add_argument('--dry', dest="dry", action="store_true", default=False,
                    help="Dry run, without running the limits. Assume the output files of combine exist, so just make the plots.")
parser.add_argument('-t','--type', dest="toyType", default='Gauss', type=str,
                    choices=['Gauss', 'CB', '2D'], help = "Choose the type of toy limits to run")
opt = parser.parse_args()

# ------
# Create a workspace to by used for toy limits.
# It will contain 2 PDFs: bkg and sig, and a data_obs DataSet for fake data. 
#-------
w = RooWorkspace("w")
w.factory('Exponential::bkg(x[100,180],tau[-0.04])')
#w.factory('Chebychev::bkg(x,{-0.4,0.1})')
if opt.toyType=="Gauss":
    w.factory('Gaussian::sig(x, 124.8, sigma[1])')

elif opt.toyType=="CB":
    # The parameters of th Double Sided Crystal Ball are taken from the analysis workspace (see readWS.py script)
    # This one is after the bug fix:
    w.factory('DoubleCB::sig(x, 124.78, 1.33, 1.18, 5.25, 1.65, 9.38)')

    # This one is before the bug fix
    # w.factory('DoubleCB::sig(x, 124.85, 0.79, 1.02, 3.75, 1.68, 4.63)')
    
elif opt.toyType=="2D":
    print "This is not yet implemented"
    # w.factory('DoubleCB::sig(x, 124.78, 1.33, 1.18, 5.25, 1.65, 9.38)')
    
x = w.var('x')
sigPdf    = w.pdf("sig")
bkgPdf    = w.pdf("bkg")
fakeData  = bkgPdf.generate(RooArgSet(x), 120)
getattr(w,'import')(fakeData,RooFit.Rename("data_obs"))
w.Print() 
w.writeToFile('ws.root')

# Workspace is saved


# Now let's make a plot
c = TCanvas("c","c",0,0,900,600)
c.cd()

testFrame = x.frame()
fakeData.plotOn(testFrame, RooFit.Binning(80), RooFit.Name('data2'))
bkgPdf.plotOn(testFrame, RooFit.Name("Bkg 2"), RooFit.LineColor(kBlue), RooFit.LineWidth(2), RooFit.FillColor(kCyan-6))
sigPdf.plotOn(testFrame, RooFit.Name("fake sig"), RooFit.LineColor(kRed+2), RooFit.LineWidth(2), RooFit.Normalization(0.10))
testFrame.Draw()

c.SaveAs('tmpfig_gen_bkg.png')


# Now we can submit jobs to run combine many times:

Sigmas = np.array([0.6, 0.8, 1.0, 1.4, 1.6, 2.0])
if opt.toyType != "Gauss":
    Sigmas = np.array([1.0]) # This is just to avoid running multiple times for CB. 
    
Nbkg = np.array([100, 120, 140])
Nsig = np.array([2, 4, 6])
toyCard2 = 'card_shape_n_roll.txt'


processes = []
for s in Sigmas:     
    with open(toyCard2, 'r') as templateCard:
        tmpCard = templateCard.read()
        for si in Nsig:
            for N in Nbkg:
                with tempfile.NamedTemporaryFile(dir='/tmp/'+username, prefix='toyCard_N'+str(N), suffix='.txt', delete=False) as card:
                    newCard = tmpCard.replace('%BKG%',str(N)).replace('%SIG%',str(si)).replace('%SIGMA%',str(s))
                    
                    card.write(newCard)
                    if not opt.dry:
                        print card.name
                        pro = subprocess.Popen(['combine','-M','Asymptotic', '-m', '125','-n','_Shape_NBKG_'+str(N)+'_NSIG_'+str(si)+'_sigma_'+str(s), card.name])
                        processes.append(pro)
                        time.sleep(0.2)
                
                    
# Wait till all processes are finished:
for p in processes:
    print 'waiting for', p
    p.wait()



limsAll = {}
                
for si in Nsig:
    for N in Nbkg:
        lims = []
        for s in Sigmas:
            fname = 'higgsCombine_Shape_NBKG_'+str(N)+'_NSIG_'+str(si)+'_sigma_'+str(s)+'.Asymptotic.mH125.root'
            rfile = TFile(fname, 'READ')
            tree = rfile.Get("limit")
            tree.Draw("limit:quantileExpected", "", "goff")
            
            lims.append(tree.GetV1()[2])
            rfile.Close()
            # print si, N, s,  lims[-1]
        limsAll['NBKG_'+str(N)+'_NSIG_'+str(si)] = lims

        #print 'Expected limit at r=0.5: ', lims
print limsAll


import itertools
marker = itertools.cycle((',', '+', '.', 'o', '*', '^'))

plt.clf()
for si in Nsig:
    for N in Nbkg:
        plt.plot(np.sqrt(Sigmas), limsAll['NBKG_'+str(N)+'_NSIG_'+str(si)], marker = marker.next(), label='$N_{bkg} = '+str(N)+'; N_{sig} = '+str(si)+'$')
plt.axis([0., 2.5, 0, 8])
plt.title('Limits from toy shapes')
plt.xlabel(r'$\sqrt{\sigma_{eff}(signal)}}$')
plt.ylabel('Median expected Limit')
plt.legend()
plt.savefig('tmpfig_scale_with_sigma.png')


tmp = []
for si in Nsig:
    for N in Nbkg:
        if opt.toyType=="Gauss":
            # print "Nbkg = ", N, 'NSig = ', si
            lim_s1p0 = limsAll['NBKG_'+str(N)+'_NSIG_'+str(si)][2] # For sigma = 1.0 GeV
            lim_s1p6 = limsAll['NBKG_'+str(N)+'_NSIG_'+str(si)][4] # For sigma = 1.6 GeV
            
            exp_diff = np.sqrt(1.6)/np.sqrt(1) - 1
            diff = (lim_s1p6 - lim_s1p0)/lim_s1p0
            print "| %i | %i | %.2f | %.2f | %.3f | %.3f |" % (N, si, lim_s1p0, lim_s1p6, diff, exp_diff)        

        elif opt.toyType=="CB":
            lim_0 = limsAll['NBKG_'+str(N)+'_NSIG_'+str(si)][0]
            print "| %i | %i | %.2f |" % (N, si, lim_0)        

            tmp.append(lim_0)
print tmp
