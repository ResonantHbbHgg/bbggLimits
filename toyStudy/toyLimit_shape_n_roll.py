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
                    choices=['Gauss', 'CB', '2D', '2D-CB'], help = "Choose the type of toy limits to run")

parser.add_argument('--mjjCut', dest="mjjCut", default=None, type=str,
                    help = "A cut on mjj for the case of 1D fit. Must be a valid string like this: 'mjj>100 && mjj<190 \
                    Also accepted strings: eff_sigma_1, eff_sigma_2")

opt = parser.parse_args()

# ------
# Create a workspace to be used for toy limits.
# It will contain 2 PDFs: bkg and sig, and a data_obs DataSet of fake data.
#-------

w = RooWorkspace("w")

w.factory('Exponential::bkg_mgg(mgg[100,180],tau_mgg[-0.04])')
w.factory('Exponential::bkg_mjj(mjj[70,190],tau_mjj[-0.03])')

#w.factory('Chebychev::bkg_mgg(mgg[100,180], {-1.0, 0.2, 0.1})')
#w.factory('Bernstein::bkg_mjj(mjj[70,190], {4.32, -0.83, -1.75})')

if opt.toyType=="Gauss" or opt.toyType=="2D":
    w.factory('Gaussian::sig_mgg(mgg, 124.8, sigma_mgg[1.6])')
    w.factory('Gaussian::sig_mjj(mjj, 123.0, sigma_mjj[15])')

elif opt.toyType=="CB" or opt.toyType=="2D-CB":
    # The parameters of th Double Sided Crystal Ball are taken from the analysis workspace (see readWS.py script)
    # This one is after the bug fix:
    w.factory('DoubleCB::sig_mgg(mgg, 124.78, sigma_mgg[1.33], 1.18, 5.25, 1.65, 9.38)')
    # This one is before the bug fix
    # w.factory('DoubleCB::sig(x, 124.85, 0.79, 1.02, 3.75, 1.68, 4.63)')

    # With parameters from our signal MC sample:
    w.factory('DoubleCB::sig_mjj(mjj, 123.53, sigma_mjj[12], 0.70, 2.06, 1.81, 0.95)')
    # Reduced the width:
    # w.factory('DoubleCB::sig_mjj(mjj, 123.53, sigma_mjj[9.6], 1.20, 0.96, 1.81, 0.95)')

w.factory('PROD::bkg_2D(bkg_mgg, bkg_mjj)')
w.factory('PROD::sig_2D(sig_mgg, sig_mjj)')

# ---------
# Re-naming. The final PDFs should be named sig and bkg (used in the datacards)
if '2D' in opt.toyType:
    w.factory('EDIT::bkg(bkg_2D, tau_mgg=tau_mgg)')
    w.factory('EDIT::sig(sig_2D, sigma_mgg=sigma_mgg)')
else:
    w.factory('EDIT::bkg(bkg_mgg, tau_mgg=tau_mgg)')
    w.factory('EDIT::sig(sig_mgg, sigma_mgg=sigma_mgg)')
    
x = w.var('mgg')
y = w.var('mjj')
sigPdf    = w.pdf("sig")
bkgPdf    = w.pdf("bkg")

# Generate fake dataset:
fakeData  = bkgPdf.generate(RooArgSet(x,y), 150)

# Save all to the warkspace root file
getattr(w,'import')(fakeData,RooFit.Rename("data_obs"))
w.Print()
w.writeToFile('ws.root')

# Workspace is saved


# Now let's make some plots
c = TCanvas("c","c",0,0,900,600)
c.cd()

for v in [x,y]:
    testFrame = v.frame()
    fakeData.plotOn(testFrame, RooFit.Binning(80), RooFit.Name('Toy data'))
    bkgPdf.plotOn(testFrame, RooFit.Name("Toy Bkg"), RooFit.LineColor(kBlue+1), RooFit.LineWidth(2))
    sigPdf.plotOn(testFrame, RooFit.Name("Toy Sig"), RooFit.LineColor(kRed+2),  RooFit.LineWidth(2), RooFit.Normalization(0.10))
    testFrame.Draw()

    c.SaveAs('tmpfig_gen_bkg_'+v.GetName()+'.png')

    
# ----------------
# Here we define some scans: over sigma values, number of background and signal events:
# ----------------
Sigmas = np.array([0.6, 1.0, 1.6, 2.0])
if opt.toyType in ["CB","2D-CB"]:
    Sigmas = np.array([1.0, 1.33, 2.0]) # This is for CB shape of the mgg in signal

Nbkg = np.array([100, 120, 250])
Nsig = np.array([2, 4, 9], dtype=float)
toyCard2 = 'card_shape_n_roll.txt'

# ----- 
# If there is an mjj cut, we need to ajust the number of events in the datacard.
# This is to compare limits of 2D fit with a limit of 1D fit after an mjj cut
if opt.mjjCut!=None:
    print "The mjj cut is:", opt.mjjCut
    mjjCut = opt.mjjCut
    if 'eff_sigma' in opt.mjjCut:
        from eff_sigma import getEffSigma
        eff_sigma_mjj = getEffSigma(w.var('mjj'), w.pdf('sig_mjj'), 70, 170)
        print eff_sigma_mjj
        if opt.mjjCut=='eff_sigma_1' or opt.mjjCut=='eff_sigma':
            print 'cut at +/- 1 effective sigmas'
            mjjCut = "mjj > %.2f && mjj < %.2f" %(eff_sigma_mjj[1], eff_sigma_mjj[2])
        if opt.mjjCut=='eff_sigma_2':
            print 'cut at +/- 2 effective sigmas'
            mjjCut = "mjj > %.2f && mjj < %.2f" %(eff_sigma_mjj[1]-eff_sigma_mjj[0], eff_sigma_mjj[2]+eff_sigma_mjj[0])
    print mjjCut
    
    for i,si in enumerate(Nsig):
        tmpSig  = w.pdf('sig_2D').generate(RooArgSet(x,y), si*1000)
        reduced = tmpSig.reduce(RooArgSet(x,y), mjjCut).sumEntries()
        print si, tmpSig.sumEntries(), reduced
        Nsig[i] = reduced/1000
    
    for i,N in enumerate(Nbkg):
        tmpBkg  = w.pdf('bkg_2D').generate(RooArgSet(x,y), N)
        reduced = tmpBkg.reduce(RooArgSet(x,y), mjjCut).sumEntries()
        print N, tmpBkg.sumEntries(), reduced
        Nbkg[i] = reduced
        
print Nsig
print Nbkg

# sys.exit()

# ------------
# Now we can submit jobs to run combine many times:

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
        if opt.toyType in ["Gauss", "2D"]:
            # print "Nbkg = ", N, 'NSig = ', si
            lim_s1p0 = limsAll['NBKG_'+str(N)+'_NSIG_'+str(si)][1] # For sigma = 1.0 GeV
            lim_s1p6 = limsAll['NBKG_'+str(N)+'_NSIG_'+str(si)][2] # For sigma = 1.6 GeV

            exp_diff = np.sqrt(1.6)/np.sqrt(1) - 1
            diff = (lim_s1p6 - lim_s1p0)/lim_s1p0
            print "| %.1f| %.1f | %.2f | %.2f | %.3f | %.3f |" % (N, si, lim_s1p0, lim_s1p6, diff, exp_diff)

        elif opt.toyType in ["CB","2D-CB"]:
            lim_0 = limsAll['NBKG_'+str(N)+'_NSIG_'+str(si)][1]
            print "| %.1f | %.1f | %.2f |" % (N, si, lim_0)

            tmp.append(lim_0)
print tmp
