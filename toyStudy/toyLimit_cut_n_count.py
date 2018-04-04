#! /usr/bin/env python
import os,sys,time
import tempfile, getpass, subprocess
import numpy as np
import matplotlib.pyplot as plt

username = getpass.getuser()
from ROOT import *
gROOT.SetBatch()

Nbkg = np.array([2, 4, 6, 6.4, 8, 10, 15, 20, 30, 50])
Nsig = np.array([2, 4, 6])

toyCard1 = 'card_cut_n_count.txt'

import argparse
parser =  argparse.ArgumentParser(description='Cut and count toy game')
parser.add_argument('--dry', dest="dry", action="store_true", default=False,
                    help="Dry run, without running the limits. Assume the output files of combine exist, so just make the plots.")
opt = parser.parse_args()

processes = []
with open(toyCard1, 'r') as templateCard:
    tmpCard = templateCard.read()
    for si in Nsig:
        for N in Nbkg:
            with tempfile.NamedTemporaryFile(dir='/tmp/'+username, prefix='toyCard_N'+str(N), suffix='.txt', delete=False) as card:
                newCard = tmpCard.replace('%BKG%',str(N)).replace('%SIG%',str(si))
                # print N
                #print newCard

                card.write(newCard)
                if not opt.dry:
                    print card.name
                    pro = subprocess.Popen(['combine','-M','Asymptotic', '-m', '125','-n','_CutCount_NBKG_'+str(N)+'_NSIG_'+str(si), card.name])
                    processes.append(pro)
                    time.sleep(0.2)

# Wait till all processes are finished:
for p in processes:
    print 'waiting for', p
    p.wait()

# Now open output files and read the limit values:
    
limsAll = {}

for si in Nsig:
    lims = []
    for N in Nbkg:
        fname = 'higgsCombine_CutCount_NBKG_'+str(N)+'_NSIG_'+str(si)+'.Asymptotic.mH125.root'
        rfile = TFile(fname, 'READ')
        tree = rfile.Get("limit")
        tree.Draw("limit:quantileExpected", "", "goff")
        
        #for v in range(0, tree.GetSelectedRows()):
        #    print v, tree.GetV1()[v], tree.GetV2()[v]
        # Quantiles are: [0]=0.025, [1]=0.16, [2]=0.5, [3]0.84, [4]=0.875, [5]=-1 (observed)
        lims.append(tree.GetV1()[2])
        rfile.Close()
        print N, lims[-1]
    limsAll[str(si)] = np.array(lims)

    print 'Expected limit at r=0.5: ', lims
print limsAll

plt.plot(np.sqrt(Nbkg), limsAll['2'], 'ro', label="N(sig) = 2 ev")
plt.plot(np.sqrt(Nbkg), limsAll['4'], 'bs', label="N(sig) = 4 ev")
plt.plot(np.sqrt(Nbkg), limsAll['6'], 'g^', label="N(sig) = 6 ev")
plt.axis([0, 8, 0, 8])
plt.title('Limits from toy counting experiments')
plt.xlabel(r'$\sqrt{N_{bkg}}$')
plt.ylabel('Median expected Limit')
plt.legend()

plt.savefig('tmpfig_cut_n_count.png')

for si in Nsig:
    print "Nsig = ", si
    # Pick the numbers close to our real analysis 
    lim_1 = limsAll[str(si)][1] # For Nbkg = 4.0 ev
    lim_2 = limsAll[str(si)][3] # For Nbkg = 6.4 ev

    exp_diff = np.sqrt(6.4)/np.sqrt(4) - 1
    diff = (lim_2 - lim_1)/lim_1
    print "| %i | %.2f | %.2f | %.3f | %.3f |" % (si, lim_1, lim_2, diff, exp_diff)
    
    
    # Pick the numbers for large background
    lim_1 = limsAll[str(si)][-2] # For Nbkg = 30 ev
    lim_2 = limsAll[str(si)][-1] # For Nbkg = 50 ev

    exp_diff = np.sqrt(50)/np.sqrt(30) - 1
    diff = (lim_2 - lim_1)/lim_1
    print "| %i | %.2f | %.2f | %.2f | %.3f |" % (si, lim_1, lim_2, diff, exp_diff)
    
