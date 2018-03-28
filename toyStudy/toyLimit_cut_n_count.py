#! /usr/bin/env python
import os,sys, argparse, time
import tempfile, getpass, subprocess
import numpy as np
import matplotlib.pyplot as plt

username = getpass.getuser()
from ROOT import *
gROOT.SetBatch()

Nbkg = np.array([2,2.5,3,3.5,4,4.5,5,6,8,10,15,20,30])
Nsig = np.array([2, 4, 6])

toyCard1 = 'card_cut_n_count.txt'
processes = []

with open(toyCard1, 'r') as templateCard:
    tmpCard = templateCard.read()
    for si in Nsig:
        for N in Nbkg:
            with tempfile.NamedTemporaryFile(dir='/tmp/'+username, prefix='toyCard_N'+str(N), suffix='.txt', delete=False) as card:
                newCard = tmpCard.replace('%BKG%',str(N)).replace('%SIG%',str(si))
                print N
                #print newCard

                card.write(newCard)
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
    limsAll[str(si)] = lims

    print 'Expected limit at r=0.5: ', lims
print limsAll

plt.plot(np.sqrt(Nbkg), limsAll['2'], 'ro', label="Sig = 2 Ev")
plt.plot(np.sqrt(Nbkg), limsAll['4'], 'bs', label="Sig = 4 Ev")
plt.plot(np.sqrt(Nbkg), limsAll['6'], 'g^', label="Sig = 6 Ev")
plt.axis([0, 8, 0, 8])
plt.title('Limits from toy counting experiments')
plt.xlabel(r'$\sqrt{N_{bkg}}$')
plt.ylabel('Median expected Limit')
plt.legend()

plt.savefig('fig_cut_n_count.png')
