from ROOT import *
import argparse, os
from HiggsAnalysis.bbggLimits.DefineScans import *

parser =  argparse.ArgumentParser(description='Benchmark plot maker')
parser.add_argument("limdir")
parser.add_argument('-s','--symmetry', dest="symmetry", action="store_true", default=False,
                    help="For negative kt values take the positive ones: they are symmetric over (kl,kt)=(0,0)")
opt = parser.parse_args()

quantiles = ['0.025', '0.160', '0.500', '0.840', '0.975', '-1']
lims = {}
for qt in quantiles:
  lims[qt] = []

myKl = []

listFileName = opt.limdir+'/KLKT_Scan_List'
worked = open(listFileName+".txt", "w+")
not_worked = open(listFileName + "_not_worked.txt", "w+")

counter = 0
#npts = 0

myScan = scan_2d

for kl in myScan['kl']:
  for kt in myScan['kt']:
    if opt.symmetry:
      if kt<=0: continue
      
    #    if npts > 50: break 
    #    npts += 1
    pointStr = ('ARW_kl_' + str(float(kl)) + '_kt_'+str(float(kt))+'_cg_0p0_c2_0p0_c2g_0p0').replace('.', 'p').replace('-', 'm')
    fname = opt.limdir + '/CombinedCard_' + pointStr + '/higgsCombine_' + pointStr + '.Asymptotic.mH125_1.root'

    tfile = TFile(fname, "READ")
    thesePoints = str(kl) + ' ' + str(kt) + ' 0.0 0.0 0.0'
    if tfile.IsZombie() == 1:
      not_worked.write(thesePoints+'\t isZombie \n')
      continue
    tree = tfile.Get('limit')
    if tree == None:
      not_worked.write(thesePoints+'\t tree is None \n')
      continue
    if tree.GetEntries() < 2:
      not_worked.write(thesePoints+'\t tree not filled correctly \n')
      continue
    myKl.append(kl)
    for qt in quantiles:
      tree.Draw("limit", "quantileExpected>"+str(float(qt)-0.001) + ' && quantileExpected < ' +str(float(qt)+0.001), "goff")
      lims[qt].append(tree.GetV1()[0])
    tfile.Close()
    if lims['0.500'][-1] == 0 or lims['-1'][-1] == 0:
      print '#######################################################', lims['0.500'][-1], lims['-1'][-1]
      not_worked.write(thesePoints+'\t Expected or Observed limit is 0 \n')
      continue
    # ipt kl kt exp obs exp_p1s exp_m1s exp_p2s exp_m2s

    str_worked = ' '.join([str(counter), str(kl), str(kt), str(lims['0.500'][-1]), str(lims['-1'][-1]),
                           str(lims['0.840'][-1]), str(lims['0.160'][-1]), str(lims['0.975'][-1]), str(lims['0.025'][-1]), '\n'])
    worked.write(str_worked)
    counter += 1

    if opt.symmetry:
      str_worked = ' '.join([str(counter), str(-kl), str(-kt), str(lims['0.500'][-1]), str(lims['-1'][-1]),
                             str(lims['0.840'][-1]), str(lims['0.160'][-1]), str(lims['0.975'][-1]), str(lims['0.025'][-1]), '\n'])
      worked.write(str_worked)
      counter += 1
      

    
worked.close()
not_worked.close()

print 'The list of points we used to make the plot is here:', listFileName+'.txt'
print 'You should use it as input to MakeKLKTplot.py script:'
print '   python scripts/MakeKLKTplot.py -l '+listFileName+'.txt'
