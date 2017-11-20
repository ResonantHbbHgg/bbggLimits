from ROOT import *
import argparse, os
from HiggsAnalysis.bbggLimits.DefineScans import *

outf = "KLKT_Scan_List.txt"

parser =  argparse.ArgumentParser(description='Benchmark plot maker')
parser.add_argument("limdir")
opt = parser.parse_args()


quantiles = ['0.025', '0.160', '0.500', '0.840', '0.975', '-1']
lims = {}
for qt in quantiles:
  lims[qt] = []

myKl = []

worked = open(outf+".txt", "w+")
not_worked = open(outf + "_not_worked.txt", "w+")

counter = 0
#npts = 0

myScan = scan_2d

for kl in myScan['kl']:
  for kt in myScan['kt']:
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
    if lims['0.500'][counter] == 0 or lims['-1'][counter] == 0:
      print '#######################################################', lims['0.500'][counter], lims['-1'][counter]
      not_worked.write(thesePoints+'\t Expected or Observed limit is 0 \n')
      continue
    # ipt kl kt exp obs exp_p1s exp_m1s exp_p2s exp_m2s
    str_worked = str(counter) + ' '
    str_worked += str(kl) + ' '
    str_worked += str(kt) + ' '
    str_excl = str_worked
    str_worked += str(lims['0.500'][counter]) + ' '
    str_worked += str(lims['-1'][counter]) + ' ' 
    str_worked += str(lims['0.840'][counter]) + ' ' 
    str_worked += str(lims['0.160'][counter]) + ' '
    str_worked += str(lims['0.975'][counter]) + ' ' 
    str_worked += str(lims['0.025'][counter]) + '\n'
    worked.write(str_worked)

    counter += 1

worked.close()
not_worked.close()
