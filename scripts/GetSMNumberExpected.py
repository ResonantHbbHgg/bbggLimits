from ROOT import *
from os import listdir
from os.path import isfile, join
import argparse
parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument("-f", "--folder", dest="folder", type=str)
opt = parser.parse_args()

onlyfiles = [f for f in listdir(opt.folder) if isfile(join(opt.folder, f))]

mysamps = {
'LT_DoubleEG.root':"Observed", 
'LT_output_GluGluToHHTo2B2G_node_SM_13TeV-madgraph.root':"SM $HH\\rightarrow b\\bar{b}\gamma\gamma$", 
'LT_output_GluGluHToGG_M-125_13TeV_powheg_pythia8.root':"SM $H\\rightarrow\gamma\gamma$ (Gluon fusion)", 
'LT_output_VBFHToGG_M-125_13TeV_powheg_pythia8.root':"SM $H\\rightarrow\gamma\gamma$ (VBF)",
'LT_output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root':"SM $H\\rightarrow\gamma\gamma$ (VH)",
'LT_output_ttHToGG_M125_13TeV_powheg_pythia8_v2.root':"SM $H\\rightarrow\gamma\gamma$ ($t\\bar{t}$H)", 
'LT_output_bbHToGG_M-125_13TeV_amcatnlo.root':"SM $H\\rightarrow\gamma\gamma$ ($b\\bar{b}$H)"
}

myyields = {}

print opt.folder, onlyfiles

for f in onlyfiles:
  if f not in mysamps: continue
  tf = TFile(opt.folder+'/'+f)
  tt = tf.Get("TCVARS")
  h_cat0 = TH1F("h_cat0", "", 10, -1, 3)
  h_cat1 = TH1F("h_cat1", "", 10, -1, 3)
  tt.Draw("cut_based_ct>>h_cat0", "evWeight*(cut_based_ct==0)", "goff")
  tt.Draw("cut_based_ct>>h_cat1", "evWeight*(cut_based_ct==1)", "goff")
  cat0 = h_cat0.Integral()
  cat1 = h_cat1.Integral()
  myyields[mysamps[f]] = [cat0,cat1]
  smw = 8.697000e-02
  if 'node_SM' in f:
    myyields[mysamps[f]] = [cat0*smw,cat1*smw]

print '\\begin{tabular}{c | c | c}'
print 'Sample \t&\t HPC \t&\t MPC \\\\ \\hline'
totsm0 = 0
totsm1 = 0
for ff in mysamps:
  f = mysamps[ff]
  if 'HH' in f: continue
  print f, '\t&\t', "{:.3f}".format(myyields[f][0]), '\t&\t', "{:.3f}".format(myyields[f][1]), '\t\t \\\\'
  if "Obs" not in f:
    totsm0 += myyields[f][0]
    totsm1 += myyields[f][1]
print '\\hline'
print 'Total SM H \t&\t', "{:.3f}".format(totsm0), '\t&\t', "{:.3f}".format(totsm1), '\t\t \\\\'
print '\\hline'
print "SM $HH\\rightarrow b\\bar{b}\gamma\gamma$ \t&\t", "{:.3f}".format(myyields["SM $HH\\rightarrow b\\bar{b}\gamma\gamma$"][0]), '\t&\t', "{:.3f}".format(myyields["SM $HH\\rightarrow b\\bar{b}\gamma\gamma$"][1]), '\t\t \\\\'
print '\\end{tabular}'
