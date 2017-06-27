from ROOT import *
from HiggsAnalysis.bbggLimits.AnalyticalReweighting import *
import argparse, sys

parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument('-f', '--inputFile', dest="fname", type=str, default=None, required=True,
                    help="Json config file")
parser.add_argument('-o', '--outDir', dest="outDir", type=str, default=".",
                    help="Output directory (will be created).")
parser.add_argument('--kl', dest='ARW_kl', type=float, default=1.0)
parser.add_argument('--kt', dest='ARW_kt', type=float, default=1.0)
parser.add_argument('--cg', dest='ARW_cg', type=float, default=0.0)
parser.add_argument('--c2', dest='ARW_c2', type=float, default=0.0)
parser.add_argument('--c2g', dest='ARW_c2g', type=float, default=0.0)

opt = parser.parse_args()

TempSignalFile = opt.fname
TempFile = TFile(TempSignalFile)
TempTree = TempFile.Get("bbggSelectionTree")
if TempTree == None :
  TempTree = TempFile.Get("TCVARS")
if TempTree == None : sys.exit()

OutFile = opt.outDir+"/LT_NR_Nodes_All_merged_kl_"+str(opt.ARW_kl).replace(".", "p")+"_kt_"+str(opt.ARW_kt).replace(".", "p")+"_cg_"+str(opt.ARW_cg).replace(".", "p")+"_c2_"+str(opt.ARW_c2).replace(".", "p")+"_c2g_"+str(opt.ARW_c2g).replace(".", "p")+".root"

AddReWeightBranch(TempTree, OutFile, opt.ARW_kl, opt.ARW_kt, opt.ARW_cg, opt.ARW_c2g, opt.ARW_c2)
TempFile.Close()
