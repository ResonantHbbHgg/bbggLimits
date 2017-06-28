import os,sys
import argparse
import glob
from ROOT import *
parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument("-f", "--folder", dest="folder", default="900", type=str)
parser.add_argument('--plotOnly', dest='plotOnly', action='store_true', default=False)
parser.add_argument('--unblind', dest='ub', action='store_true', default=False)
opt = parser.parse_args()

#build new workspace with masking
#../LIMS_LT_350_HMHPC_990_HMMPC_680_LMHPC_980_LMMPC_475_v66/CombinedCard_Node_SMLT_350_HMHPC_990_HMMPC_680_LMHPC_980_LMMPC_475/hhbbgg_13TeV_DataCard.txt
dcName = glob.glob(opt.folder+"/hhbbgg_13TeV_DataCard.txt")[0]
fdName = opt.folder#glob.glob(opt.folder+"/Comb*")[0]

rtName = dcName.replace(".txt", ".root")
print rtName

if not opt.plotOnly:
  command = 'text2workspace.py '+ dcName +' --channel-masks'
  os.system(command)
  #mask_ch1_cat2,mask_ch1_cat3,mask_ch2_cat0,mask_ch2_cat1
  command = 'combine -n LMHPC -M Asymptotic --X-rtd TMCSO_AdaptivePseudoAsimov=50 --setPhysicsModelParameters mask_ch2_cat0=0,mask_ch2_cat1=1,mask_ch1_cat2=1,mask_ch1_cat3=1 --datacard ' + rtName
  os.system(command)
  command = 'mv higgsCombineLMHPC.Asymptotic.mH120.root ' + fdName + '/higgsCombineLMHPC.Asymptotic.mH120.root'
  os.system(command)

  command = 'combine -n LMMPC -M Asymptotic --X-rtd TMCSO_AdaptivePseudoAsimov=50 --setPhysicsModelParameters mask_ch2_cat0=1,mask_ch2_cat1=0,mask_ch1_cat2=1,mask_ch1_cat3=1 --datacard ' + rtName
  os.system(command)
  command = 'mv higgsCombineLMMPC.Asymptotic.mH120.root ' + fdName + '/higgsCombineLMMPC.Asymptotic.mH120.root'
  os.system(command)

  command = 'combine -n HMHPC -M Asymptotic --X-rtd TMCSO_AdaptivePseudoAsimov=50 --setPhysicsModelParameters mask_ch2_cat0=1,mask_ch2_cat1=1,mask_ch1_cat2=0,mask_ch1_cat3=1 --datacard ' + rtName
  os.system(command)
  command = 'mv higgsCombineHMHPC.Asymptotic.mH120.root ' + fdName + '/higgsCombineHMHPC.Asymptotic.mH120.root'
  os.system(command)

  command = 'combine -n HMMPC -M Asymptotic --X-rtd TMCSO_AdaptivePseudoAsimov=50 --setPhysicsModelParameters mask_ch2_cat0=1,mask_ch2_cat1=1,mask_ch1_cat2=1,mask_ch1_cat3=0 --datacard ' + rtName
  os.system(command)
  command = 'mv higgsCombineHMMPC.Asymptotic.mH120.root ' + fdName + '/higgsCombineHMMPC.Asymptotic.mH120.root'
  os.system(command)

'''
parser.add_argument('-f', '--inputFolders', dest="folders", default=None, type=str, nargs='+', required=True,
                    help="input folders")
parser.add_argument('-n', '--inputNames', dest="names", default=None, type=str, nargs='+', required=True,
                    help="input folders cat names")
parser.add_argument('-l', '--lumi', dest='lumi', default='36.5', type=str, help='Integrated luminosoty')
parser.add_argument('--label', dest='label', default='', type=str, help='Label')
parser.add_argument('--log', dest='log', action='store_true', default=False)
parser.add_argument('--isAsymptotic', dest='asymp', action='store_true', default=False)
parser.add_argument('--normSM', dest='normSM', action='store_true', default=False)
parser.add_argument('--unblind', dest='unblind', action='store_true', default=False)
'''

lmmpc = fdName + '/higgsCombineLMMPC '
lmhpc = fdName + '/higgsCombineLMHPC '
lmcom = opt.folder.replace('CombinedCard', 'LowMass')+'/datacards/ '
hmmpc = fdName + '/higgsCombineHMMPC '
hmhpc = fdName + '/higgsCombineHMHPC '
hmcom = opt.folder.replace('CombinedCard', 'HighMass')+'/datacards/ '
allco = fdName + '/higgsCombineCombinedCard_Node '

unblind = ''
if opt.ub: unblind = ' --unblind '

command = 'scripts/MakeSMHHCatsPlot.py -f ' + lmmpc + lmhpc + lmcom + hmmpc + hmhpc + hmcom + allco + unblind + ' -n LM-MP LM-HP "LM Comb." HM-MP HM-HP "HM Comb." "All Comb." -l 35.9 --log --isAsymptotic '
os.system(command)
command = 'mv test_SM.pdf ' + opt.folder + '/HHSM_CategoryBreakdown.pdf'
os.system(command)

command = 'scripts/MakeSMHHCatsPlot.py -f ' + lmmpc + lmhpc + lmcom + hmmpc + hmhpc + hmcom + allco + unblind + ' -n LM-MP LM-HP "LM Comb." HM-MP HM-HP "HM Comb." "All Comb." -l 35.9 --log --isAsymptotic --normSM '
os.system(command)
command = 'mv test_SM.pdf ' + opt.folder + '/HHSM_CategoryBreakdown_norm.pdf'
os.system(command)


