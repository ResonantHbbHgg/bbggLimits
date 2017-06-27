from ROOT import *
import os,sys, os.path
import getpass
from HiggsAnalysis.bbggLimits.DefineScans import *
username = getpass.getuser()


HERE = os.environ['PWD']

scan = scan_2d

FOLDER = 'LIMS_LT_350_HMHPC_970_HMMPC_600_LMHPC_985_LMMPC_600_v66'
EOSLOC = 'root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/resonant_HH/RunII/LimitTrees/2016/'
L_OUTLOC_org = '/tmp/'+username+'/'

counter = 0
nbad = 0
for kl in scan['kl']:
 for kt in scan['kt']:
  for cg in scan['cg']:
   for c2 in scan['c2']:
    for c2g in scan['c2g']:
#     print counter
     counter += 1

     pointStr = 'kl'+str(kl).replace('.', 'p').replace('-', 'm')+'_kt'+str(kt).replace('.', 'p').replace('-', 'm')+'_cg'+str(cg).replace('.', 'p').replace('-', 'm')+'_c2'+str(c2).replace('.', 'p').replace('-', 'm')+'_c2g'+str(c2g).replace('.', 'p').replace('-', 'm')

     hm_f = FOLDER+"/HighMass_Node_SM"+pointStr
     lm_f = FOLDER+"/LowMass_Node_SM"+pointStr
     co_f = FOLDER+"/CombinedCard_Node_SM"+pointStr
     hm_f_b = 1
     lm_f_b = 1
     co_f_b = 1
     if os.path.isdir(hm_f) == False : hm_f_b=0
     if os.path.isdir(lm_f) == False : lm_f_b=0
     if os.path.isdir(co_f) == False : co_f_b=0
     if hm_f_b == 0 or lm_f_b == 0 or co_f_b == 0:
       print kl, kt, cg, c2, c2g,'High mass dir', hm_f_b, 'Low mass dir', lm_f_b, 'Comb dir', co_f_b
       nbad += 1
       continue

     lm_c = lm_f + '/datacards/higgsCombineLowMass_Node_SM'+pointStr+'.Asymptotic.mH125.root'
     hm_c = hm_f + '/datacards/higgsCombineHighMass_Node_SM'+pointStr+'.Asymptotic.mH125.root'
     co_c = co_f + '/higgsCombine_Node_SM'+pointStr+'.Asymptotic.mH125.root'
     hm_c_b = 1
     lm_c_b = 1
     co_c_b = 1
     if os.path.isfile(hm_c) == False : hm_c_b=0
     if os.path.isfile(lm_c) == False : lm_c_b=0
     if os.path.isfile(co_c) == False : co_c_b=0
     if hm_c_b == 0 or lm_c_b == 0 or co_c_b == 0:
       print kl, kt, cg, c2, c2g,'High mass lim file', hm_c_b, 'Low mass lim file', lm_c_b, 'Comb lim file', co_c_b
       nbad += 1
       continue

print "Number of bad points", nbad
