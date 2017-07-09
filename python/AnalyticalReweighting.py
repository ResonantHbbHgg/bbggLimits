from ROOT import *
from HHStatAnalysis.AnalyticalModels.NonResonantModel import NonResonantModel
from array import array
import os

###
## Code adapted from: https://github.com/acarvalh/HHStatAnalysis/blob/templatesForAABBandTTBB/AnalyticalModels/test/nonResonant_test_v0.py
## Thanks to Alexandra!
###

def AddReWeightBranch(intree, newfilename, KL, CT, CG, C2G, C2T):

  cmssw_base =  os.environ['CMSSW_BASE']

  outfile = TFile(newfilename, "RECREATE")

  model = NonResonantModel()
#  dumb = model.ReadCoefficients(cmssw_base+"/src/HHStatAnalysis/AnalyticalModels/data/coefficientsByBin_extended_3M.txt")
  dumb = model.ReadCoefficients(cmssw_base+"/src/HHStatAnalysis/AnalyticalModels/data/coefficientsByBin_extended_3M_costHHSim_19-4.txt")

#  histfilename=cmssw_base+"/src/HiggsAnalysis/bbggLimits/data/NR_AnalyticalWeight/Hist2DSum_V0_SM_box.root"
  histfilename=cmssw_base+"/src/HHStatAnalysis/Support/NonResonant/HistSum2D_4b_rebin_SimCostHH_19-4.root"
#  histtitle= "SumV0_AnalyticalBinExt"
  histtitle= "SumV0_AnalyticalBinExtSimCostHH"
  fileHH=TFile(histfilename)
  sumHAnalyticalBin = fileHH.Get(histtitle)
  print KL, CT, CG, C2G, C2T
  #need to turn kappa top = 0 to kappa top = very small, otherwise weights are 0
  if abs(CT) < 0.00001: CT = 0.00001
  calcSumOfWeights = model.getNormalization(KL,CT,C2T,CG,C2G,histfilename,histtitle)

  outtree = intree.CloneTree(0)
  new_evWeight = array('f', [0])
  NR_weight = array('f', [0])
  NR_kl = array('f', [0])
  NR_ct = array('f', [0])
  NR_cg = array('f', [0])
  NR_c2g = array('f', [0])
  NR_c2t = array('f', [0])
  _new_evWeight = outtree.Branch('new_evWeight', new_evWeight, 'new_evWeight/F')
  _NR_weight = outtree.Branch('NR_weight', NR_weight, 'NR_weight/F')
  _NR_kl = outtree.Branch('NR_kl', NR_kl, 'NR_kl/F')
  _NR_ct = outtree.Branch('NR_ct', NR_ct, 'NR_ct/F')
  _NR_cg = outtree.Branch('NR_cg', NR_cg, 'NR_cg/F')
  _NR_c2t = outtree.Branch('NR_c2t', NR_c2t, 'NR_c2t/F')
  _NR_c2g = outtree.Branch('NR_c2g', NR_c2g, 'NR_c2g/F')
  nentries = intree.GetEntries()

  for i in range(0, nentries):
    if i%1000 == 0: print i
    intree.GetEntry(i)

    NR_kl[0] = KL
    NR_ct[0] = CT
    NR_cg[0] = CG
    NR_c2t[0] = C2T
    NR_c2g[0] = C2G

    mhh = intree.gen_mHH
    cost = intree.gen_cosTheta
    mhhcost = [mhh,cost,0,0]

    bmhh = sumHAnalyticalBin.GetXaxis().FindBin(mhh)
    bcost = sumHAnalyticalBin.GetYaxis().FindBin(abs(cost))

    effSumV0 = sumHAnalyticalBin.GetBinContent(bmhh,bcost)

    #factor of 6 comes from the fact that the normalization of the weights has been calculated
    #with 300k events per node point, in bbgg we only have 6, so we are off by a factor of 30/5=6
    weight = 6*model.getScaleFactor(mhh , cost,KL, CT , C2T , CG , C2G , effSumV0 , calcSumOfWeights)

    NR_weight[0] = weight
    if 'bb' in str(intree.GetName()):
      new_evWeight[0] = (intree.genTotalWeight)*weight
    else:
      new_evWeight[0] = (intree.evWeight)*weight*50000. #need to multiply by 50k since, in the limit trees, evWeight has been normalized by the #evs in MC

    outtree.Fill()

  outfile.cd()
  outtree.Write()
  outfile.Close()

