from ROOT import *
from HiggsAnalysis.bbggLimits.DefineScans import *

gROOT.SetBatch()

case = ["High", "Low"]

newF = "/eos/cms/store/group/phys_higgs/resonant_HH/RunII/LimitTrees/ForApproval_2016/LT_350_HMHPC_970_HMMPC_600_LMHPC_985_LMMPC_600_CASEMass/"
oldF = "/eos/cms/store/group/phys_higgs/resonant_HH/RunII/LimitTrees/2016/LT_350_HMHPC_970_HMMPC_600_LMHPC_985_LMMPC_600_CASEMass/"

for cc in case:
  for ii in range(0, len(klJHEP)):
    kl = klJHEP[ii]
    kt = ktJHEP[ii]
    cg = cgJHEP[ii]
    c2 = c2JHEP[ii]
    c2g =c2gJHEP[ii]

    pointStr = 'kl_'+ str(kl).replace('.', 'p') + '_kt_' + str(kt).replace('.', 'p') + '_cg_' + str(cg).replace('.', 'p') + '_c2_' + str(c2).replace('.', 'p') + '_c2g_' + str(c2g).replace('.', 'p')

    fname = 'LT_NR_Nodes_All_merged_' + pointStr + '.root'

    fNew = TFile(newF.replace("CASE", cc)+fname, "READ")
    nTree = fNew.Get("TCVARS")

    fOld = TFile(oldF.replace("CASE", cc)+fname, "READ")
    oTree = fOld.Get("TCVARS")

    if "High" in cc:
      newHist = TH1F("newHist", ";mtot;", 100, 350, 2000)
      oldHist = TH1F("oldHist", ";mtot;", 100, 350, 2000)
    else:
      newHist = TH1F("newHist", ";mtot;", 100, 250, 350)
      oldHist = TH1F("oldHist", ";mtot;", 100, 250, 350)

    oldHist.SetLineColor(kRed)

    nTree.Draw("mtot>>newHist", "(NR_weight)*(cut_based_ct>-1)")
    oTree.Draw("mtot>>oldHist", "(NR_weight)*(cut_based_ct>-1)")

    c = TCanvas("c", "c", 800,600)
    newHist.Draw()
    oldHist.Draw("same")
    c.SaveAs(cc+'_'+pointStr+'.pdf')
