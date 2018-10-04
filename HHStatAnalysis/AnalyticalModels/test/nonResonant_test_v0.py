#! /usr/bin/env python
# Analytical reweighting implementation for H->4b
# This file is part of https://github.com/cms-hh/HHStatAnalysis.
# python nonResonant_test_v0.py --kl 1 --kt 1 
# compiling
from optparse import OptionParser
import ROOT
import numpy as np
from HHStatAnalysis.AnalyticalModels.NonResonantModel import NonResonantModel

parser = OptionParser()
parser.add_option("--kl", type="float", dest="kll", help="Multiplicative factor in the H trilinear wrt to SM")
parser.add_option("--kt", type="float", dest="ktt", help="Multiplicative factor in the H top Yukawa wrt to SM")

parser.add_option("--c2", type="float", dest="c22", help="ttHH with triangle loop", default=0)
parser.add_option("--cg", type="float", dest="cgg", help="HGG contact", default=0)
parser.add_option("--c2g", type="float", dest="c2gg", help="HHGG contact", default=0)
parser.add_option("--doPlot", action='store_true', default=False, dest='doPlot', 
                    help='calculate the limit in the benchmark poin specified')

(options, args) = parser.parse_args()
print " "
kl = options.kll
kt = options.ktt
c2 = options.c22
cg = options.cgg
c2g = options.c2gg

#print "events for V0 (the same of the fullsim version of Moriond 2016) \n We sum SM + box + the benchmarks from 2-13"
#if c2 != 0 or cg != 0 or c2g != 0 :  print "The analytical function is not yet implemented"

###########################################################
# read events and apply weight
###########################################################
def main():
  #if 1 > 0 :
  # declare the 2D ===> should be global variable
  model = NonResonantModel()
  # obtaining BSM/SM coeficients
  #dumb = model.ReadCoefficients("../data/coefficientsByBin_A1A3A7.txt") 
  dumb = model.ReadCoefficients("../data/coefficientsByBin_extended_3M_costHHSim_19-4.txt") 
  #dumb = model.ReadCoefficients("../data/coefficientsByBin_extended_3M.txt") 
  counteventSM=0
  sumWeight=0
  # now loop over events, calculate weights using the coeffitients and  plot histograms
  # We sum SM + box + the benchmarks from 2-13 
  # read the 2D histo referent to the sum of events
  #fileHH=ROOT.TFile("../../../Analysis/Support/NonResonant/Hist2DSum_V0_SM_box.root")
  #sumHAnalyticalBin = fileHH.Get("SumV0_AnalyticalBin")
  histfilename="../../../Analysis/Support/NonResonant/HistSum2D_4b_rebin_SimCostHH_19-4.root"
  #histtitle= "SumV0_AnalyticalBinExt" #
  histtitle= "SumV0_AnalyticalBinExtSimCostHH" #
  fileHH=ROOT.TFile(histfilename)
  sumHAnalyticalBin = fileHH.Get(histtitle)
  calcSumOfWeights = model.getNormalization(kl, kt,c2,cg,c2g,histfilename,histtitle)  # this input is flexible, tatabb may have only the SM
  # print "sum of weights calculated" , calcSumOfWeights 
  # read the events
  pathBenchEvents="/eos/user/a/acarvalh/asciiHH_tofit/GF_HH_BSM/" #"/eos/user/a/acarvalh/asciiHH_tofit/GF_HH_BSM/" # events to reweight   
  file=ROOT.TFile(pathBenchEvents+"events_SumV0.root")
  tree=file.Get("treeout")
  nev = tree.GetEntries()
  # declare the histograms 
  CalcMhh = np.zeros((nev))
  CalcCost = np.zeros((nev))
  CalcPtH = np.zeros((nev))
  CalcPtHH = np.zeros((nev))
  CalcWeight = np.zeros((nev))

  CalcPtHgg = np.zeros((nev))
  CalcMhhReco = np.zeros((nev))
  CalcMXReco = np.zeros((nev))

  countevent = 0
  #for kll in range(-5,5) : model.getNormalization(kll, kt,sumHBenchBin)
  for iev in range(0,nev) :
      tree.GetEntry(iev)
      mhh = tree.Genmhh
      cost = tree.GenHHCost

      mhhcost= [mhh,cost,0,0] # to store [mhh , cost] of that event
      # find the Nevents from the sum of events on that bin
      bmhh = sumHAnalyticalBin.GetXaxis().FindBin(mhh)
      bcost = sumHAnalyticalBin.GetYaxis().FindBin(abs(cost))
      #print (tree.GenHHCost,cost,bcost)
      effSumV0 = sumHAnalyticalBin.GetBinContent(bmhh,abs(bcost))  # quantity of simulated events in that bin (without cuts)
      weight = model.getScaleFactor(mhh , cost,kl, kt,c2,cg,c2g, effSumV0 , calcSumOfWeights)   # model.effSM,model.MHH,model.COSTS,model.A1,model.A3,model.A7, effSumV0) 
      #############################################
      # fill histograms to test
      #############################################
      if weight > 0: 
               #print countevent
               CalcMhh[countevent] = float(mhhcost[0]) 
               CalcCost[countevent] = float(abs(mhhcost[1])) 
               CalcPtH[countevent] = float(mhhcost[2]) 
               CalcPtHH[countevent] = float(mhhcost[3]) 
               CalcWeight[countevent] = weight 
               countevent+=1
               sumWeight+=weight
  print "plotted histogram reweighted from ",countevent," events, ", float(100*(nev-countevent)/nev)," % of the events was lost in empty bins in SM simulation"
  print "sum of weights",sumWeight
  ############################################################################################################################
  # Draw test histos
  ###############################################
  drawtest =0 
  nevtest=50000
  if kl == 1 and kt == 1 and c2 ==0 and cg == 0 and c2g ==0 : 
     filne = "/eos/user/a/acarvalh/asciiHH_tofit/GF_HH_BSM/GF_HH_0.lhe.decayed"    # 0 is SM
     nevtest = 100000
     drawtest = 1 
  # BSM events
  pathBSMtest="/afs/cern.ch/work/a/acarvalh/generateHH/asciiHH_tofit/GF_HH_toRecursive/" # events of file to superimpose a test
  # see the translation of coefficients for this last on: If you make this script smarter (to only read files we ask to test) you can implement more
  # https://github.com/acarvalh/generateHH/blob/master/fit_GF_HH_lhe/tableToFitA3andA7.txt
  if kl == -10 and kt == 0.5 and c2 ==0 and cg == 0 and c2g ==0 :
     drawtest =1
     filne = pathBSMtest+"GF_HH_42.lhe.decayed"
  if kl == 0.0001 and kt == 2.25 and c2 ==0 and cg == 0 and c2g ==0  :
     drawtest =1
     filne = pathBSMtest+"GF_HH_9.lhe.decayed"
  if kl == 2.5 and kt == 1.0 and c2 ==0 and cg == 0 and c2g ==0  :
     drawtest =1
     filne = pathBSMtest+"GF_HH_60.lhe.decayed"
  klJHEP=[1.0, 7.5,  1.0,  1.0,  -3.5, 1.0, 2.4, 5.0, 15.0, 1.0, 10.0, 2.4, 15.0]
  ktJHEP=[1.0, 1.0,  1.0,  1.0,  1.5,  1.0, 1.0, 1.0, 1.0,  1.0, 1.5,  1.0, 1.0]
  c2JHEP=[0.0, -1.0, 0.5, -1.5, -3.0,  0.0, 0.0, 0.0, 0.0,  1.0, -1.0, 0.0, 1.0]
  cgJHEP=[0.0, 0.0, 0.6,  0.0, 0.0,   0.8, 0.2, 0.2, -1.0, -0.6, 0.0, 1.0, 0.0]
  c2gJHEP=[0.0, 0.0, 1.0, -0.8, 0.0, -1.0, -0.2,-0.2,  1.0,  0.6, 0.0, -1.0, 0.0]
  # python recastHH.py --kl 7.5 --kt 1 --c2 -1
  # python recastHH.py --kl 1.0 --kt 1.0 --c2 0.5 --cg -0.8 --c2g 0.6
  # python recastHH.py --kl 1.0 --kt 1.0 --c2 -1.5 --cg 0.0 --c2g -0.8
  # python recastHH.py --kl 1.0 --kt 1.0 --c2 0.0 --cg 0.8 --c2g -1.0
  # python recastHH.py --kl 1.0 --kt 1.0 --c2 1.0 --cg -0.6 --c2g 0.6
  for sam in range(0,13):
    #print (sam, ktJHEP[sam] , kt , klJHEP[sam] , c2 ,c2JHEP[sam] , cg , cgJHEP[sam] , c2g , c2gJHEP[sam])
    if kl == klJHEP[sam] and kt == ktJHEP[sam] and c2 ==c2JHEP[sam] and cg == cgJHEP[sam] and c2g ==c2gJHEP[sam] : 
       print sam
       filne="/eos/user/a/acarvalh/asciiHH_tofit/GF_HH_BSM/GF_HH_"+str(sam)+".lhe.decayed"
       nevtest=100000
       drawtest = sam 
  ############################################################################################################################
  CalcMhhTest = np.zeros((nevtest))
  CalcCostTest = np.zeros((nevtest))
  CalcPtHTest = np.zeros((nevtest))
  CalcPtHHTest = np.zeros((nevtest))

  CalcPtHggTest = np.zeros((nevtest))
  CalcMhhRecoTest = np.zeros((nevtest))
  CalcMXRecoTest = np.zeros((nevtest))
  if options.doPlot : 
    print "draw plain histogram to test"
    model.LoadTestEvents(CalcMhhTest,CalcCostTest,CalcPtHTest,CalcPtHHTest,filne)  
    model.plotting(kl,kt,c2,cg,c2g,CalcMhh,CalcCost,CalcPtHgg,CalcPtHH,CalcMhhReco,CalcMXReco,CalcWeight,CalcMhhTest,CalcCostTest,CalcPtHggTest,CalcPtHHTest,CalcMhhRecoTest,CalcMXRecoTest,drawtest)

##########################################
if __name__ == "__main__":  
   main()


#print len(MHH),A1[0][0]

#options.kll, options.ktt)


# obtaining BSM/SM scale factors

#canvas = ROOT.TCanvas("")
#scaleFactors.Draw('colz')
#canvas.SaveAs("{}.pdf".format(scaleFactors.GetName()))
