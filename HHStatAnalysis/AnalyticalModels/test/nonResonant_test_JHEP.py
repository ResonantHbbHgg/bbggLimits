#! /usr/bin/env python
# Analytical reweighting implementation for H->4b
# This file is part of https://github.com/cms-hh/HHStatAnalysis.
# python nonResonant_test_JHEP.py  --kl 1 --kt 1
# python nonResonant_test_JHEP.py  --kl 0.0001 --kt 0.0001 --c2 -3.0 --cg 0.0 --c2g -1.5
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

print "Weights calculated from the 12 benchmarks defined in 1507.02245v4 (JHEP version) each one with 100k events "
#if c2 != 0 or cg != 0 or c2g != 0 :  print "The analytical function is not yet implemented"

###########################################################
# read events and apply weight
###########################################################
def main():
  model = NonResonantModel()
  # obtaining BSM/SM coeficients
  #dumb = model.ReadCoefficients("../data/coefficientsByBin_A1A3A7.txt") 
  #dumb = model.ReadCoefficients("../data/coefficientsByBin_extended_3M.txt") 
  dumb = model.ReadCoefficients("../data/coefficientsByBin_extended_3M_costHHSim_19-4.txt") 
  # We sum SM + box + the benchmarks from 2-13 
  # read the 2D histo referent to the sum of events
  """
  fileHH=ROOT.TFile("../../../Support/NonResonant/Distros_5p_SM3M_sumBenchJHEP_13TeV.root")
  sumJHEPAnalyticalBin = fileHH.Get("H1bin2")
  SMAnalyticalBin = fileHH.Get("H0bin2")
  """
  fileHH=ROOT.TFile("../../../Analysis/Support/NonResonant/Distros_5p_SM3M_sumBenchJHEP_13TeV_19-4.root") #Distros_5p_SM3M_sumBenchJHEP_13TeV.root") #Distros_5p_SM3M_rebin_sumBenchJHEP_5D_13TeV.root") #
  #fileHH=ROOT.TFile("../../../Analysis/Support/NonResonant/Distros_5p_SM3M_rebin_sumBenchJHEP_5D_13TeV.root") 
  #sumJHEPAnalyticalBin = fileHH.Get("H1bin3")
  histfile = "../../../Analysis/Support/NonResonant/Distros_5p_SM3M_sumBenchJHEP_13TeV_19-4.root" #Distros_5p_SM3M_sumBenchJHEP_13TeV.root" # #Distros_5p_SM3M_rebin_sumBenchJHEP_5D_13TeV.root"
  #histfile = "../../../Analysis/Support/NonResonant/Distros_5p_SM3M_rebin_sumBenchJHEP_5D_13TeV.root"
  histtitle = "H1bin4"
  sumJHEPAnalyticalBin = fileHH.Get(histtitle)
  SMAnalyticalBin = fileHH.Get("H0bin4")
  
  #fileHHname = "../../../Analysis/Support/NonResonant/Distros_5p_SM3M_rebin_sumBenchJHEP_5D_13TeV.root"
  calcSumOfWeights = model.getNormalization(kl, kt,c2,cg,c2g,histfile,histtitle)  # this input is flexible, tatabb may have only the SM
  #print ("normalization is: ",calcSumOfWeights)
  #xaxis = sumJHEPAnalyticalBin.GetXaxis()
  #yaxis = sumJHEPAnalyticalBin.GetYaxis()
  #print "Sum hist ",sumJHEPAnalyticalBin.GetNbinsX(),sumJHEPAnalyticalBin.GetNbinsY(),sumJHEPAnalyticalBin.Integral(),sumJHEPAnalyticalBin.GetXaxis().GetBinLowEdge(1),sumJHEPAnalyticalBin.GetXaxis().GetBinUpEdge(xaxis.GetNbins())
  #print SMAnalyticalBin.GetBinContent(4,4)
  # now loop over events, calculate weights using the coeffitients and  plot histograms
  # events to reweights, in text format (for testing only)
  pathBenchEvents="/eos/user/a/acarvalh/asciiHH_tofit/GF_HH_BSM/" # events to reweight
  # declare the histograms 
  CalcMhh = np.zeros((1200000))
  CalcCost = np.zeros((1200000))
  CalcPtH = np.zeros((1200000))
  CalcPtHgg = np.zeros((1200000))
  CalcPtHH = np.zeros((1200000))
  CalcWeight = np.zeros((1200000))
  CalcMhhReco = np.zeros((1200000))
  CalcMXReco = np.zeros((1200000))
  ##########################################
  # initialize tables of coefficients by bins
  # calculate mhh and cost* and find the bin
  # initialize events reading 
  countline=0
  countevent=0
  counteventSM=0
  # read events as text files for events to test 
  # particuliarity of the text file with events = each 2 lines are one event there
  # save the information of the two Higgses
  Px = np.zeros((2)) 
  Py = np.zeros((2)) 
  Pz = np.zeros((2)) 
  En = np.zeros((2))  
  for sam in  range(1,13): # read events from the list of 12 benchmarks
       filne = pathBenchEvents+"GF_HH_"+str(sam)+".lhe.decayed"    # 0 is SM = here it does not enter in the events to be reweighted in this version
       f = open(filne, 'r+')
       lines = f.readlines() # get all lines as a list (array)
       countline = 0 # particuliarity of the text file with events = each 2 lines are one event there
       for line in  lines:
          model.ReadLine(line, countline,Px,Py,Pz,En)
          #print countline
          countline+=1
          mhhcost= [0,0,0,0] # to store [mhh , cost,ptH , ptHH] of that event
          if countline==2 : # if read 2 lines 
            model.CalculateMhhCost(mhhcost,countline,Px,Py,Pz,En) # ==> adapt to your input 
            bmhh = sumJHEPAnalyticalBin.GetXaxis().FindBin(mhhcost[0])
            bcost = sumJHEPAnalyticalBin.GetYaxis().FindBin(abs(mhhcost[1]))
            #print (mhhcost[1],bcost)
            effSumV0 = sumJHEPAnalyticalBin.GetBinContent(bmhh,bcost)  # quantity of simulated events in that bin (without cuts)
            #weight = model.getScaleFactor(mhhcost,kl, kt,0,model.effSM,model.effSum,model.MHH,model.COSTS,model.A1,model.A3,model.A7,0)  
            #print mhhcost[1],bcost
            #print effSumV0
            weight = model.getScaleFactor(mhhcost[0],mhhcost[1],kl, kt,c2,cg,c2g, effSumV0,calcSumOfWeights) 
            countline=0
            #############################################
            # fill histograms
            #############################################
            if weight > 0: 
               #print countevent
               CalcMhh[countevent] = float(mhhcost[0]) 
               CalcCost[countevent] = float(mhhcost[1]) 
               CalcPtH[countevent] = float(mhhcost[2]) 
               CalcPtHH[countevent] = float(mhhcost[3]) 
               CalcWeight[countevent] = weight 
               countevent+=1
            #else : print ("neg weight ",weight,mhhcost[0],mhhcost[1],bmhh,bcost,kl, kt,c2,cg,c2g,effSumV0,calcSumOfWeights)
       f.close()
  print "plotted hostogram reweighted from ",countevent," events, ", float(100*(1200000-countevent)/1200000)," % of the events was lost in empty bins in SM simulation"
  print kl, kt,c2,cg,c2g, "Sum of weights:", CalcWeight.sum()," correction ",calcSumOfWeights
  ############################################################################################################################
  # Draw test histos
  ###############################################
  drawtest =-1 
  nevtest=50000
  #if kl == 1 and kt == 1 and c2 ==0 and cg == 0 and c2g ==0 : 
  #   filne = pathBenchEvents+"GF_HH_0.lhe.decayed"    # 0 is SM
  #   nevtest = 100000
  #   drawtest = 1 
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
  if kl == 0.0001 and kt == 1.0 and c2 ==0 and cg == 0 and c2g ==0  :
     drawtest =1
     filne = pathBSMtest+"GF_HH_4.lhe.decayed"
  if kl == 2.5 and kt == 1.0 and c2 ==0 and cg == 0 and c2g ==0  :
     drawtest =1
     filne = pathBSMtest+"GF_HH_60.lhe.decayed"
  if kl == -15 and kt == 0.5 and c2 ==0 and cg == 0 and c2g ==0  :
     drawtest =1
     filne = pathBSMtest+"GF_HH_40.lhe.decayed"
  if kl == 5 and kt == 1.5 and c2 ==0 and cg == 0 and c2g ==0  :
     drawtest =1
     filne = pathBSMtest+"GF_HH_74.lhe.decayed"
  if kl == 7.5 and kt == 2.0 and c2 ==0 and cg == 0 and c2g ==0  :
     drawtest =1
     filne = pathBSMtest+"GF_HH_88.lhe.decayed"
  if kl == 0.0001 and kt == 0.0001 and c2 ==-3.0 and cg == 0.0 and c2g ==-1.5  :
     drawtest =1
     filne = pathBSMtest+"GF_HH_281.lhe.decayed"
  if kl == -10 and kt == 0.0 and c2 ==1.0 and cg == 1.0 and c2g ==1.0  :
     drawtest =1
     filne = pathBSMtest+"GF_HH_280.lhe.decayed"
  #########################################################
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
  if options.doPlot and drawtest >0 : 
     print "draw plain histogram to test"
     model.LoadTestEvents(CalcMhhTest,CalcCostTest,CalcPtHTest,CalcPtHHTest,filne)  
  model.plotting(kl,kt,c2,cg,c2g,CalcMhh,CalcCost,CalcPtHgg,CalcPtHH,CalcMhhReco,CalcMXReco,CalcWeight,CalcMhhTest,CalcCostTest,CalcPtHggTest,CalcPtHHTest,CalcMhhRecoTest,CalcMXRecoTest,drawtest)
###############################################################################################################################
if __name__ == "__main__":  
   main()


#print len(MHH),A1[0][0]

#options.kll, options.ktt)


# obtaining BSM/SM scale factors

#canvas = ROOT.TCanvas("")
#scaleFactors.Draw('colz')
#canvas.SaveAs("{}.pdf".format(scaleFactors.GetName()))
