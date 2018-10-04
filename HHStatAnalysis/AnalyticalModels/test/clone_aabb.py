#!/usr/bin/env python

histdone=1 # 0 = no
# option 0 is  to make the 2D histos to normalize, and a second one with 1 to make the file with events and weights to test

# good old python modules
import json
import os
import importlib
from array import array
from glob import glob
import numpy as np
import os, sys, time,math
import shutil,subprocess
from HHStatAnalysis.AnalyticalModels.NonResonantModel import NonResonantModel

# ROOT imports
import ROOT
from ROOT import TChain, TH1F, TFile, vector, gROOT, TTree

##################
# read the histos tree and contruct the tree of the relevant variables 
LM=1
if LM == 0 : path = "/afs/cern.ch/user/a/andrey/public/HH/LT-Jan25-APZ_HighMass/" 
if LM == 1 : path = "/afs/cern.ch/user/a/andrey/public/HH/LT-Jan25-APZ_LowMass/" 
if histdone==0 : path = "/afs/cern.ch/user/a/andrey/public/HH/GenTrees_for_Xandra/" 
outpath="/afs/cern.ch/user/a/acarvalh/public/toAABB/"
data="../../../HHStatAnalysis/AnalyticalModels/data/"
model = NonResonantModel()
# obtaining BSM/SM coeficients
dumb = model.ReadCoefficients("../../../HHStatAnalysis/AnalyticalModels/data/coefficientsByBin_klkt.txt")  


if LM == 0 : fileoutput = "events_SumV0_afterBaseline_HM.root"
if LM == 1 : fileoutput = "events_SumV0_afterBaseline_LM.root"

files = []
endfile = "_13TeV-madgraph.root"

if histdone==1 : 
  files.append("LT_output_GluGluToHHTo2B2G_node_SM")
  files.append("LT_output_GluGluToHHTo2B2G_node_box")
  for ifile in range(2,14) : files.append("LT_output_GluGluToHHTo2B2G_node_"+str(ifile))

if histdone==0 : 
  files.append("output_GluGluToHHTo2B2G_node_SM")
  files.append("output_GluGluToHHTo2B2G_node_box")
  for ifile in range(2,14) : files.append("output_GluGluToHHTo2B2G_node_"+str(ifile))

# to analytical re
# We sum SM + box + the benchmarks from 2-13 
# read the 2D histo referent to the sum of events
fileHH=ROOT.TFile("../../../Support/NonResonant/Hist2DSum_V0_SM_box.root")
sumHAnalyticalBin = fileHH.Get("SumV0_AnalyticalBin")
sumHBenchBin = fileHH.Get("SumV0_AnalyticalBin")

countSig=np.zeros((15)) 
normSig= np.ones((15)) 
if 1 > 0 :
  ################
  sumWSM=0
  sumWSMA=0
  sumWboxA=0
  sumW1=0
  sumW2=0
  sumW3=0
  sumW4=0
  sumW5=0
  sumW6=0
  sumW7=0
  sumW8=0
  sumW9=0
  sumW10=0
  sumW11=0
  sumW12=0
  # lambda scan: 
  sumL0=0
  sumL2p4=0
  sumL3=0
  sumL4=0
  sumL5=0
  sumL7=0
  sumL10=0
  sumL12p5=0
  sumL15=0
  sumL20=0
  sumLm1=0
  sumLm2p4=0
  sumLm3=0
  sumLm4=0
  sumLm5=0
  sumLm7=0
  sumLm10=0
  sumLm12p5=0
  sumLm15=0
  sumLm20=0
  #fileHH=ROOT.TFile(outpath+"HistSum2D.root")
  sumHBenchBin = fileHH.Get("SumV0_BenchBin")
  sumHAnalyticalBin = fileHH.Get("SumV0_AnalyticalBin")
  print "Sum to Bench hist ",sumHBenchBin.GetNbinsX(),sumHBenchBin.GetNbinsY(),sumHBenchBin.Integral()
  # check mhh hist to SM and box
  binsxV0 = array( 'f',  [240.,250.,270.,290.,310.,330.,350.,390.,410.,450.,500.,550.,600.,700.,800.,1000.,1200 ])
  binsyV0 = array( 'f', [-1., -0.6,0.6,1.] )
  histmhhRe = ROOT.TH1D('MhhRe', '', 30,200.,1000.) 
  histmhhSM = ROOT.TH1D('MhhSM', '', 30,200.,1000.) 
  histmhhReA = ROOT.TH1D('MhhSMA', '', 30,200.,1000.) 

  histmhhRecoRe = ROOT.TH1D('MhhRecoRe', '', 30,200.,1000.) 
  histmhhRecoSM = ROOT.TH1D('MhhRecoSM', '', 30,200.,1000.) 
  histmhhRecoReA = ROOT.TH1D('MhhRecoSMA', '', 30,200.,1000.) 

  histmhhBox = ROOT.TH1D('MhhBox', '', 30,200.,1000.) 
  histmhhBoxA = ROOT.TH1D('MhhBoxA', '', 30,200.,1000.) 
  histmhhL15 = ROOT.TH1D('MhhL15', '', 30,200.,1000.) 
  histmhhL2p4 = ROOT.TH1D('MhhL2p4', '', 30,200.,1000.) 
  histmhhBench1 = ROOT.TH1D('MhhBench1', '', 30,0.,1800.) 
  histmhhBench2 = ROOT.TH1D('MhhBench2', '', 30,0.,1800.) 
  histmhhBench6 = ROOT.TH1D('MhhBench6', '', 30,0.,1800.) 
  ###############################################
  # Read histograms with JHEP benchmarks and SM
  fileH=ROOT.TFile("../../../Support/NonResonant/Distros_5p_500000ev_12sam_13TeV_JHEP_500K.root")
  bench = []
  for ibench in range(0,12) : bench.append(fileH.Get(str(ibench)+"_bin1")) # in the old binning
  print "Bench hist ",bench[0].GetNbinsX(),bench[0].GetNbinsY(),bench[0].Integral()
  fileSM=ROOT.TFile("../../../Support/NonResonant/Distros_5p_SM3M_sumBenchJHEP_13TeV.root")
  histSM = fileSM.Get("H0bin1") # fine binning (the H0bin2 is with the bin used to analytical)
  ##############################################
  # do one file with events to test
  fileout=ROOT.TFile(outpath+fileoutput,"recreate")
  treeout = TTree('treeout', 'treeout')
  Genmhh = np.zeros(1, dtype=float)
  mhh = np.zeros(1, dtype=float)
  mbb = np.zeros(1, dtype=float)
  GenHHCost = np.zeros(1, dtype=float) 
  # weight benchmarks JHEP
  effSumV0AnalyticalBin = np.zeros(1, dtype=float)
  weightSM = np.zeros(1, dtype=float)
  weightSMA = np.zeros(1, dtype=float)
  weight1 = np.zeros(1, dtype=float)
  weight2 = np.zeros(1, dtype=float)
  weight3 = np.zeros(1, dtype=float)
  weight4 = np.zeros(1, dtype=float)
  weight5 = np.zeros(1, dtype=float)
  weight6 = np.zeros(1, dtype=float)
  weight7 = np.zeros(1, dtype=float)
  weight8 = np.zeros(1, dtype=float)
  weight9 = np.zeros(1, dtype=float)
  weight10 = np.zeros(1, dtype=float)
  weight11 = np.zeros(1, dtype=float)
  weight12 = np.zeros(1, dtype=float)
  # lamda scan
  weightL0 = np.zeros(1, dtype=float)
  weightL0p5 = np.zeros(1, dtype=float)
  weightL2 = np.zeros(1, dtype=float)
  weightL2p4 = np.zeros(1, dtype=float) 
  weightL3 = np.zeros(1, dtype=float) 
  weightL4 = np.zeros(1, dtype=float) 
  weightL5 = np.zeros(1, dtype=float) 
  weightL7 = np.zeros(1, dtype=float) 
  weightL10 = np.zeros(1, dtype=float) 
  weightL12p5 = np.zeros(1, dtype=float) 
  weightL15 = np.zeros(1, dtype=float) 
  weightL20 = np.zeros(1, dtype=float) 
  weightLm1 = np.zeros(1, dtype=float) 
  weightLm2p4 = np.zeros(1, dtype=float) 
  weightLm3 = np.zeros(1, dtype=float) 
  weightLm4 = np.zeros(1, dtype=float) 
  weightLm5 = np.zeros(1, dtype=float) 
  weightLm7 = np.zeros(1, dtype=float) 
  weightLm10 = np.zeros(1, dtype=float) 
  weightLm12p5 = np.zeros(1, dtype=float) 
  weightLm15 = np.zeros(1, dtype=float) 
  weightLm20 = np.zeros(1, dtype=float) 
  sumL0=0
  sumL0p5=0
  sumL2=0
  sumL2p4=0
  sumL3=0
  sumL4=0
  sumL5=0
  sumL7=0
  sumL10=0
  sumL12p5=0
  sumL15=0
  sumL20=0
  sumLm1=0
  sumLm2p4=0
  sumLm3=0
  sumLm4=0
  sumLm5=0
  sumLm7=0
  sumLm10=0
  sumLm12p5=0
  sumLm15=0
  sumLm20=0
  treeout.Branch('Genmhh', Genmhh, 'Genmhh/D')
  treeout.Branch('mhh', mhh, 'mhh/D')
  treeout.Branch('mbb', mbb, 'mbb/D')
  treeout.Branch('GenHHCost', GenHHCost, 'GenHHCost/D')
  treeout.Branch('effSumV0AnalyticalBin', effSumV0AnalyticalBin, 'effSumV0AnalyticalBin/D')
  treeout.Branch('weightSM', weightSM, 'weightSM/D')
  treeout.Branch('weightSMA', weightSMA, 'weightSMA/D')
  treeout.Branch('weight1', weight1, 'weight1/D')
  treeout.Branch('weight2', weight2, 'weight2/D')
  treeout.Branch('weight3', weight3, 'weight3/D')
  treeout.Branch('weight4', weight4, 'weight4/D')
  treeout.Branch('weight5', weight5, 'weight5/D')
  treeout.Branch('weight6', weight6, 'weight6/D')
  treeout.Branch('weight7', weight7, 'weight7/D')
  treeout.Branch('weight8', weight8, 'weight8/D')
  treeout.Branch('weight9', weight9, 'weight9/D')
  treeout.Branch('weight10', weight10, 'weight10/D')
  treeout.Branch('weight11', weight11, 'weight11/D')
  treeout.Branch('weight12', weight12, 'weight12/D')
  # lamda scan
  treeout.Branch('weightL0', weightL0, 'weightL0/D')
  treeout.Branch('weightL0p5', weightL0p5, 'weightL0p5/D')
  treeout.Branch('weightL2', weightL2, 'weightL2/D')
  treeout.Branch('weightL2p4', weightL2p4, 'weightL2p4/D')
  treeout.Branch('weightL3', weightL3, 'weightL3/D')
  treeout.Branch('weightL4', weightL4, 'weightL4/D')
  treeout.Branch('weightL5', weightL5, 'weightL5/D')
  treeout.Branch('weightL7', weightL7, 'weightL7/D')
  treeout.Branch('weightL10', weightL10, 'weightL10/D')
  treeout.Branch('weightL12p5', weightL12p5, 'weightL12p5/D')
  treeout.Branch('weightL15', weightL15, 'weightL15/D')
  treeout.Branch('weightL20', weightL20, 'weightL20/D')
  treeout.Branch('weightLm1', weightLm1, 'weightLm1/D')
  treeout.Branch('weightLm2p4', weightLm2p4, 'weightLm2p4/D')
  treeout.Branch('weightLm3', weightLm3, 'weightLm3/D')
  treeout.Branch('weightLm4', weightLm4, 'weightLm4/D')
  treeout.Branch('weightLm5', weightLm5, 'weightLm5/D')
  treeout.Branch('weightLm7', weightLm7, 'weightLm7/D')
  treeout.Branch('weightLm10', weightLm10, 'weightLm10/D')
  treeout.Branch('weightLm12p5', weightLm12p5, 'weightLm12p5/D')
  treeout.Branch('weightLm15', weightLm15, 'weightLm15/D')
  treeout.Branch('weightLm20', weightLm20, 'weightLm20/D')
#########################################################
normBench = [8425.67378246,8465.07111168,8447.02627838,96.8417106507,8443.43311691,8440.44181998,8445.16043015,8440.88468543,8447.01718131,8430.17825354,8444.32422433,8451.44169753]
normSM = 50711.7506188
normSManal =  1
normBox = 1
listLam=[0.0001,0.5,2, 2.5, 3.0, 4.0, 4.0, 5.0, 7.0, 10.0, 12.5, 15.0, 20.0, -1.0,-2.4, -3.0, -4.0, -5.0, -7.0, -10.0, -12.5, -15.0, -20.0]
normL =[] 
for ii in range(0,len(listLam)): normL.append(model.getNormalization(float(listLam[ii]),1.0,sumHAnalyticalBin)/6)
# the factor 6 is from the ratio of the nevents used to normalize to SM
# 6 = 300k / 50k
##############################################################
# loop in all events
countevent=0
for ifile in range(0,14) : # len(files)  
  print path+files[ifile]+endfile
  file=ROOT.TFile(path+files[ifile]+endfile)
  if histdone==1 : tree=file.Get("TCVARS")
  if histdone==0 : tree=file.fsDir.Get("GenTree")
  nev = tree.GetEntries()
  print nev
  counter=0
  #if ifile==16 : treeout.Fill()
  for iev in range(0,nev) :    
    tree.GetEntry(iev)

    #countSig[ifile]+=(w_oneInvFb)/normSig[ifile]
    # make tree for test
    if 1 > 0 :
      if histdone==0 : 
        Genmhh[0] = tree.mHH
        GenHHCost[0] = tree.cosTheta
      elif histdone==1 : 
        Genmhh[0] = tree.gen_mHH
        GenHHCost[0] = tree.gen_cosTheta
      # find the bin the event belong
      bmhh = histSM.GetXaxis().FindBin(Genmhh[0])
      bcost = histSM.GetYaxis().FindBin(GenHHCost[0])
      mergecostSum = 0
      # to the benchmarks sum reduce the binning in cosTheta*
      for ii in range(1,11) : mergecostSum+= sumHBenchBin.GetBinContent(bmhh,ii) 
      if mergecostSum >0 : 
         weightSM[0] = (histSM.GetBinContent(bmhh,bcost) / mergecostSum)/normSM 
         weight1[0] = (bench[0].GetBinContent(bmhh,bcost) / mergecostSum)/normBench[0] 
         weight2[0] = (bench[1].GetBinContent(bmhh,bcost) / mergecostSum)/normBench[1]
         weight3[0] = (bench[2].GetBinContent(bmhh,bcost) / mergecostSum)/normBench[2]   
         weight4[0] = (bench[3].GetBinContent(bmhh,bcost) / mergecostSum)/normBench[3]   
         weight5[0] = (bench[4].GetBinContent(bmhh,bcost) / mergecostSum)/normBench[4] 
         weight6[0] = (bench[5].GetBinContent(bmhh,bcost) / mergecostSum)/normBench[5]   
         weight7[0] = (bench[6].GetBinContent(bmhh,bcost) / mergecostSum)/normBench[6]  
         weight8[0] = (bench[7].GetBinContent(bmhh,bcost) / mergecostSum)/normBench[7]   
         weight9[0] = (bench[8].GetBinContent(bmhh,bcost) / mergecostSum)/normBench[8]  
         weight10[0] = (bench[9].GetBinContent(bmhh,bcost) / mergecostSum)/normBench[9]  
         weight11[0] = (bench[10].GetBinContent(bmhh,bcost) / mergecostSum)/normBench[10] 
         weight12[0] = (bench[11].GetBinContent(bmhh,bcost) / mergecostSum)/normBench[11]    
         sumWSM+=weightSM[0]
         sumW1+=weight1[0]
         sumW2+=weight2[0]
         sumW3+=weight3[0]
         sumW4+=weight4[0]
         sumW5+=weight5[0]
         sumW6+=weight6[0]
         sumW7+=weight7[0]
         sumW8+=weight8[0]
         sumW9+=weight9[0]
         sumW10+=weight10[0]
         sumW11+=weight11[0]
         sumW12+=weight12[0]
      # weight analytical 
      mhhcost= [Genmhh[0],GenHHCost[0]] # to store [mhh , cost] of that event
      bmhh = sumHAnalyticalBin.GetXaxis().FindBin(mhhcost[0])
      bcost = sumHAnalyticalBin.GetYaxis().FindBin(mhhcost[1])
      if sumHAnalyticalBin.GetBinContent(bmhh,bcost) >0 : # to be done with all events
         # find the Nevents from the sum of events on that bin
         effSumV0 = sumHAnalyticalBin.GetBinContent(bmhh,bcost)  # quantity of simulated events in that bin (without cuts)
         weightSMA[0] = model.getScaleFactor(mhhcost,1.0, 1.0, effSumV0)/ (model.getNormalization(1.0, 1.0,sumHAnalyticalBin)/6) #normSManal
         sumWSMA+=weightSMA[0]
         # lambda scan
         weightL0[0] = model.getScaleFactor(mhhcost, 0.0001,1.0, effSumV0) / normL[0]
         weightL0p5[0] = model.getScaleFactor(mhhcost, 0.5,1.0, effSumV0) / normL[1]
         weightL2[0] = model.getScaleFactor(mhhcost, 2,1.0, effSumV0) / normL[2]
         weightL2p4[0] = model.getScaleFactor(mhhcost, 2.5,1.0, effSumV0) / normL[3]
         weightL3[0] = model.getScaleFactor(mhhcost, 3.0,1.0, effSumV0) / normL[4]
         weightL4[0] = model.getScaleFactor(mhhcost, 4.0,1.0, effSumV0) / normL[5]
         weightL5[0] = model.getScaleFactor(mhhcost, 5.0,1.0, effSumV0) / normL[6]
         weightL7[0] = model.getScaleFactor(mhhcost, 7.0,1.0, effSumV0) / normL[7]
         weightL10[0] = model.getScaleFactor(mhhcost, 10.0,1.0, effSumV0) / normL[8]
         weightL12p5[0] = model.getScaleFactor(mhhcost, 12.5,1.0, effSumV0) / normL[9]
         weightL15[0] = model.getScaleFactor(mhhcost, 15.0,1.0, effSumV0) / normL[10]
         weightL20[0] = model.getScaleFactor(mhhcost, 20.0,1.0, effSumV0) / normL[11]
         weightLm1[0] = model.getScaleFactor(mhhcost, -1.0,1.0, effSumV0) / normL[12]
         weightLm2p4[0] = model.getScaleFactor(mhhcost,-2.4,1.0, effSumV0) / normL[13] 
         weightLm3[0] = model.getScaleFactor(mhhcost, -3.0,1.0, effSumV0) / normL[14]
         weightLm4[0] = model.getScaleFactor(mhhcost, -4.0,1.0, effSumV0) / normL[15]
         weightLm5[0] = model.getScaleFactor(mhhcost, -5.0,1.0, effSumV0) / normL[16]
         weightLm7[0] = model.getScaleFactor(mhhcost, -7.0,1.0, effSumV0) / normL[17]
         weightLm10[0] = model.getScaleFactor(mhhcost, -10.0,1.0, effSumV0) / normL[18]
         weightLm12p5[0] = model.getScaleFactor(mhhcost, -12.5,1.0, effSumV0) / normL[19]
         weightLm15[0] = model.getScaleFactor(mhhcost, -15.0,1.0, effSumV0) / normL[20]
         weightLm20[0] = model.getScaleFactor(mhhcost, -20.0,1.0, effSumV0) / normL[21]
         sumL0+=weightL0[0]
         sumL0p5+=weightL0p5[0]
         sumL2+=weightL2[0]
         sumL2p4+=weightL2p4[0]
         sumL3+=weightL3[0]
         sumL4+=weightL4[0]
         sumL5+=weightL5[0]
         sumL7+=weightL7[0]
         sumL10+=weightL10[0]
         sumL12p5+=weightL12p5[0]
         sumL15+=weightL15[0]
         sumL20+=weightL20[0]
         sumLm1+=weightLm1[0]
         sumLm2p4+=weightLm2p4[0]
         sumLm3+=weightLm3[0]
         sumLm4+=weightLm4[0]
         sumLm5+=weightLm5[0]
         sumLm7+=weightLm7[0]
         sumLm10+=weightLm10[0]
         sumLm12p5+=weightLm12p5[0]
         sumLm15+=weightLm15[0]
         sumLm20+=weightLm20[0]
      treeout.Fill()
      # make histogram to test
      histmhhReA.Fill(Genmhh[0],weightSMA[0])
      histmhhRe.Fill(Genmhh[0],weightSM[0])

      histmhhRecoReA.Fill(Genmhh[0],weightSMA[0])
      histmhhRecoRe.Fill(Genmhh[0],weightSM[0])

      histmhhBoxA.Fill(Genmhh[0],weightL0[0])
      histmhhBench1.Fill(Genmhh[0],weight1[0])
      histmhhBench2.Fill(Genmhh[0],weight2[0])
      histmhhBench6.Fill(Genmhh[0],weight6[0])

      histmhhL15.Fill(Genmhh[0],weightL15[0])
      histmhhL2p4.Fill(Genmhh[0],weightL2p4[0])
      if ifile ==0 : 
        histmhhSM.Fill(Genmhh[0])
        histmhhRecoSM.Fill(Genmhh[0])
      if ifile ==1 : histmhhBox.Fill(Genmhh[0])
      countevent+=1
  print counter
print countevent
# save tree of events
if 1 > 0 :
  print "W1",sumW1
  print "W2",sumW2
  print "W3",sumW3
  print "W4",sumW4
  print "W5",sumW5
  print "W6",sumW6
  print "W7",sumW7
  print "W8",sumW8
  print "W9",sumW9
  print "W10",sumW10
  print "W11",sumW11
  print "W12",sumW12
  print "WSM",sumWSM
  print "WSMA",sumWSMA
  print "lambda scan:" 
  print 0,sumL0
  print 0.5,sumL0p5
  print 2,sumL2
  print 2.5,sumL2p4
  print 3,sumL3
  print 4,sumL4
  print 5,sumL5
  print 6,sumL7
  print 10,sumL10
  print 12.5,sumL12p5
  print 15,sumL15
  print 20,sumL20
  print -1,sumLm1
  print -2.4,sumLm2p4
  print -3,sumLm3
  print -4,sumLm4
  print -5,sumLm5
  print -7,sumLm7
  print -10,sumLm10
  print -12.5,sumLm12p5
  print -15,sumLm15
  print -20,sumLm20
  fileout.Write()
  fileout.Close()
  cs=ROOT.TCanvas("cs","cs",10,10,500,500)
  leg = ROOT.TLegend(0.5,0.60,0.99,0.99);
  histmhhRe.Scale(1./histmhhRe.Integral())
  histmhhRe.SetLineWidth(2)
  histmhhRe.SetLineColor(1)
  histmhhRe.Draw("hist")
  if LM == 1 : leg.AddEntry(histmhhRe,"reweigted (from 210k events, hist)")
  if LM == 0 : leg.AddEntry(histmhhRe,"reweigted (from 74k events, hist)")
  histmhhReA.Scale(1./histmhhReA.Integral())
  histmhhReA.SetLineWidth(2)
  histmhhReA.SetLineColor(2)
  histmhhReA.Draw("same,hist")
  if LM == 1 : leg.AddEntry(histmhhReA,"reweigted (from 210k events, from analitical)")
  if LM == 0 : leg.AddEntry(histmhhReA,"reweigted (from 74k events, from analitical)")
  histmhhSM.Scale(1./histmhhSM.Integral())
  histmhhSM.SetLineWidth(2)
  histmhhSM.SetLineColor(8)
  histmhhSM.Draw("same,hist")
  if LM == 1 : leg.AddEntry(histmhhSM,"SM (from ratio 15028 events)")
  if LM == 0 : leg.AddEntry(histmhhSM,"SM (from ratio 2744 events)")
  leg.Draw("same")
  if LM == 1 : cs.SaveAs("SMtest_afterCuts_aabb_HM.png") 
  if LM == 0 : cs.SaveAs("SMtest_afterCuts_aabb_LM.png")
  cs.Clear()
  leg.Clear()
  ###########################
  histmhhRecoRe.Scale(1./histmhhRecoRe.Integral())
  histmhhRecoRe.SetLineWidth(2)
  histmhhRecoRe.SetLineColor(1)
  histmhhRecoRe.Draw("hist")
  leg.AddEntry(histmhhRecoRe,"reweigted (from hist)")
  histmhhRecoReA.Scale(1./histmhhRecoReA.Integral())
  histmhhRecoReA.SetLineWidth(2)
  histmhhRecoReA.SetLineColor(2)
  histmhhRecoReA.Draw("same,hist")
  leg.AddEntry(histmhhRecoReA,"reweigted (from analitical)")
  histmhhRecoSM.Scale(1./histmhhRecoSM.Integral())
  histmhhRecoSM.SetLineWidth(2)
  histmhhRecoSM.SetLineColor(8)
  histmhhRecoSM.Draw("same,hist")
  if LM == 1 : leg.AddEntry(histmhhRecoSM,"SM (from ratio 15028 events)")
  if LM == 0 : leg.AddEntry(histmhhRecoSM,"SM (from ratio 2744 events)")
  leg.Draw("same")
  if LM == 1 : cs.SaveAs("SMtestReco_afterCuts_aabb_HM.png") 
  if LM == 0 : cs.SaveAs("SMtestReco_afterCuts_aabb_LM.png") 
  cs.Clear()
  leg.Clear()
  ###########################
  histmhhBox.Scale(1./histmhhBox.Integral())
  histmhhBox.SetLineWidth(2)
  histmhhBox.SetLineColor(8)
  histmhhBox.Draw("hist")
  histmhhBoxA.Scale(1./histmhhBoxA.Integral())
  histmhhBoxA.SetLineWidth(2)
  histmhhBoxA.SetLineColor(1)
  histmhhBoxA.Draw("same,hist")
  if LM == 1 : leg.AddEntry(histmhhBoxA,"#kappa_{#lambda} = 0 (from 15028 events)")
  if LM == 0 : leg.AddEntry(histmhhBoxA,"#kappa_{#lambda} = 0 (from 4503 events)")
  leg.AddEntry(histmhhBox,"reweigted (from analitical)")
  leg.Draw("same")
  if LM == 0 : cs.SaveAs("Boxtest_afterCuts_aabb_HM.png") 
  if LM == 1 : cs.SaveAs("Boxtest_afterCuts_aabb_LM.png") 
  cs.Clear()
  leg.Clear()
  #################################################"
  histmhhBench1.Scale(1./histmhhBench1.Integral())
  histmhhBench1.SetLineWidth(2)
  histmhhBench1.SetLineColor(1)
  leg.AddEntry(histmhhBench1,"BM1")
  histmhhBench1.Draw("hist")
  histmhhBench2.Scale(1./histmhhBench2.Integral())
  histmhhBench2.SetLineWidth(2)
  histmhhBench2.SetLineColor(8)
  leg.AddEntry(histmhhBench2,"BM2")
  histmhhBench2.Draw("same,hist")
  histmhhBench6.Scale(1./histmhhBench6.Integral())
  histmhhBench6.SetLineWidth(2)
  histmhhBench6.SetLineColor(2)
  leg.AddEntry(histmhhBench6,"BM6")
  histmhhBench6.Draw("same,hist")
  leg.Draw("same")
  if LM == 0 : cs.SaveAs("Bench1_afterCuts_reco_aabb_HM.png") 
  if LM == 1 : cs.SaveAs("Bench1_afterCuts_reco_aabb_LM.png") 
  cs.Clear()
  leg.Clear()
  leg.Clear()
  #################################################"
  histmhhL15.Scale(1./histmhhL15.Integral())
  histmhhL15.SetLineWidth(2)
  histmhhL15.SetLineColor(8)
  histmhhL15.Draw("hist")
  leg.AddEntry(histmhhL15,"#kappa_{#lambda} = 15")
  histmhhBoxA.Scale(1./histmhhBoxA.Integral())
  histmhhBoxA.SetLineWidth(2)
  histmhhBoxA.SetLineColor(2)
  histmhhBoxA.Draw("same,hist")
  leg.AddEntry(histmhhBoxA,"#kappa_{#lambda} = 0")
  histmhhL2p4.Scale(1./histmhhL2p4.Integral())
  histmhhL2p4.SetLineWidth(2)
  histmhhL2p4.SetLineColor(6)
  histmhhL2p4.Draw("same,hist")
  leg.AddEntry(histmhhL2p4,"#kappa_{#lambda} = 2.4")
  histmhhReA.Scale(1./histmhhReA.Integral())
  histmhhReA.SetLineWidth(2)
  histmhhReA.SetLineColor(1)
  histmhhReA.Draw("same,hist")
  leg.AddEntry(histmhhReA,"SM")
  leg.Draw("same")
  if LM == 0 : cs.SaveAs("kl15_afterCuts_reco_aabb_HM.png") 
  if LM == 1 : cs.SaveAs("kl15_afterCuts_reco_aabb_LM.png") 
print "done "
print "Sig ",countSig
