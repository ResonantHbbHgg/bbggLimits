import os,sys,json,time,re
import logging
from shutil import copy

def DataCardMaker_wHiggs(Folder, nCats, signalExp, observed, higgsExp, HType):
  if HType=="HighMass":
    inputDatacardName = os.getenv("CMSSW_BASE")+'/src/HiggsAnalysis/bbggLimits/LimitSetting/Models/NonResDatacardModel_HM_wHiggs.txt'
  else:
    inputDatacardName = os.getenv("CMSSW_BASE")+'/src/HiggsAnalysis/bbggLimits/LimitSetting/Models/NonResDatacardModel_LM_wHiggs.txt'

  inputDatacard = open(inputDatacardName, 'r')
  outputDatacard = open(Folder+'/datacards/hhbbgg_13TeV_DataCard.txt', 'w')
  outToWrite = ''
  for line in inputDatacard:
    outTemp = line
    ##workspaces
    outTemp = outTemp.replace("INPUTBKGLOC", Folder + '/workspaces/hhbbgg.inputbkg_13TeV.root')
    outTemp = outTemp.replace("INPUTSIGLOC", Folder + '/workspaces/hhbbgg.mH125_13TeV.inputsig.root')
    ##observed
    outTemp = outTemp.replace("OBSCAT0", '{:.0f}'.format(float(str(observed.split(',')[0]))))
    outTemp = outTemp.replace("OBSCAT1", '{:.0f}'.format(float(str(observed.split(',')[1]))))
    ##expected signal
    outTemp = outTemp.replace("SIGCAT0", str(signalExp.split(',')[0]))
    outTemp = outTemp.replace("SIGCAT1", str(signalExp.split(',')[1]))
    ## higgs
    for hty in higgsExp:
      upper_hty = hty.upper()
      #location
      outTemp = outTemp.replace("INPUT"+upper_hty+"LOC", Folder + '/workspaces/hhbbgg.'+hty+'.inputhig.root')
      #exp
      outTemp = outTemp.replace(upper_hty+"C0", str(higgsExp[hty][0]))
      outTemp = outTemp.replace(upper_hty+"C1", str(higgsExp[hty][1]))
    outToWrite += outTemp
  outputDatacard.write(outToWrite)
  outputDatacard.close()
      

def DataCardMaker(Folder, nCats, signalExp, observed, isRes = 0, HType=0):
  if isRes == 0 and nCats == 1:
    print 'Resonant needs two cats!'
    sys.exit(2)

  if nCats == 2:
    inputDatacardName = os.getenv("CMSSW_BASE")+'/src/HiggsAnalysis/bbggLimits/LimitSetting/Models/LowMassResDatacardModel.txt'
    if isRes == 0 and HType=="HighMass":
      inputDatacardName = os.getenv("CMSSW_BASE")+'/src/HiggsAnalysis/bbggLimits/LimitSetting/Models/NonResDatacardModel_HM.txt'
    if isRes == 0 and HType=="LowMass":
      inputDatacardName = os.getenv("CMSSW_BASE")+'/src/HiggsAnalysis/bbggLimits/LimitSetting/Models/NonResDatacardModel_HM.txt'

    inputDatacard = open(inputDatacardName, 'r')
    outputDatacard = open(Folder+'/datacards/hhbbgg_13TeV_DataCard.txt', 'w')
    outToWrite = ''
    for line in inputDatacard:
      outTemp = line.replace("INPUTBKGLOC", Folder+'/workspaces/hhbbgg.inputbkg_13TeV.root')
      outTemp2 = outTemp.replace("INPUTSIGLOC", Folder+'/workspaces/hhbbgg.mH125_13TeV.inputsig.root')
      outTemp3 = outTemp2.replace("OBSCAT0", '{:.0f}'.format(float(str(observed.split(',')[0]))))
      outTemp4 = outTemp3.replace("OBSCAT1", '{:.0f}'.format(float(str(observed.split(',')[1]))))
      outTemp5 = outTemp4.replace("SIGCAT0", str(signalExp.split(',')[0]))
      outTemp6 = outTemp5.replace("SIGCAT1", str(signalExp.split(',')[1]))
      if float(observed.split(',')[0]) < 11:
        newTemp1 = outTemp6
        newTemp2 = newTemp1.replace('CMS_hhbbgg_13TeV_mgg_bkg_slope3_cat0', '### CMS_hhbbgg_13TeV_mgg_bkg_slope3_cat0')
        outTemp6 = newTemp2.replace('CMS_hhbbgg_13TeV_mjj_bkg_slope3_cat0', '### CMS_hhbbgg_13TeV_mjj_bkg_slope3_cat0')
      if float(observed.split(',')[1]) < 11:
        newTemp1 = outTemp6
        newTemp2 = newTemp1.replace('CMS_hhbbgg_13TeV_mgg_bkg_slope3_cat1', '### CMS_hhbbgg_13TeV_mgg_bkg_slope3_cat1')
        outTemp6 = newTemp2.replace('CMS_hhbbgg_13TeV_mjj_bkg_slope3_cat1', '### CMS_hhbbgg_13TeV_mjj_bkg_slope3_cat1')
      outToWrite += outTemp6
    outputDatacard.write(outToWrite)
    outputDatacard.close()

  if nCats == 1 and isRes == 1:
    inputDatacardName = 'Models/HighMassResDatacardModel.txt'
    inputDatacard = open(inputDatacardName, 'r')
    outputDatacard = open(Folder+'/datacards/hhbbgg_13TeV_DataCard.txt', 'w')
    outToWrite = ''
    for line in inputDatacard:
      outTemp = line.replace("INPUTBKGLOC", Folder+'/workspaces/hhbbgg.inputbkg_13TeV.root')
      outTemp2 = outTemp.replace("INPUTSIGLOC", Folder+'/workspaces/hhbbgg.mH125_13TeV.inputsig.root')
      outTemp3 = outTemp2.replace("OBSCAT0", str(observed))
      outTemp5 = outTemp3.replace("SIGCAT0", str(signalExp))
      if float(observed.split(',')[0]) < 11:
        newTemp1 = outTemp5
        newTemp2 = newTemp1.replace('CMS_hhbbgg_13TeV_mgg_bkg_slope3_cat0', '### CMS_hhbbgg_13TeV_mgg_bkg_slope3_cat0')
        outTemp5 = newTemp2.replace('CMS_hhbbgg_13TeV_mjj_bkg_slope3_cat0', '### CMS_hhbbgg_13TeV_mjj_bkg_slope3_cat0')
      outToWrite += outTemp5
    outputDatacard.write(outToWrite)
    outputDatacard.close()

