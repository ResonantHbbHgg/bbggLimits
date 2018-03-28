import os,sys,string

def DataCardMaker_wHiggs(Folder, nCats, signalExp, observed, higgsExp, HType):
  if HType=="HighMass":
    inputDatacardName = os.getenv("CMSSW_BASE")+'/src/HiggsAnalysis/bbggLimits/Models/NonResDatacardModel_HM_wHiggs.txt'
  else:
    inputDatacardName = os.getenv("CMSSW_BASE")+'/src/HiggsAnalysis/bbggLimits/Models/NonResDatacardModel_LM_wHiggs.txt'

  outToWrite=''
  with open(inputDatacardName, 'r') as cardTemp:
    outToWrite = cardTemp.read()

  #print outToWrite
  outToWrite = outToWrite.replace("INPUTBKGLOC", str(Folder + '/workspaces/hhbbgg.inputbkg_13TeV.root'))
  outToWrite = outToWrite.replace("INPUTSIGLOC", str(Folder + '/workspaces/hhbbgg.mH125_13TeV.inputsig.root'))
  ##observed
  #print observed.split(',')
  outToWrite = outToWrite.replace("OBSCAT0", str('%.0f'%float(observed.split(',')[0])))
  outToWrite = outToWrite.replace("OBSCAT1", '%.0f'%float(observed.split(',')[1]))
  #print outToWrite
  ##expected signal
  outToWrite = outToWrite.replace("SIGCAT0", str(signalExp.split(',')[0]))
  outToWrite = outToWrite.replace("SIGCAT1", str(signalExp.split(',')[1]))
  ## higgs
  for hty in higgsExp:
    upper_hty = hty.upper()
    #location
    outToWrite = outToWrite.replace("INPUT"+upper_hty+"LOC", Folder + '/workspaces/hhbbgg.'+hty+'.inputhig.root')
    #exp
    outToWrite = outToWrite.replace(upper_hty+"C0", str(higgsExp[hty][0]))
    outToWrite = outToWrite.replace(upper_hty+"C1", str(higgsExp[hty][1]))
    
  with open(Folder+'/datacards/hhbbgg_13TeV_DataCard.txt', 'w') as outputDatacard:
    outputDatacard.write(outToWrite)
      
def DataCardMaker(Folder, nCats, signalExp, observed, isRes = 0, HType=0):
  if isRes==0:
    print 'This fucntion is to be used for Resonant only. Use DataCardMaker_wHiggs() for Non-resonant'
    sys.exit(2)

  if nCats == 2:
    inputDatacardName = os.getenv("CMSSW_BASE")+'/src/HiggsAnalysis/bbggLimits/Models/LowMassResDatacardModel.txt'

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

