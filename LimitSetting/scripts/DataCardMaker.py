import sys, getopt, os

def main(argv):
	Folder = ''
	nCats = -99
	signalExp = ''
	observed = ''
	isRes = 0
	try:
		opts, args = getopt.getopt(argv,"f:n:s:b:o:r",["Folder=", "nCats=", "signalExp=", "observed=", "isRes"])
	except getopt.GetoptError:
		print 'DataCardMaker.py -f <Folder> -n <nCats> -s <sigExpCat0,sigExpCat1,...>  -o <observedCat0,observedCat1,...> [if resonant: -r]'
		sys.exit(2)
	for opt, arg in opts:
		if opt == "-m":
			mass = int(arg)
		if opt == "-f":
			Folder = arg
		if opt == "-n":
			nCats = int(arg)
		if opt == "-s":
			signalExp = arg
		if opt == "-o":
			observed = arg
		if opt == "-r":
			isRes = 1

	if Folder == '' or nCats == -99 or signalExp == '' or observed == '':
		print 'DataCardMaker.py -f <Folder> -n <nCats> -s <sigExpCat0,sigExpCat1,...> -o <observedCat0,observedCat1,...> [if resonant: -r]'
		sys.exit(2)

	if isRes == 1 and nCats == 1:
		print 'Resonant needs two cats!'
		sys.exit(2)	

	if nCats == 2:
		inputDatacardName = os.getenv("CMSSW_BASE")+'/src/HiggsAnalysis/bbggLimits/LimitSetting/Models/LowMassResDatacardModel.txt'
		if isRes == 0:
			inputDatacardName = os.getenv("CMSSW_BASE")+'/src/HiggsAnalysis/bbggLimits/LimitSetting/Models/NonResDatacardModel.txt'
                        
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
				

if __name__ == "__main__":
	main(sys.argv[1:])
