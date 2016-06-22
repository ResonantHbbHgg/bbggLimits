import sys, getopt, os

def main(argv):
	folder = ""
	try:
		opts, args = getopt.getopt(argv,"f:",["folders="])
	except getopt.GetoptError:
		print 'NonResonantOrganizer.py -f <folder>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == "-f":
			folder = arg
	if folder == "":
		print 'NonResonantOrganizer.py -f <folder>'
		sys.exit(2)
	
	prevDir = ""
	nodesList = []
	nodeSeparateLocation = folder+"/SeparateCategories/"
	if not os.path.isdir(nodeSeparateLocation):
		prevDir = os.listdir(folder)
		os.mkdir(folder+"/SeparateCategories/")
		for f in prevDir:
			if ".json" in f:
				continue
			os.rename(folder+"/"+f, folder+"/SeparateCategories/"+f)
			node = f.replace("HighMass_", "").replace("LowMass_", "")
			if node not in nodesList:
				nodesList.append(node)
		
	if os.path.isdir(nodeSeparateLocation):
		prevDir = os.listdir(folder+"/SeparateCategories/")
		for f in prevDir:
			if ".json" in f:
				continue
			node = f.replace("HighMass_", "").replace("LowMass_", "")
			if node not in nodesList:
				nodesList.append(node)

	if len(nodesList) < 1:
		print "EMPTY LIST OF NODES!!"
		sys.exit(2)
	
	for n in nodesList:
		if not os.path.isdir(folder+"/"+n): os.mkdir(folder+"/"+n)
		if not os.path.isdir(folder+"/"+n+"/datacards/"): os.mkdir(folder+"/"+n+"/datacards/")
		HMdatacard = open(nodeSeparateLocation+"/HighMass_"+n+"/datacards/hgg.mH125_8TeV.txt", "r")
		HMdatacardLines = []
		for line in HMdatacard:
			HMdatacardLines.append(line)
		LMdatacard = open(nodeSeparateLocation+"/LowMass_"+n+"/datacards/hgg.mH125_8TeV.txt", "r")
		LMdatacardLines = []
		for line in LMdatacard:
			LMdatacardLines.append(line)
		combinedDatacard = open(folder+"/"+n+"/datacards/hgg.mH125_8TeV.txt", "w")
		if len(LMdatacardLines) != len(HMdatacardLines):
			print "SOMETHING's WRONG! DATACARDS HAVE DIFFERENT NUMBER OF LINES", nodeSeparateLocation/"HighMass_"+n+"/datacards/hgg.mH125_8TeV.txt", nodeSeparateLocation/"LowMass_"+n+"/datacards/hgg.mH125_8TeV.txt"
			sys.exit(2)
		for l in range(0, len(LMdatacardLines)):
			if "imax" in LMdatacardLines[l]:
				combinedDatacard.write("imax 4\n")
				continue
			if "jmax" in LMdatacardLines[l]:
				combinedDatacard.write("jmax *\n")
				continue
			if "shapes" in LMdatacardLines[l]:
				combinedDatacard.write(LMdatacardLines[l].replace("LowMass", "SeparateCategories/LowMass"))
				HMNewCat2 = HMdatacardLines[l].replace(" cat0 ", " cat2 ")
				HMNewCat3 = HMNewCat2.replace(" cat1 ", " cat3 ")
				combinedDatacard.write(HMNewCat3.replace("HighMass", "SeparateCategories/HighMass"))
				continue
			if "bin cat0 cat0" in LMdatacardLines[l]:
				combinedDatacard.write("bin cat0 cat0 cat1 cat1 cat2 cat2 cat3 cat3\n")
				continue
			if "bin cat0 cat1" in LMdatacardLines[l]:
				combinedDatacard.write("bin cat0 cat1 cat2 cat3\n")
				continue
			if "observation" in LMdatacardLines[l]:
				toWrite = "observation "
				toWrite += LMdatacardLines[l].replace("observation ", "").replace("\n", "")
				toWrite += HMdatacardLines[l].replace("observation", "")
				combinedDatacard.write(toWrite)
				continue
			if "process Sig Bkg" in LMdatacardLines[l]:
				combinedDatacard.write("process Sig Bkg Sig Bkg Sig Bkg Sig Bkg\n")
				continue
			if "process 0 1" in LMdatacardLines[l]:
				combinedDatacard.write("process 0 1 0 1 0 1 0 1\n")
				continue
			if "rate" in LMdatacardLines[l]:
				toWrite = "rate "
				toWrite += LMdatacardLines[l].replace("rate ", "").replace("\n", "")
				toWrite += HMdatacardLines[l].replace("rate", "")
				combinedDatacard.write(toWrite)
				continue
			if "CMS_" in LMdatacardLines[l] and "cat0" in LMdatacardLines[l]:
				combinedDatacard.write(LMdatacardLines[l])
				toWrite = LMdatacardLines[l].replace("cat0", "cat2")
				combinedDatacard.write(toWrite)
				continue
			if "CMS_" in LMdatacardLines[l] and "cat1" in LMdatacardLines[l]:
				combinedDatacard.write(LMdatacardLines[l])
				toWrite = LMdatacardLines[l].replace("cat1", "cat3")
				combinedDatacard.write(toWrite)
				continue
			systs = ["lumi_13TeV", "DiphoTrigger", "CMS_hgg_eff_g", "btag_eff", "maajj_cut_acceptance", "PDF", "QCD_scale", "gg_migration", "gluonSplitting", "pTj_cut_acceptance"]
			skip = 0
			for sys in systs:
				if sys in LMdatacardLines[l]:
					toWrite = LMdatacardLines[l].replace("\n", "")
					toWrite += HMdatacardLines[l].replace(sys, "").replace("lnN", "")
					combinedDatacard.write(toWrite)
					skip = 1
			if skip: continue
			if LMdatacardLines[l] == HMdatacardLines[l]:
				combinedDatacard.write(LMdatacardLines[l])
				continue
			
		combinedDatacard.close()
		LMdatacard.close()
		HMdatacard.close()
				

if __name__ == "__main__":
	main(sys.argv[1:])
