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
		HMdatacard = nodeSeparateLocation+"/HighMass_"+n+"/datacards/hgg.mH125_8TeV.txt"
		LMdatacard = nodeSeparateLocation+"/LowMass_"+n+"/datacards/hgg.mH125_8TeV.txt"
		os.system("combineCards.py "+ HMdatacard + " " + LMdatacard + " > " + folder+"/"+n+"/datacards/hgg.mH125_8TeV.txt")
		combinedCardName = folder+"/"+n+"/datacards/hgg.mH125_8TeV.txt"
		#CombineDatacards_MediumBTag_v44/SeparateCategories//HighMass_Node0/datacards/
		combinedCard = open(combinedCardName, "r")
		lines = []
		for line in combinedCard:
			if str(folder+"/SeparateCategories//HighMass_"+str(n)+"/datacards/") in line:
				lines.append(line.replace(folder+"/SeparateCategories//HighMass_"+str(n)+"/datacards/"+folder, folder+"/SeparateCategories/"))
				continue
			if str(folder+"/SeparateCategories//LowMass_"+str(n)+"/datacards/") in line:
				lines.append(line.replace(folder+"/SeparateCategories//LowMass_"+str(n)+"/datacards/"+folder, folder+"/SeparateCategories/"))
				continue
			else:
				lines.append(line)
		combinedCard.close()
		os.system("rm "+combinedCardName)
		combinedCard = open(combinedCardName, "w")
		for line in lines:
			combinedCard.write(line)
		combinedCard.close()
				

if __name__ == "__main__":
	main(sys.argv[1:])
