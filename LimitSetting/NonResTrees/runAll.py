from ROOT import *
import os


nodes = [ ["box", 49600], ["SM", 50000], [2, 50000], [3, 47600], [4, 50000], [5, 50000], [6, 50000], [7, 49800], [8, 50000], [9, 50000], [10, 50000], [11, 50000], [12, 50000], [13, 50000] ]

Signals = "/afs/cern.ch/work/r/rateixei/work/DiHiggs/flg76X/CMSSW_7_6_3/src/flashgg/bbggTools/test/RunJobs/EGMMVA_CorrPreSel/Hadd/output_GluGluToHHTo2B2G_node_THENODE_13TeV-madgraph.root"
Data = "/afs/cern.ch/work/r/rateixei/work/DiHiggs/flg76X/CMSSW_7_6_3/src/flashgg/bbggTools/test/RunJobs/EGMMVA_CorrPreSel/Hadd/DoubleEG.root"

#dirPrefix = "EGM_MVA_2cat_MixedMod4_low0.46_high0.936"
dirPrefix = "EGM_MVA_2cat_MediumWP"

postFix = " -MX -btagWP 0.8"
#postFix = " -MX -doCatMixed -btagHigh 0.936 -btagLow 0.46"
#postFix = " -MX -tilt -singleCat -btagWP 0.46"

directory = dirPrefix
os.system( "mkdir " + directory+"_LowMass" )
os.system( "mkdir " + directory+"_HighMass" )

print "DOING LowMassCat Data"
command = "LimitTreeMaker -inputFile " + Data + " -o " +   directory+"_LowMass" + " -min 0 -max 350 -scale 1." + postFix
os.system(command)
print "DOING HighMassCat Data"
command = "LimitTreeMaker -inputFile " + Data + " -o " +   directory+"_HighMass" + " -min 350 -max 35000 -scale 1." + postFix
os.system(command)

for MM in nodes:
	i = MM[0]
	sigScale = 2.7/float(MM[1])
	print "DOING LowMassCat Signal, node ", i
	command = "LimitTreeMaker -inputFile " + Signals.replace("THENODE", str(i)) + " -o " + directory+"_LowMass" + " -min 0 -max 350 -scale " + str(sigScale) + postFix
	print command
	os.system(command)
	print "DOING HighMassCat Signal, node ", i
	command = "LimitTreeMaker -inputFile " + Signals.replace("THENODE", str(i)) + " -o " + directory+"_HighMass" + " -min 350 -max 35000 -scale " + str(sigScale) + postFix
	print command
	os.system(command)
	
