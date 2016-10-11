from ROOT import *
import os


nodes = [ ["box", 49600], ["SM", 50000], [2, 50000], [3, 47600], [4, 50000], [5, 50000], [6, 50000], [7, 49800], [8, 50000], [9, 50000], [10, 50000], [11, 50000], [12, 50000], [13, 50000] ]

Signals = "/tmp/rateixei/eos/cms/store/user/rateixei/HHbbgg/FlatTrees/ICHEP_Regressed4b/output_GluGluToHHTo2B2G_node_THENODE_13TeV-madgraph.root"
Data = "/tmp/rateixei/eos/cms/store/user/rateixei/HHbbgg/FlatTrees/ICHEP_Regressed4b/DoubleEG.root"

dirPrefix = "NonResAnalysis_ICHEP"

postFix = " -MX -btagWP 0.8 "

SFs = " -doBVariation 0 -doPhoVariation 0"

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
	command = "LimitTreeMaker -inputFile " + Signals.replace("THENODE", str(i)) + " -o " + directory+"_LowMass" + " -min 0 -max 350 -scale " + str(sigScale) + postFix + SFs
	print command
	os.system(command)
#	continue
	print "DOING HighMassCat Signal, node ", i
	command = "LimitTreeMaker -inputFile " + Signals.replace("THENODE", str(i)) + " -o " + directory+"_HighMass" + " -min 350 -max 35000 -scale " + str(sigScale) + postFix + SFs
	print command
	os.system(command)
	
