from ROOT import *
import os


massesRad = [ [250, 49800], [260, 50000], [270, 48400], [280, 50000], [300, 50000], [320, 50000], [340, 50000], [350, 50000], [400, 50000], [450, 50000], [500, 50000], [550, 49200],[600, 50000], [650, 49200], [700, 50000], [750, 50000], [800, 48400], [900, 50000] ]

massesGrav = [ [250, 50000], [260, 49200], [270, 50000], [280, 484800], [300, 50000], [320, 49200], [340, 49200], [350, 49200], [400, 50000], [450, 50000], [500, 50000], [550, 50000],  [600, 50000], [650, 50000], [700, 50000], [750, 50000], [800, 48200], [900, 50000] ]


masses = massesRad

SignalFiles = "radion.txt"
DataFiles = "data.txt"

#Foldername
dirPrefix = "HighMassAnalysis"

#Postfix to produce low mass samples
#postFix = " -MX -tilt -doCatMixed -btagHigh 0.936 -btagLow 0.46"

#Postfix to produce high mass samples
postFix = " -MX -tilt -singleCat -btagWP 0.46"

#Apply central scalefactors (only MC)
SFs = " -doBVariation 0 -doPhoVariation 0"

#pol1
p0 = -34.9426
p1 = 0.17724


#Mass window definitions
def minMass(m):
        width = p0 + p1*m# + p2*m*m + p3*m*m*m
        mmass = m - (width)/2.
        return mmass

def maxMass(m):
        width = p0 + p1*m# + p2*m*m + p3*m*m*m
        mmass = m + (width)/2.
        return mmass


for MM in masses:
	i = MM[0]
	sigScale = 2.7/float(MM[1])
	directory = dirPrefix + "_" + str(i)
	print "Processing mass:", i
	print "----  mass window:", minMass(i)," - ",maxMass(i)
	os.system( "mkdir " + directory )
	print "DOING SIGNAL"
	command = "LimitTreeMaker -i "+SignalFiles+" -o " + directory + " -min " + str(minMass( i)) + " -max " + str(maxMass( i)) + " -scale " + str(sigScale) + postFix + SFs
	print command
	os.system(command)
#	continue
	print "DOING DATA"
	command = "LimitTreeMaker -i "+DataFiles+" -o " + directory + " -min " + str(minMass( i)) + " -max " + str(maxMass( i)) + " -scale 1." + postFix
	print command
	os.system(command)
	continue
	print "DOING DIPHOTON"
	command = "LimitTreeMaker -i dipho.txt -o " + directory + " -min " + str(minMass( i)) + " -max " + str(maxMass( i)) + " -scale 0.0106508" + postFix + SFs
	print command
	os.system(command)
	print "DOING GJ20"
	command = "LimitTreeMaker -i g20.txt -o " + directory + " -min " + str(minMass( i)) + " -max " + str(maxMass( i)) + " -scale 1.83018" + postFix + SFs
	print command
	os.system(command)
	print "DOING GJ40"
	command = "LimitTreeMaker -i g40.txt -o " + directory + " -min " + str(minMass( i)) + " -max " + str(maxMass( i)) + " -scale 0.255013" + postFix + SFs
	print command
	os.system(command)
