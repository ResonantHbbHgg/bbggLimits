from ROOT import *
import os


masses = [ [250, 49600], [260, 50000], [280, 50000], [300, 47600], [320, 50000], [350, 50000], [400, 50000], [450, 49800], [500, 50000], [600, 50000], [700, 50000], [800, 50000], [900, 50000] ]

dirPrefix = "EGM_MVA_2cat_MixedMod4_low0.46_high0.936_PCR"

#postFix = " -MX -tilt -btagWP 0.46 -photonCR"
postFix = " -MX -tilt -doCatMixed -btagHigh 0.936 -btagLow 0.46 -photonCR"
#postFix = " -MX -tilt -singleCat -btagWP 0.46"

#pol1
p0 = -34.9426
p1 = 0.17724

def minMass(m):
        width = p0 + p1*m# + p2*m*m + p3*m*m*m
        mmass = m - (width)/2.
        return mmass
#       return (m-10)

def maxMass(m):
        width = p0 + p1*m# + p2*m*m + p3*m*m*m
        mmass = m + (width)/2.
        return mmass
#       return (m+10)


for MM in masses:
	i = MM[0]
	sigScale = 2.7/float(MM[1])
	directory = dirPrefix + "_" + str(i)
	print "Processing mass:", i
	print "----  mass window:", minMass(i)," - ",maxMass(i)
	os.system( "mkdir " + directory )
	print "DOING SIGNAL"
	command = "LimitTreeMaker -i Radion.txt -o " + directory + " -min " + str(minMass( i)) + " -max " + str(maxMass( i)) + " -scale " + str(sigScale) + postFix
	print command
#	os.system(command)
	print "DOING DATA"
	command = "LimitTreeMaker -i data.txt -o " + directory + " -min " + str(minMass( i)) + " -max " + str(maxMass( i)) + " -scale 1." + postFix
	print command
	os.system(command)
	continue
	print "DOING DIPHOTON"
	command = "LimitTreeMaker -i dipho.txt -o " + directory + " -min " + str(minMass( i)) + " -max " + str(maxMass( i)) + " -scale 0.0106508" + postFix
	print command
	os.system(command)
	print "DOING GJ20"
	command = "LimitTreeMaker -i g20.txt -o " + directory + " -min " + str(minMass( i)) + " -max " + str(maxMass( i)) + " -scale 1.83018" + postFix
	print command
	os.system(command)
	print "DOING GJ40"
	command = "LimitTreeMaker -i g40.txt -o " + directory + " -min " + str(minMass( i)) + " -max " + str(maxMass( i)) + " -scale 0.255013" + postFix
	print command
	os.system(command)
