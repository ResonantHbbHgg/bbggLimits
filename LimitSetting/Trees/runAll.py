from ROOT import *
import os

masses = [ 250, 260, 300, 320, 350, 400, 450, 500, 600, 700, 800, 900 ]
dirPrefix = "trees_mass_20GeV"

p0 = -118.815
p1 = 0.722824
p2 = -0.000994524
p3 = 5.32791e-07

def minMass(m):
#	width = p0 + p1*m + p2*m*m + p3*m*m*m
#	mmass = m - (width)/2.
#	return mmass
	return (m-10)

def maxMass(m):
#	width = p0 + p1*m + p2*m*m + p3*m*m*m
#	mmass = m + (width)/2.
#	return mmass
	return (m+10)

for i in masses:
	directory = dirPrefix + "_" + str(i)
	print "Processing mass:", i
	print "----  mass window:", minMass(i)," - ",maxMass(i)
	os.system( "mkdir " + directory )
	print "DOING SIGNAL"
#	command = "LimitTreeMaker -i Radion.txt -o " + directory + " -min " + str(i-5) + " -max " + str(i+5) + "-scale 5.4e-05 -MX"
	command = "LimitTreeMaker -i Radion.txt -o " + directory + " -min " + str(minMass( i)) + " -max " + str(maxMass( i)) + " -scale 5.4e-05 -MX"
	print command
	os.system(command)
	print "DOING DATA"
	command = "LimitTreeMaker -i data.txt -o " + directory + " -min " + str(minMass( i)) + " -max " + str(maxMass( i)) + " -scale 1. -MX"
#	command = "LimitTreeMaker -i data.txt -o " + directory + " -min " + str(i-5) + " -max " + str(i+5) + "-scale 1. -MX"
	print command
	os.system(command)
#	LimitTreeMaker -i data.txt -o ${dirPrefix}_$i -min $((i-5)) -max $((i+5)) -scale 1. -MX
	print "DOING DIPHOTON"
#	command = "LimitTreeMaker -i dipho.txt -o " + directory + " -min " + str(i-5) + " -max " + str(i+5) + "-scale 0.0106508 -MX"
	command = "LimitTreeMaker -i dipho.txt -o " + directory + " -min " + str(minMass( i)) + " -max " + str(maxMass( i)) + " -scale 0.0106508 -MX"
	print command
	os.system(command)
#	LimitTreeMaker -i dipho.txt -o ${dirPrefix}_$i -min $((i-5)) -max $((i+5)) -scale 0.0106508 -MX
	print "DOING GJ20"
#	command = "LimitTreeMaker -i g20.txt -o " + directory + " -min " + str(i-5) + " -max " + str(i+5) + "-scale 1.83018 -MX"
	command = "LimitTreeMaker -i g20.txt -o " + directory + " -min " + str(minMass( i)) + " -max " + str(maxMass( i)) + " -scale 1.83018 -MX"
	print command
	os.system(command)
#	LimitTreeMaker -i gj20.txt -o ${dirPrefix}_$i -min $((i-5)) -max $((i+5)) -scale 1.83018 -MX
	print "DOING GJ40"
#	LimitTreeMaker -i gj40.txt -o ${dirPrefix}_$i -min $((i-5)) -max $((i+5)) -scale 0.255013 -MX
#	command = "LimitTreeMaker -i g40.txt -o " + directory + " -min " + str(i-5) + " -max " + str(i+5) + "-scale 0.255013 -MX"
	command = "LimitTreeMaker -i g40.txt -o " + directory + " -min " + str(minMass( i)) + " -max " + str(maxMass( i)) + " -scale 0.255013 -MX"
	print command
	os.system(command)
