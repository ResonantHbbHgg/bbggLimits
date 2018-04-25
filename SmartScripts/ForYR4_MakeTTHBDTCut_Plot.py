#!/usr/bin/env python

from ROOT import *

gOpt = TGraph()
iCounter = 0

gROOT.SetBatch()

for x in range(-10,6):
    
    iCounter = iCounter+1
    cut=-0.1+x*0.02
    print cut
    if cut < 0:
        file = "outDir_m"+str(abs(cut*10))+"/CombinedCard_ARW_/result_1.log"
    else:
        file = "outDir_"+str(abs(cut*10))+"/CombinedCard_ARW_/result_1.log"

    print file

    with open(file) as f:
        for line in f:
            if line.find("50.0%") > -1:
                print line
                lhs, rhs = line.split("<",1)
                print lhs
                limit = float(rhs)
                print "limit = "+str(limit)
                gOpt.SetPoint(iCounter, cut, limit)
                
    
c = TCanvas("c", "c", 800, 600)
plotter = TH1F("","", 1, -0.3, 0.3)
plotter.Draw()
plotter.SetMaximum(1.8)
plotter.SetMinimum(1.5)
#c.SetGrid()
gOpt.SetMarkerStyle(25)

gOpt.Draw("P")

c.SaveAs("Optimisation.png")
