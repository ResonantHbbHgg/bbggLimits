from ROOT import *
gROOT.SetBatch(True)
import sys, getopt, os
from array import array

def main(argv):
	folder = ""
	lumi = ""
	try:
		opts, args = getopt.getopt(argv,"f:l:",["folder=", "lumi="])
	except getopt.GetoptError:
		print 'NonResonantOrganizer.py -f <folder> -l lumi'
		sys.exit(2)
	for opt, arg in opts:
		if opt == "-f":
			folder = arg
		if opt == "-l":
			lumi = arg
	if folder == "":
		print 'NonResonantOrganizer.py -f <folder>'
		sys.exit(2)

	if not os.path.isdir(folder):
		print "FOLDER DOESN'T EXIST!"
		sys.exit(2)

	dirs = os.listdir(folder)
	nodesList = []
	centralVal = []
	s1_up = []
	s1_down = []
	s2_up = []
	s2_down = []
	zeros = []
	Maximum = 0
	for d in dirs:
		if "Node" not in d:
			continue
		nodeNumber = int(d.replace("Node", ""))
		nodesList.append( nodeNumber )
		limFile = open(folder+"/"+d+"/logs/higgsCombineTest.Asymptotic.mH125.0._m"+str(nodeNumber)+"_higgs.txt")
		cenVal = -1
		for line in limFile:
			if "50.0%" in line:
				cenVal = float(line.split(" < ")[1])
				centralVal.append( cenVal )
		if cenVal == -1:
			print "CENTRAL VALUE NOT FOUND!"
			sys.exit(2)
		limFile = open(folder+"/"+d+"/logs/higgsCombineTest.Asymptotic.mH125.0._m"+str(nodeNumber)+"_higgs.txt")
		zeros.append(0)
		for line in limFile:
			if "2.5%" in line:
				val = cenVal - float(line.split(" < ")[1])
				s2_down.append( val )
				continue
			if "16.0%" in line:
				val = cenVal - float(line.split(" < ")[1])
				s1_down.append( val )
				continue
			if "84.0%" in line:
				val = float(line.split(" < ")[1]) - cenVal
				s1_up.append( val )
				continue
			if "97.5%" in line:
				val = float(line.split(" < ")[1]) - cenVal
				s2_up.append( val )
				if val > Maximum:
					Maximum = val
				continue
	print s1_down, s2_down
	gr_centralVal_1s = TGraphAsymmErrors(len(nodesList), array('d', nodesList), array('d', centralVal), array('d', zeros), array('d', zeros), array('d', s1_down), array('d', s1_up))
	gr_centralVal_2s = TGraphAsymmErrors(len(nodesList), array('d', nodesList), array('d', centralVal), array('d', zeros), array('d', zeros), array('d', s2_down), array('d', s2_up))
	
	leg = TLegend(0.65, 0.7, 0.89, 0.89)
	leg.SetFillStyle(0);
	leg.SetLineWidth(0);
	leg.SetBorderSize(0);
	
	c0 = TCanvas("c", "c", 800, 600)
	gr_centralVal_2s.SetMaximum(Maximum*2)
	gr_centralVal_2s.Draw("AP")
	gr_centralVal_2s.SetLineColor(kGreen+1)
	gr_centralVal_2s.SetLineWidth(4)
	gr_centralVal_2s.SetTitle("")
	gr_centralVal_2s.GetXaxis().SetLimits(-1, 14)
	gr_centralVal_2s.GetXaxis().SetTitle("Benchmark Points")
	gr_centralVal_2s.GetYaxis().SetTitle("#sigma(pp#rightarrowHH#rightarrowbb#gamma#gamma) [fb]")
	c0.Update()
	gr_centralVal_1s.Draw("EP")
	gr_centralVal_1s.SetMarkerStyle(21)
	gr_centralVal_1s.SetMarkerColor(kBlue+1)
	gr_centralVal_1s.SetLineWidth(5)
	gr_centralVal_1s.SetLineColor(kYellow)
	
	leg.AddEntry(gr_centralVal_1s, "Expected 95% upper limit", "p")
	leg.AddEntry(gr_centralVal_1s, "Expected limit #pm 1#sigma", "l")
	leg.AddEntry(gr_centralVal_2s, "Expected limit #pm 2#sigma", "l")
	leg.Draw("same")
	
	tlatex = TLatex()
	tlatex.SetNDC()
	tlatex.SetTextAngle(0)
	tlatex.SetTextColor(kBlack)
	tlatex.SetTextFont(63)
	tlatex.SetTextAlign(11)
	tlatex.SetTextSize(25)
	tlatex.DrawLatex(0.11, 0.91, "CMS")
	tlatex.SetTextFont(53)
	tlatex.DrawLatex(0.18, 0.91, "Preliminary")
	tlatex.SetTextFont(43)
	tlatex.DrawLatex(0.6, 0.91, "#sqrt{s} = 13 TeV, L = " + str(lumi) + " fb^{-1}")
	c0.SaveAs(folder+"/NonResPlot.pdf")
	
	gr_centralVal_1s_hh = TGraphAsymmErrors(len(nodesList), array('d', nodesList), array('d', [i/2.6 for i in centralVal]), array('d', zeros), array('d', zeros), array('d', [i/2.6 for i in s1_down]), array('d', [i/2.6 for i in s1_up]))
	gr_centralVal_2s_hh = TGraphAsymmErrors(len(nodesList), array('d', nodesList), array('d', [i/2.6 for i in centralVal]), array('d', zeros), array('d', zeros), array('d', [i/2.6 for i in s2_down]), array('d', [i/2.6 for i in s2_up]))

	c0 = TCanvas("c", "c", 800, 600)
	gr_centralVal_2s_hh.SetMaximum(Maximum*2/2.6)
	gr_centralVal_2s_hh.Draw("AP")
	gr_centralVal_2s_hh.SetLineColor(kGreen+1)
	gr_centralVal_2s_hh.SetLineWidth(4)
	gr_centralVal_2s_hh.SetTitle("")
	gr_centralVal_2s_hh.GetXaxis().SetLimits(-1, 14)
	gr_centralVal_2s_hh.GetXaxis().SetTitle("Benchmark Points")
	gr_centralVal_2s_hh.GetYaxis().SetTitle("#sigma(pp#rightarrowHH)/Br(HH#rightarrowbb#gamma#gamma) [pb]")
	c0.Update()
	gr_centralVal_1s_hh.Draw("EP")
	gr_centralVal_1s_hh.SetMarkerStyle(21)
	gr_centralVal_1s_hh.SetMarkerColor(kBlue+1)
	gr_centralVal_1s_hh.SetLineWidth(5)
	gr_centralVal_1s_hh.SetLineColor(kYellow)

	tlatex.SetTextAngle(0)
	tlatex.SetTextColor(kBlack)
	tlatex.SetTextFont(63)
	tlatex.SetTextAlign(11)
	tlatex.SetTextSize(25)
	tlatex.DrawLatex(0.11, 0.91, "CMS")
	tlatex.SetTextFont(53)
	tlatex.DrawLatex(0.18, 0.91, "Preliminary")
	tlatex.SetTextFont(43)
	tlatex.DrawLatex(0.6, 0.91, "#sqrt{s} = 13 TeV, L = " + str(lumi) + " fb^{-1}")
	leg.Draw("same")
	c0.SaveAs(folder+"/NonResPlot_hh.pdf")


if __name__ == "__main__":
	main(sys.argv[1:])
