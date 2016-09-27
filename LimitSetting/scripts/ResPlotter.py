from ROOT import *
gROOT.SetBatch(True)
import sys, getopt, os
from array import array

def main(argv):
	folder = ""
	lumi = ""
	obs = 0
	try:
		opts, args = getopt.getopt(argv,"f:l:o",["folder=", "lumi=", "observed"])
	except getopt.GetoptError:
		print 'ResPlotter.py -f <folder> -l lumi'
		sys.exit(2)
	for opt, arg in opts:
		if opt == "-f":
			folder = arg
		if opt == "-l":
			lumi = arg
		if opt == "-o":
			obs = 1
	if folder == "":
		print 'ResPlotter.py -f <folder>'
		sys.exit(2)

	if not os.path.isdir(folder):
		print "FOLDER DOESN'T EXIST!"
		sys.exit(2)

	dirs = os.listdir(folder)
	massList = []
	centralVal = []
	s1_up = []
	s1_down = []
	s2_up = []
	s2_down = []
	zeros = []
	Maximum = 0
	observed = []
	for d in dirs:
		isRad = 0
		isGrav = 0
		if "Radion" in d:
			isRad = 1
		if "Graviton" in d:
			isGrav = 1
		if isRad == 0 and isGrav == 0:
			continue
		ffold = "Radion_M"
		if isGrav == 1:
			ffold = "Graviton_M"
		mass = int(d.replace(ffold, ""))
		massList.append( mass )
		limFile = open(folder+"/"+d+"/logs/higgsCombineTest.Asymptotic.mH125.0._m"+str(mass)+"_higgs.txt")
		cenVal = -1
		for line in limFile:
			if "50.0%" in line:
				cenVal = float(line.split(" < ")[1])
				centralVal.append( cenVal )
		if cenVal == -1:
			print "CENTRAL VALUE NOT FOUND!"
			sys.exit(2)
		limFile = open(folder+"/"+d+"/logs/higgsCombineTest.Asymptotic.mH125.0._m"+str(mass)+"_higgs.txt")
		zeros.append(0)
		for line in limFile:
			if "Observed Limit: r <" in line:
				val = float(line.split(" < ")[1])
				observed.append(val)
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

	gr_centralVal_1s = TGraphAsymmErrors(len(massList), array('d', massList), array('d', centralVal), array('d', zeros), array('d', zeros), array('d', s1_down), array('d', s1_up))
	gr_centralVal_2s = TGraphAsymmErrors(len(massList), array('d', massList), array('d', centralVal), array('d', zeros), array('d', zeros), array('d', s2_down), array('d', s2_up))
	gr_observed = TGraph(len(massList), array('d', massList), array('d',observed))
	gr_observed.SetLineWidth(3)
	gr_observed.SetLineColor(kBlack)

	print "Mass", "Expected", "Observed"
	for jj,m in enumerate(massList):
		if obs: print m, centralVal[jj], observed[jj]
		if not obs: print m, centralVal[jj]

	leg = TLegend(0.65, 0.7, 0.89, 0.89)
	leg.SetFillStyle(0);
	leg.SetLineWidth(0);
	leg.SetBorderSize(0);
	
	c0 = TCanvas("c", "c", 800, 600)
	gr_centralVal_2s.SetMaximum(Maximum*2)
	gr_centralVal_2s.Draw("AP3")
	gr_centralVal_2s.SetLineColor(0)
	gr_centralVal_2s.SetFillColor(kGreen+1)
	gr_centralVal_2s.SetLineWidth(4)
	gr_centralVal_2s.SetTitle("")
	gr_centralVal_2s.GetXaxis().SetLimits(massList[0] - 20, massList[len(massList)-1]+20)
	gr_centralVal_2s.GetXaxis().SetTitle("M(X) [GeV]")
	gr_centralVal_2s.GetYaxis().SetTitle("#sigma(pp#rightarrowX#rightarrowHH#rightarrowbb#gamma#gamma) [fb]")
	c0.Update()
	gr_centralVal_1s.Draw("3same")
	gr_centralVal_1s.SetMarkerStyle(21)
	gr_centralVal_1s.SetMarkerColor(kBlue+1)
	gr_centralVal_1s.SetLineWidth(3)
	gr_centralVal_1s.SetLineColor(kBlue+1)
	gr_centralVal_1s.SetFillColor(kYellow)
	c0.Update()
	gr_centralVal_1s.Draw("XLsame")
	sigma1 = gr_centralVal_1s.Clone()
	sigma1.SetLineColor(0)
	if(obs): gr_observed.Draw("Lsame")
	if(obs): leg.AddEntry(gr_observed, "Observed 95% upper limit", "l")
	leg.AddEntry(gr_centralVal_1s, "Expected 95% upper limit", "l")
	leg.AddEntry(sigma1, "Expected limit #pm 1#sigma", "f")
	leg.AddEntry(gr_centralVal_2s, "Expected limit #pm 2#sigma", "f")
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
	c0.SaveAs(folder+"/ResPlot.pdf")
	
	gr_centralVal_1s_hh = TGraphAsymmErrors(len(massList), array('d', massList), array('d', [i/2.6 for i in centralVal]), array('d', zeros), array('d', zeros), array('d', [i/2.6 for i in s1_down]), array('d', [i/2.6 for i in s1_up]))
	gr_centralVal_2s_hh = TGraphAsymmErrors(len(massList), array('d', massList), array('d', [i/2.6 for i in centralVal]), array('d', zeros), array('d', zeros), array('d', [i/2.6 for i in s2_down]), array('d', [i/2.6 for i in s2_up]))

	c0 = TCanvas("c", "c", 800, 600)
	gr_centralVal_2s_hh.SetMaximum(Maximum*2/2.6)
	gr_centralVal_2s_hh.Draw("AP")
	gr_centralVal_2s_hh.SetLineColor(kGreen+1)
	gr_centralVal_2s_hh.SetLineWidth(4)
	gr_centralVal_2s_hh.SetTitle("")
	gr_centralVal_2s_hh.GetXaxis().SetLimits(massList[0] - 20, massList[len(massList)-1]+20)
	gr_centralVal_2s_hh.GetXaxis().SetTitle("Benchmark Points")
	gr_centralVal_2s_hh.GetYaxis().SetTitle("#sigma(pp#rightarrowX#rightarrowHH)/Br(HH#rightarrowbb#gamma#gamma) [pb]")
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
	c0.SaveAs(folder+"/ResPlot_hh.pdf")


if __name__ == "__main__":
	main(sys.argv[1:])
