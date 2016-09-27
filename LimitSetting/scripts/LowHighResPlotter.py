from ROOT import *
gROOT.SetBatch(True)
import sys, getopt, os
from array import array
from ResonantCrossSections import *

col1 = TColor(11111, 249./255., 243./255., 33./255.)

def main(argv):
	lowfolder = ""
	highfolder = ""
	lumi = ""
	doObs = 0
	radion = 0
	graviton = 0
	try:
		opts, args = getopt.getopt(argv,"L:H:l:oRG",["LowFolder=", "HighFolder=", "lumi=", "observed", "radion", "graviton"])
	except getopt.GetoptError:
		print 'ResPlotter.py -f <folder> -l lumi'
		sys.exit(2)
	for opt, arg in opts:
		if opt == "-f":
			folder = arg
		if opt == "-l":
			lumi = arg
		if opt == "-L":
			lowfolder = arg
		if opt == "-H":
			highfolder = arg
		if opt == "-o":
			doObs = 1
		if opt == "-R":
			radion = 1
		if opt == "-G":
			graviton = 1

	if graviton == 1 and radion == 1:
		print "Either choose -R or -G..."
		sys.exit(2)

	if lowfolder == "" or highfolder == "":
		print 'ResPlotter.py -f <folder>'
		sys.exit(2)

	if not os.path.isdir(lowfolder) or not os.path.isdir(highfolder):
		print "FOLDER DOESN'T EXIST!"
		sys.exit(2)

	dirs = os.listdir(lowfolder)
	massList = []
	centralVal = []
	s1_up = []
	s1_down = []
	s2_up = []
	s2_down = []
	zeros = []
	observed = []
	Maximum = 0
	for d in dirs:
		if "Radion" not in d and "Graviton" not in d:
			continue
		if(radion): mass = int(d.replace("Radion_M", ""))
		if(graviton): mass = int(d.replace("Graviton_M", ""))
		if mass > 500: continue
		massList.append( mass )
		limFile = open(lowfolder+"/"+d+"/logs/higgsCombineTest.Asymptotic.mH125.0._m"+str(mass)+"_higgs.txt")
		cenVal = -1
		for line in limFile:
			if "50.0%" in line:
				cenVal = float(line.split(" < ")[1])
				centralVal.append( cenVal )
		if cenVal == -1:
			print "CENTRAL VALUE NOT FOUND!"
			sys.exit(2)
		limFile = open(lowfolder+"/"+d+"/logs/higgsCombineTest.Asymptotic.mH125.0._m"+str(mass)+"_higgs.txt")
		zeros.append(0)
		for line in limFile:
			if "Observed" in line:
				val = float( line.split(" < ")[1] )
				observed.append(val)
				continue
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

	gr_observed = TGraph(len(massList), array('d', massList), array('d', observed))

	highdirs = sorted(os.listdir(highfolder))
	h_massList = []
	h_centralVal = []
	h_s1_up = []
	h_s1_down = []
	h_s2_up = []
	h_s2_down = []
	h_zeros = []
	h_observed = []
	for d in highdirs:
		if "Radion" not in d and "Graviton" not in d:
			continue
		if(radion): mass = int(d.replace("Radion_M", ""))
		if(graviton): mass = int(d.replace("Graviton_M", ""))
		if int(mass) < 500: continue
		h_massList.append( mass )
		limFile = open(highfolder+"/"+d+"/logs/higgsCombineTest.Asymptotic.mH125.0._m"+str(mass)+"_higgs.txt")
		cenVal = -1
		for line in limFile:
			if "50.0%" in line:
				cenVal = float(line.split(" < ")[1])
				h_centralVal.append( cenVal )
		if cenVal == -1:
			print "CENTRAL VALUE NOT FOUND!"
			sys.exit(2)
		limFile = open(highfolder+"/"+d+"/logs/higgsCombineTest.Asymptotic.mH125.0._m"+str(mass)+"_higgs.txt")
		h_zeros.append(0)
		for line in limFile:
			if "Observed" in line:
				val = float(line.split(" < ")[1])
				h_observed.append( val )
			if "2.5%" in line:
				val = cenVal - float(line.split(" < ")[1])
				h_s2_down.append( val )
				continue
			if "16.0%" in line:
				val = cenVal - float(line.split(" < ")[1])
				h_s1_down.append( val )
				continue
			if "84.0%" in line:
				val = float(line.split(" < ")[1]) - cenVal
				h_s1_up.append( val )
				continue
			if "97.5%" in line:
				val = float(line.split(" < ")[1]) - cenVal
				h_s2_up.append( val )
				continue

	h_gr_centralVal_1s = TGraphAsymmErrors(len(h_massList), array('d', h_massList), array('d', h_centralVal), array('d', h_zeros), array('d', h_zeros), array('d', h_s1_down), array('d', h_s1_up))
	h_gr_centralVal_2s = TGraphAsymmErrors(len(h_massList), array('d', h_massList), array('d', h_centralVal), array('d', h_zeros), array('d', h_zeros), array('d', h_s2_down), array('d', h_s2_up))

	h_gr_observed = TGraph(len(h_massList), array('d', h_massList), array('d', h_observed))

	h_gr_observed.SetLineWidth(3)
	h_gr_observed.SetLineColor(kBlack)
	h_gr_observed.SetLineStyle(kDashed)
	h_gr_observed.SetMarkerStyle(20)
	gr_observed.SetLineWidth(3)
	gr_observed.SetLineColor(kBlack)
	gr_observed.SetLineStyle(kDashed)
	gr_observed.SetMarkerStyle(20)
	
	leg = ""
	if(graviton or radion): leg = TLegend(0.5, 0.7, 0.89, 0.87)
	if(0): leg = TLegend(0.5, 0.51, 0.89, 0.7)
	leg.SetFillStyle(0)
	leg.SetLineWidth(0)
	leg.SetBorderSize(0)
	
	c0 = TCanvas("c", "c", 800, 600)
	c0.SetLogy()
#	gr_centralVal_2s.SetMaximum(Maximum*2.5)
	gr_centralVal_2s.SetMaximum(400)
	gr_centralVal_2s.SetMinimum(0.8)
	gr_centralVal_2s.Draw("AP3")
	gr_centralVal_2s.SetLineColor(0)
#	gr_centralVal_2s.SetFillColor(kGreen+1)
	gr_centralVal_2s.SetFillColor(11111)
	gr_centralVal_2s.SetLineWidth(4)
	gr_centralVal_2s.SetTitle("")
#	gr_centralVal_2s.GetXaxis().SetLimits(massList[0] - 20, massList[len(massList)-1]+20)
	gr_centralVal_2s.GetXaxis().SetLimits(massList[0] - 20, h_massList[len(h_massList)-1]+20)
	gr_centralVal_2s.GetXaxis().SetTitle("M_{X} [GeV]")
	gr_centralVal_2s.GetXaxis().SetTitleSize(0.045)
	gr_centralVal_2s.GetYaxis().SetTitleSize(0.045)
	gr_centralVal_2s.GetYaxis().SetTitle("#sigma(pp#rightarrowX#rightarrowHH#rightarrowbb#gamma#gamma) [fb]")
	c0.Update()
	gr_centralVal_1s.Draw("3same")
	gr_centralVal_1s.SetMarkerStyle(21)
	gr_centralVal_1s.SetMarkerColor(kBlue+1)
	gr_centralVal_1s.SetLineWidth(3)
	gr_centralVal_1s.SetLineColor(kBlue+1)
#	gr_centralVal_1s.SetFillColor(kYellow)
#	gr_centralVal_1s.SetFillColor(11111)
	gr_centralVal_1s.SetFillColor(kGreen+1)
	c0.Update()
	gr_centralVal_1s.Draw("XLsame")
	sigma1 = gr_centralVal_1s.Clone()
	sigma1.SetLineColor(0)

#	h_gr_centralVal_2s.SetFillColor(kGreen+1)
	h_gr_centralVal_2s.SetFillColor(11111)
	h_gr_centralVal_2s.SetLineColor(0)
	h_gr_centralVal_2s.Draw("3same")
#	h_gr_centralVal_1s.SetFillColor(kYellow)
#	h_gr_centralVal_1s.SetFillColor(11111)
	h_gr_centralVal_1s.SetFillColor(kGreen+1)
	h_gr_centralVal_1s.Draw("3same")
	c0.Update()
	h_gr_centralVal_1s.SetLineColor(kBlue+1)
	h_gr_centralVal_1s.SetLineWidth(3)
	h_gr_centralVal_1s.Draw("XLsame")
	h_sigma1 = h_gr_centralVal_1s.Clone()
	h_sigma1.SetLineColor(0)

	if(doObs): h_gr_observed.Draw("PLsame")
	if(doObs): gr_observed.Draw("PLsame")
	if(doObs): leg.AddEntry(gr_observed, "Observed 95% upper limit", "lp")

	leg.SetTextSize(0.035)

	if(graviton or radion):	legRes = TLegend(0.13, 0.7, 0.5, 0.87)
	if(0):	legRes = TLegend(0.5, 0.7, 0.89, 0.89)
	legRes.SetFillStyle(0)
	legRes.SetLineWidth(0)
	legRes.SetBorderSize(0)
	if(radion):
		legRes.SetHeader("#splitline{pp#rightarrowX#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}{Spin-0 Resonance}")
		legRes.AddEntry(gr_rad, "Bulk Radion, #Lambda_{R} = 1 TeV", "l")
		gr_rad.Draw("Lsame")
	if(graviton):
		legRes.SetHeader("#splitline{pp#rightarrowX#rightarrowHH#rightarrowb#bar{b}#gamma#gamma}{Spin-2 Resonance}")
		legRes.AddEntry(gr_rad, "Bulk Graviton, k/#bar{M}_{Pl} = 1", "l")
		gr_grav.Draw("Lsame")
	legRes.SetTextSize(0.038)
	legRes.Draw("same")

	leg.AddEntry(gr_centralVal_1s, "Expected 95% upper limit", "l")
	leg.AddEntry(sigma1, "Expected limit #pm 1#sigma", "f")
	leg.AddEntry(gr_centralVal_2s, "Expected limit #pm 2#sigma", "f")
#	leg.AddEntry(h_gr_centralVal_1s, "High Mass Expected 95% upper limit", "l")
#	leg.AddEntry(h_sigma1, "High Mass Expected limit #pm 1#sigma", "f")
#	leg.AddEntry(h_gr_centralVal_2s, "High Mass Expected limit #pm 2#sigma", "f")
	leg.Draw("same")

	line = ""
	if(0): line = TLine(500, 0, 500, gr_centralVal_2s.GetYaxis().GetXmax())
	if(graviton or radion): line = TLine(500, 0, 500, gr_centralVal_2s.GetYaxis().GetXmax()*0.15)
	line.SetLineStyle(2)
	line.Draw("same")
	
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
	tlatex.SetTextSize(23)
#	tlatex.DrawLatex(0.6, 0.91, "#sqrt{s} = 13 TeV, L = " + str(lumi) + " fb^{-1}")
	tlatex.DrawLatex(0.65, 0.91,"L = " + str(lumi) + " fb^{-1} (13 TeV)")
	c0.SaveAs(lowfolder+"/Low_High_ResPlot.pdf")

	sys.exit(2)	

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
	tlatex.SetTextSize(23)
	tlatex.DrawLatex(0.11, 0.91, "CMS")
	tlatex.SetTextFont(53)
	tlatex.DrawLatex(0.18, 0.91, "Preliminary")
	tlatex.SetTextFont(43)
	tlatex.DrawLatex(0.6, 0.91, "#sqrt{s} = 13 TeV, L = " + str(lumi) + " fb^{-1}")
	leg.Draw("same")
	c0.SaveAs(folder+"/ResPlot_hh.pdf")


if __name__ == "__main__":
	main(sys.argv[1:])
