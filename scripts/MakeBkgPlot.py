#!/usr/bin/env python

from ROOT import *
from HiggsAnalysis.bbggLimits.BkgPlotter import *
import sys, getopt, os

def main(argv):
	gSystem.Load("libHiggsAnalysisCombinedLimit.so")

	wfile = ""
	cat = -1
	obs = ""
	lumi = ""
	doBands = 1
	analysis = ""
	Blinded = 0
	myLabel = ""
	bins = []
	try:
		opts, args = getopt.getopt(argv,"w:c:o:l:a:b:BL:",["workspace=", "cat=", "observable=","lumi=","analysis=","bins=", "Blind", "Label"])
	except getopt.GetoptError:
		print 'MakeBkgPlot.py -w <workspace file> -c <cat> -o <obs1,obs2,obs3> -l <lumi> -a <analysis title> -b <binobs1,binobs2> -B (blinded)'
		sys.exit(2)
	for opt, arg in opts:
		if opt == "-w":
			wfile = arg
		if opt == "-c":
			cat = arg
		if opt == "-o":
			obs = str(arg).split(',')
			if 'mjj' not in obs and 'mgg' not in obs:
				print 'Observable must be either mjj or mgg'
				sys.exit(2)
		if opt == "-b":
			arr = str(arg).split(',')
			for b in arr: bins.append(int(b))
		if opt == "-l":
			lumi = arg
		if opt == "-a":
			analysis = arg
		if opt == "-B":
			Blinded = 1
		if opt == "-L":
			myLabel = arg
	if wfile == "" or cat == -1 or obs == "" or lumi == "" or analysis == "":
		print 'MakeBkgPlot.py -w <workspace file> -c <cat> -o <observable> -l <lumi> -a <analysis title>'
		sys.exit(2)


	wroot = TFile(wfile, "READ")
	workspace = wroot.Get("w_all")

	blindedRegions = {}
	blindedRegions['mgg'] = [120, 130]
	blindedRegions['mjj'] = [80, 130]

	CAT = cat
	print CAT
	if int(CAT) == -1:
		CAT = "0"
	for i,ob in enumerate(obs):
		data2D = workspace.data("data_obs_cat"+str(CAT))
		data2D.Print()
		#mjjBkgTmpBer1_cat0_CMS_Bkg_cat0
		ccat = CAT
		if int(CAT) == 2: ccat = 0
		if int(CAT) == 3: ccat = 1
		pdf = workspace.pdf(ob+"BkgTmpBer1_cat"+str(ccat)+"_CMS_Bkg_cat"+str(CAT))
#		pdf = workspace.pdf("BkgPdf_cat"+str(cat))
		var = workspace.var(ob)
		data = data2D.reduce(RooArgSet(var))

		label = "m_{jj} [GeV]"
		if 'mgg' in ob:
			label = "m_{#gamma#gamma} [GeV]"

		MakeBkgPlot(data, pdf, var, label, lumi, cat, analysis, doBands, myLabel+"background_fit_"+ob, bins[i], Blinded)

if __name__ == "__main__":
	main(sys.argv[1:])

