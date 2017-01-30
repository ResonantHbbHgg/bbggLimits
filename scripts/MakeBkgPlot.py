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
	bins = []
	try:
		opts, args = getopt.getopt(argv,"w:c:o:l:a:b:",["workspace=", "cat=", "observable=","lumi=","analysis=","bins="])
	except getopt.GetoptError:
		print 'MakeBkgPlot.py -w <workspace file> -c <cat> -o <obs1,obs2,obs3> -l <lumi> -a <analysis title> -b <binobs1,binobs2>'
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
	if wfile == "" or cat == -1 or obs == "" or lumi == "" or analysis == "":
		print 'MakeBkgPlot.py -w <workspace file> -c <cat> -o <observable> -l <lumi> -a <analysis title>'
		sys.exit(2)


	wroot = TFile(wfile, "READ")
	workspace = wroot.Get("w_all")

	CAT = cat
	print CAT
	if int(CAT) == -1:
		CAT = "0"
	for i,ob in enumerate(obs):
		data2D = workspace.data("data_obs_cat"+str(CAT))
		data2D.Print()
		pdf = workspace.pdf(ob+"BkgTmpBer1_cat"+str(CAT))
#		pdf = workspace.pdf("BkgPdf_cat"+str(cat))
		var = workspace.var(ob)
		data = data2D.reduce(RooArgSet(var))

		label = "M(jj) [GeV]"
		if 'mgg' in ob:
			label = "M(#gamma#gamma) [GeV]"	

		MakeBkgPlot(data, pdf, var, label, lumi, cat, analysis, doBands, "background_fit_"+ob, bins[i])

if __name__ == "__main__":
	main(sys.argv[1:])

