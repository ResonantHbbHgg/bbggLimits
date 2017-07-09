#!/usr/bin/env python

from ROOT import *
from HiggsAnalysis.bbggLimits.SigPlotter import *
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
	Label = ""
        DSCB = False
        iHiggs = 'ggh'
	xmax = {'mgg':140, 'mjj': 190}
        xmin = {'mgg':118, 'mjj': 70}
	try:
		opts, args = getopt.getopt(argv,"w:c:o:l:a:b:L:i:D",["workspace=", "cat=", "observable=","lumi=","analysis=","bins=", "Label=", "iHiggs", "DSCB"])
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
		if opt == "-L":
			Label = arg
		if opt == "--DSCB":
			DSCB = True
		if opt == "-i" or opt == "--iHiggs":
			iHiggs = str(arg)
	if wfile == "" or cat == -1 or obs == "" or lumi == "" or analysis == "":
		print 'MakeBkgPlot.py -w <workspace file> -c <cat> -o <observable> -l <lumi> -a <analysis title>'
		sys.exit(2)


	wroot = TFile(wfile, "READ")
	workspace = wroot.Get("w_all")

	CAT = cat
	print CAT
	if int(CAT) == -1:
		CAT = "0"
	ccat = CAT
	if int(CAT) == 2: ccat = 0
	if int(CAT) == 3: ccat = 1
        
	for i,ob in enumerate(obs):
		print "Dataset name:", "Hig_"+iHiggs+"_cat"+str(CAT)
		data2D = workspace.data("Hig_"+iHiggs+"_cat"+str(CAT))
		data2D.Print()
		print "##############################################"
		print "######## Number of weighted events:", "Hig_"+iHiggs+"_cat"+str(CAT), Label, data2D.sumEntries()
		print "##############################################"
		print 'mggHig_ggh_cat0_CMS_hig_ggh_cat2'
		pdfname = ob+"Hig_"+iHiggs+"_cat"+str(ccat)+"_CMS_hig_"+iHiggs+'_cat'+str(CAT)
                print 'mjjHig_ggh_cat0_CMS_hig_ggh_cat2'
		print pdfname
#		workspace.Print()
		pdf = workspace.pdf(pdfname)#ob+"Hig_"+iHiggs+"_cat"+str(ccat)+"_CMS_hig_"+iHiggs+'cat'+str(CAT))
#		print ob+"Hig_"+iHiggs+"_cat"+str(CAT)
		pdf.Print()
#		pdf = workspace.pdf("BkgPdf_cat"+str(cat))
		var = workspace.var(ob)
		data = data2D.reduce(RooArgSet(var))

		label = "M(jj) [GeV]"
		if 'mgg' in ob:
			label = "M(#gamma#gamma) [GeV]"	

		MakeSigPlot(data, pdf, var, label, lumi, cat, analysis, doBands, Label+"_signal_fit_"+ob+"_cat"+str(CAT), bins[i], xmin[ob], xmax[ob], 1, DSCB, iHiggs)

if __name__ == "__main__":
	main(sys.argv[1:])

