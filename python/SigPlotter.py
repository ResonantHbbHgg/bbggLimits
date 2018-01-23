from ROOT import *
from math import sqrt
from HiggsAnalysis.bbggLimits.NiceColors import *
from HiggsAnalysis.bbggLimits.MyCMSStyle import *

#import RooFit

def getEffSigma(mass, pdf, wmin=110., wmax=130.,  step=0.01, epsilon=1.e-4):
  cdf = pdf.createCdf(RooArgSet(mass))
  point=wmin;
  points = [];
  if wmax > 179: step = 0.1

  while (point <= wmax):
    mass.setVal(point)
    if (pdf.getVal() > epsilon):
      points.append([point,cdf.getVal()]);
    point+=step

  low = wmin;
  high = wmax;
  width = wmax-wmin;

  for i in range(0, len(points)):
    for j in range(i, len(points)): 
      wy = points[j][1] - points[i][1]
      if (abs(wy-0.683) < epsilon): 
        wx = points[j][0] - points[i][0]
        if (wx < width):
          low = points[i][0];
          high = points[j][0];
          width=wx;
  print "effSigma: [", low, "-", high, "] = ", width/2.
  return width/2.


def MakeSigPlot(data, pdf, var, label, lumi, cat, analysis, doBands, fname, binning, Xmin = -1, Xmax = -1, isSig=0, DSCB=0, iHiggs="0", outPath='./'):

	gROOT.SetBatch(kTRUE)

#	pdf.fitTo(data)

	frame = var.frame(RooFit.Title(" "),RooFit.Bins(binning))

	if not isSig: data.plotOn(frame,RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0))
	if isSig: data.plotOn(frame,RooFit.MarkerStyle(kOpenSquare), RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0))

#	SigColor = kAzure
	SigColor = TColor.GetColor(NiceBlueDark)
	if Xmin != -1 and Xmax != -1:
		var.setRange("plot", Xmin, Xmax)
		pdf.plotOn(frame,RooFit.LineColor(SigColor),RooFit.Precision(1E-5), RooFit.Range("plot"))
	else:
		pdf.plotOn(frame,RooFit.LineColor(SigColor),RooFit.Precision(1E-5))
	if not DSCB:
		pdf.plotOn(frame,RooFit.LineColor(kAzure-4),RooFit.LineStyle(kDashDotted),RooFit.Components(str(var.GetName())+"GaussSig_cat"+str(cat)), RooFit.Precision(1E-5))
		pdf.plotOn(frame,RooFit.LineColor(kAzure-9),RooFit.LineStyle(kDashed),RooFit.Components(str(var.GetName())+"CBSig_cat"+str(cat)),RooFit.Precision(1E-5))
        

	curve = frame.getObject( int(1) )
	datah = frame.getObject( int(0) )
	datah.SetLineWidth(1)
	if not DSCB:
	        gauss = frame.getObject( int(2) )
	        cbs = frame.getObject( int(3) )
		cbs.SetLineWidth(2)
#	gauss.SetLineWidth(2)
#	datah.SetMarkerStyle(20)

#	sigmas = MakeBands(data, pdf, var, frame, curve)
#	sigmas[0].SetFillColor(kBlue-5)
#	sigmas[1].SetFillColor(kCyan)

	Max = frame.GetMaximum()
	SetGeneralStyle()
	c = TCanvas("c", "c", 800, 600)
	SetPadStyle(c)
#	c.SetLogy()
	frame.Draw()
	xmax = frame.GetXaxis().GetXmax()
	xmin = frame.GetXaxis().GetXmin()

	deltabin = (xmax - xmin)/binning
	sigmas = MakeBands(data, pdf, var, frame, curve, xmin, xmax, deltabin, isSig)
		
	print xmax, xmin
	sigmas[0].SetFillColor(0)
	sigmas[0].SetLineColor(0)
	sigmas[1].SetFillColor(0)
	sigmas[1].SetLineColor(0)

	sigmas[1].Draw("AE3")
	sigmas[1].SetMaximum(Max*1.1)
	sigmas[1].SetMinimum(0.00001)
	sigmas[1].GetXaxis().SetTitle(label)
	#sigmas[1].GetXaxis().SetTitleSize(0.045)
	#sigmas[1].GetYaxis().SetTitleSize(0.045)
	#sigmas[1].GetXaxis().SetRangeUser(xmin*1.0001, xmax*0.9999)
	sigmas[1].GetYaxis().SetTitle("Normalization (arbitrary units)")
	c.Update()
	SetAxisTextSizes(sigmas[1], 0.05, 1.05)
	c.Update()
	if Xmin > 0 and Xmax > 0:
		print Xmin, Xmax
		sigmas[1].GetXaxis().SetLimits(Xmin, Xmax)
		sigmas[1].GetXaxis().SetRangeUser(Xmin, Xmax)

	c.Update()
	sigmas[0].Draw("E3same")
	frame.Draw("same")
	
	datah.Draw("EPsame")

	tlatex = TLatex()
	tlatex.SetNDC()
	tlatex.SetTextAngle(0)
	tlatex.SetTextColor(kBlack)
	tlatex.SetTextFont(63)
	tlatex.SetTextAlign(11)
	tlatex.SetTextSize(25)
	tlatex.SetTextFont(43)
	tlatex.SetTextSize(20)
	tlatex.SetTextSize(23)
	xbegin = 0.18
	ybegin = 0.91
	Cat = "High-purity Category"
	if int(cat) == 1 or int(cat) == 3:
		Cat = "Medium-purity Category"
	if int(cat) == -1:
		Cat = "High-mass (Single Cat.)"
	print cat, Cat
	if "|" in analysis:
		an = analysis.split("|")
		for iAN,AN in enumerate(an):
                  #			tlatex.DrawLatex(xbegin, ybegin-0.06*float(iAN), AN)
			tlatex.DrawLatex(0.60, ybegin-0.055*float(iAN)-0.1, AN)
                        #		tlatex.SetTextFont(63)
                        #		tlatex.DrawLatex(xbegin, ybegin, an[0])
                        #		tlatex.SetTextFont(43)
                        #		tlatex.DrawLatex(xbegin, ybegin-0.06, an[1])	
                        #		tlatex.DrawLatex(xbegin, ybegin-0.06*float(len(an)), Cat)
		tlatex.DrawLatex(0.60, ybegin-0.06*float(len(an))-0.1, Cat)
	else:
		tlatex.DrawLatex(0.65, ybegin-0.1, analysis)
		tlatex.DrawLatex(0.65, ybegin-0.16, Cat)

	leg = TLegend(0.60, ybegin-float(len(analysis.split("|"))+2)*0.06-0.33, 0.85, ybegin-float(len(analysis.split("|"))+2)*0.06)
	if '|' not in analysis:
		leg =  TLegend(0.7, ybegin-0.61, 0.935, ybegin-0.21)

	leg.SetFillStyle(0)
	leg.SetLineWidth(0)
	leg.SetBorderSize(0)

	nBkgParams = pdf.getParameters(data).getSize()
	print "Number of background parameters:", nBkgParams
	pdf.getParameters(data).Print("v")

        bkgModel = "Signal model"
	
	leg.AddEntry(datah, "Signal Simulation", "pe")
	leg.AddEntry(curve, bkgModel, "l")
	if not DSCB:
		leg.AddEntry(gauss, "Gaussian component", "l")
		leg.AddEntry(cbs, "Crystal Ball component", "l")
	meanName = str(var.GetName()) + "_sig_m0_cat"+str(cat)
	thisMean = pdf.getParameters(data).getRealValue(meanName)
	if iHiggs  != "0":
		meanName = str(var.GetName()) + "_hig_m0_"+iHiggs+"_cat"+str(cat)
		thisMean = pdf.getParameters(data).getRealValue(meanName)
	if thisMean != 0:
		leg.AddEntry(sigmas[0], "#mu = "+ str("%.1f" % thisMean)+ " GeV", "l")
		leg.AddEntry(sigmas[1], "#sigma_{Eff} = "+str("%.1f" % getEffSigma(var, pdf, Xmin, Xmax)) + " GeV", "l")
	leg.Draw()
	DrawCMSLabels(c, '')
	c.SaveAs(outPath+'/'+fname+".pdf")
	c.SaveAs(outPath+'/'+fname+".png")

def MakeBands(data, pdf, var, frame, curve, xmin, xmax, deltabin, isSig):

	onesigma = TGraphAsymmErrors()
	twosigma = TGraphAsymmErrors()

	bins = []
	bins.append([xmin*1.001, xmin, xmin+deltabin])
	for ibin in range(1,frame.GetXaxis().GetNbins()+1):
		lowedge = frame.GetXaxis().GetBinLowEdge(ibin)
		upedge  = frame.GetXaxis().GetBinUpEdge(ibin)
		center  = frame.GetXaxis().GetBinCenter(ibin)
		bins.append(  (center,lowedge,upedge) )
	bins.append([xmax*0.999, xmax-deltabin, xmax])


	allbins = []
	for ibin,bin in enumerate(bins):
		center,lowedge,upedge = bin

		nombkg = curve.interpolate(center)
		largeNum = nombkg*50#(1 + sqrt(nombkg)*2)
#		if nombkg < 1:
#			largeNum = 9999999999999999999999999999999999
		largeNum = max(1,largeNum)

		nlim = RooRealVar("nlim%s" % var.GetName(),"",0.,-largeNum,largeNum)
#		nlim.removeRange()

		onesigma.SetPoint(ibin,center,nombkg)
		twosigma.SetPoint(ibin,center,nombkg)
		if isSig: continue

		nlim.setVal(nombkg)

		if 1:
			print "computing error band ", ibin, lowedge, upedge, nombkg, 

		var.setRange("errRange",lowedge,upedge)
		errm = 0
		errp = 0
		counter = 0
#		while errm == 0 or errp == 0:
		errmSum = 0
		countermSum = 0
		errpSum = 0
		counterpSum = 0
		epdf = RooExtendPdf("epdf","",pdf,nlim, "errRange")
		nll = epdf.createNLL(data,RooFit.Extended())
		minim = RooMinimizer(nll)
		for ct in range(0, 500):
#			epdf = RooExtendPdf("epdf","",pdf,nlim, "errRange")
#			nll = epdf.createNLL(data,RooFit.Extended())
	#		nll = pdf.createNLL(data)
			minim = RooMinimizer(nll)
			minim.setMinimizerType("Minuit2")
			minim.setStrategy(0)
			minim.setPrintLevel( -1 )#if not options.verbose else 2)
			minim.migrad()
#			minim.hesse()
#			minim.migrad()
#			minim.minos()
			minim.setStrategy(0)
			minim.minos(RooArgSet(nlim))

			errm, errp = -nlim.getErrorLo(),nlim.getErrorHi()
#			if errm !=0 and errp !=0: break

			if abs(errm) > nombkg*0.01:
				errmSum += -nlim.getErrorLo()
				countermSum+=1
#				if len(countermSum) < 3: countermSum.append(-nlim.getErrorLo())
			if abs(errp) > nombkg*0.01:
				errpSum += nlim.getErrorHi()
				counterpSum+=1
#				if len(counterpSum) < 3: counterpSum.append(nlim.getErrorHi())
#			del minim
#			del nll
#			del epdf
			if counterpSum > 10 and countermSum > 10: break

#			if len(counterpSum) > 2 and len(countermSum) > 2: break
#			counter +=1
#			if counter > 10: break



		errm = errmSum/float(countermSum)
		errp = errpSum/float(counterpSum)

#		if len(countermSum) == 3:
#			errm = sorted(countermSum)[1]
#		else:
#			errm = 0
#		if len(counterpSum) == 3:
#			errp = sorted(counterpSum)[1]
#		else:
#			errp = 0

		onesigma.SetPointError(ibin,center-lowedge,upedge-center,errm,errp)

		errm = 0
		errp = 0
		counter = 0
		errmSum = 0
		countermSum = []
		errpSum = 0
		counterpSum = []
#		while errm == 0 or errp == 0:
		for ct in range(0,2000):
#			epdf = RooExtendPdf("epdf","",pdf,nlim, "errRange")
#			nll = epdf.createNLL(data,RooFit.Extended())
	#		nll = pdf.createNLL(data)
#			minim = RooMinimizer(nll)
			minim.setMinimizerType("Minuit2")
#			minim.setStrategy(2)
			minim.setPrintLevel( -1 )#if not options.verbose else 2)
			minim.setErrorLevel(2.)
			minim.migrad()
			minim.minos(RooArgSet(nlim))

			errm, errp = -nlim.getErrorLo(),nlim.getErrorHi()
#			if errm != 0 and errp !=0: break

			if abs(errm) > nombkg*0.01:
				errmSum += -nlim.getErrorLo()
				if len(countermSum) < 3: countermSum.append(-nlim.getErrorLo())
			if abs(errp) > nombkg*0.01:
				errpSum += nlim.getErrorHi()
				if len(counterpSum) < 3: counterpSum.append(nlim.getErrorHi())

#			del epdf
#			del nll
#			del minim

			print len(counterpSum)
			if len(counterpSum) > 2 and len(countermSum) > 2: break
#			counter += 1
#			if counter > 10: break

#		errm = errmSum/float(countermSum)
#		errp = errpSum/float(counterpSum)
		errm = sorted(countermSum)[1]
		errp = sorted(counterpSum)[1]
		allbins.append([errm,errp,nombkg])
		twosigma.SetPointError(ibin,center-lowedge,upedge-center,errm,errp)
		
#		del minim
#		del nll

	print allbins
	return [onesigma,twosigma]
