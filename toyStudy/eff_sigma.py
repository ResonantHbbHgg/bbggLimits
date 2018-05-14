from ROOT import *
gROOT.SetBatch()

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
  # print "effSigma: [", low, "-", high, "] = ", width/2.
  return width/2., low, high


if __name__ == "__main__":
      
  rooWsSig = TFile('hhbbgg.mH125_13TeV.inputsig.root')
  
  sigWs = rooWsSig.Get('w_all')
  # sigWs.Print()
  
  mgg = sigWs.var('mgg')
  
  pdfSig_mgg  = sigWs.pdf("mggSig_cat0_CMS_sig_cat2")
  
  print getEffSigma(mgg, pdfSig_mgg, 110, 170)
