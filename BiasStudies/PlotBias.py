from ROOT import *

nevs = [10, 50, 75, 100, 150, 500, 1500]
pts = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10']
lbs = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

fbase = 'mlfit_alpha_PT_t1000_eNV.root'

out = TFile('bias_out.root', 'RECREATE')
grs = []
for n in nevs:
  gr = TGraphErrors()
  gr.SetName('gr_err_'+str(n))
  for ip,p in enumerate(pts):
    tf = TFile(fbase.replace('PT', p).replace('NV', str(n)))
    tt = tf.Get('tree_fit_sb')
    hname = 'hist'+str(n)+str(p)
    hist = TH1F(hname, '', 100, -2, 2)
    tt.Draw('(mu-1)/muErr>>'+hname, 'numbadnll>-1','goff')
    fitres = hist.Fit('gaus', 'S')
    mean = fitres.GetParams()[1]
    e_mean = fitres.GetErrors()[1]
    hname2 = 'hist'+str(n)+str(p)+'2'
    hist2 = TH1F(hname2, '', 100, -5, 5)
    tt.Draw('(mu-1)/muErr>>'+hname2, 'numbadnll>-1','goff')
    fitres2 = hist2.Fit('gaus', 'S')
    mean2 = fitres2.GetParams()[1]
    e_mean2 = fitres2.GetErrors()[1]
    mymean = mean
    myerr = e_mean
    if abs(mean2) < abs(mean):
      mymean = mean2
      myerr = e_mean2
    print mymean, myerr
    gr.SetPoint(ip, lbs[ip], mymean)
    gr.SetPointError(ip, float(0), myerr)
  grs.append(gr)

out.cd()
for igr,gr in enumerate(grs):
  gr.Draw("A")
  gr.GetXaxis().SetRangeUser(-0.1, 1.1)
  gr.GetXaxis().SetLimits(-0.1, 1.1)
  gr.SetLineColor(igr+1)
  gr.SetMarkerColor(igr+1)
  gr.Write()
out.Close()
