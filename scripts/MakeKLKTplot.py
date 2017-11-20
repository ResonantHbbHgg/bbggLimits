import argparse, os, json, sys
from HiggsAnalysis.bbggLimits.Utils2DScan import *
from ROOT import *
from array import array
from HiggsAnalysis.bbggLimits.NiceColors import *
from HiggsAnalysis.bbggLimits.MyCMSStyle import *

gROOT.SetBatch()

parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument("-u", '--unblind', dest='unblind', action='store_true', default=False)
parser.add_argument("-l", '--limitsFile', dest="lims", type=str, nargs='+')
parser.add_argument('-o', '--outFile', dest='outf', type=str, default='KLKT')
parser.add_argument('--funkyColor', dest='funky', action='store_true', default=False)
opt = parser.parse_args()

transparency = 0.6
transparency_obs = transparency

myGraphs = []
for ff in opt.lims:
  if 'txt' not in ff:
    print 'Unknown file format, please specify either txt or json'
    sys.exit()

myGraphs =  readTxt(opt.lims)

outFile = opt.lims[0].split('/')[0]+'/'+opt.outf
outf = TFile(outFile+'.root', 'RECREATE')

outf.cd()
graphs = {}
contours = {}
xmin = -20
xmax = 20
ymin = -2.4
ymax = 2.4
npx = 50
npy = npx
for igr,gr in enumerate(myGraphs):
  if igr == len(myGraphs)-1: break
  graphs[gr.GetName()] = gr
  if 'diff' in str(gr.GetName()):
    intergr = InterpolatedGraph(gr, xmin, xmax, ymin, ymax, npx, npy)
    intergr.Write()
#    myConts = GetContours(gr)
    myConts = GetContours(intergr)
#    contours[gr.GetName()] = []
    contours[intergr.GetName()] = []
    for mc in myConts:
      mc.Write()
#      contours[gr.GetName()].append(mc)
      contours[intergr.GetName()].append(mc)
  if 'cen' in str(gr.GetName()):
    intercen = InterpolatedGraph(gr, xmin, xmax, ymin, ymax, npx, npy)
    intercen.Write()
  if 'obs' in str(gr.GetName()):
    interobs = InterpolatedGraph(gr, xmin, xmax, ymin, ymax, npx, npy)
    interobs.Write()
  gr.Write()

print contours
#1s bands:
print 'Splitting bands'
bsplit_1s = SplitBands('1s', contours)
bsplit_2s = SplitBands('2s', contours)

print bsplit_1s
print bsplit_2s

print 'Making bands'
Points = []
stepx = float(xmax-xmin)/float(npx)
stepy = float(ymax-ymin)/float(npy)
for xx in xrange(0, npx+1):
  for yy in xrange(0, npy+1):
    Points.append([xmin+xx*stepx, ymin+yy*stepx])
#Points = myGraphs[len(myGraphs)-1]
band_1s_low = MakeBands(bsplit_1s[1], bsplit_1s[3], Points, -2.5, bname='band_1s_low')
band_1s_upp = MakeBands(bsplit_1s[0], bsplit_1s[2], Points, 2.5, bname='band_1s_upp')
band_2s_low = MakeBands(bsplit_2s[1], bsplit_2s[3], Points, -2.5, bname='band_2s_low')
band_2s_upp = MakeBands(bsplit_2s[0], bsplit_2s[2], Points, 2.5, bname='band_2s_upp')

#Do observed bands
obsSplit = SplitBands('obs', contours)
obs_upper0 = TGraph()
obs_lower0 = TGraph()
for xx in xrange(0,npx+1):
  X = xmin+xx*stepx
  obs_upper0.SetPoint(xx, X, ymax)
  obs_lower0.SetPoint(xx, X, ymin)
band_obs_upp = MakeBands([obs_upper0], obsSplit[0], Points, 2.5, bname='band_obs_upp')
band_obs_low = MakeBands(obsSplit[1], [obs_lower0], Points, -2.5, bname='band_obs_low')
#band_obs_upp.SetLineColor(TColor.GetColor(NiceBlue))
#band_obs_low.SetLineColor(TColor.GetColor(NiceBlue))
band_obs_upp.SetLineColor(kBlack)
band_obs_low.SetLineColor(kBlack)
band_obs_upp.SetLineWidth(1)
band_obs_low.SetLineWidth(1)
band_obs_upp.SetFillStyle(1001)
band_obs_low.SetFillStyle(1001)
obsCol = TColor.GetColor(NiceBlue)
if opt.funky: obsCol = 6
band_obs_upp.SetFillColorAlpha(obsCol, transparency_obs)
band_obs_low.SetFillColorAlpha(obsCol, transparency_obs)

band_1s_low.Write()
band_1s_upp.Write()
band_2s_low.Write()
band_2s_upp.Write()

band_2s_low.SetFillStyle(1001)
band_2s_low.SetFillColorAlpha(kGray+1, transparency)
band_2s_upp.SetFillStyle(1001)
band_2s_upp.SetFillColorAlpha(kGray+1, transparency)
band_1s_low.SetFillStyle(1001)
band_1s_low.SetFillColorAlpha(kGray+2, transparency)
band_1s_upp.SetFillStyle(1001)
band_1s_upp.SetFillColorAlpha(kGray+2, transparency)

gr_th = MakeThPlot(myGraphs[len(myGraphs)-1])
gr_th.Write()

xsfunc = TF2('xsfunc', functionGF_klkt_Norm, -20., 20., -2.5, 2.5)
isolines = [0.01, 0.1, 1.0, 2., 5., 10., 20., 50.]
xsfunc.SetContour(len(isolines), array('d', isolines))
xsfunc.SetNpx(250)
xsfunc.SetNpy(250)
xsfunc.SetLineColor(kBlack)
xsfunc.SetLineStyle(3)
xsfunc.SetLineWidth(1)

hframe = TH2F("hframe", ";#kappa_{#lambda};#kappa_{t}", 100, xmin, xmax, 100, ymin, ymax)
pt_sm = TGraph()
pt_sm.SetPoint(0, 1, 1)
pt_sm.SetMarkerStyle(33)
pt_sm.SetMarkerSize(3)
pt_sm.SetMarkerColor(TColor.GetColor(NiceRed))

leg = TLegend(0.1, 0.05, 0.9, 0.17)
leg.SetNColumns(3)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetHeader('95% CL upper limits')

#leg.SetHeader('c_{g} = c_{2} = c_{2g} = 0')

SetGeneralStyle()
c1 = TCanvas("c1", "c1", 900, 1000)
gStyle.SetOptStat(0)
hframe.Draw()
c1.SetBottomMargin(0.25)
c1.SetTopMargin(-0.1)
SetPadStyle(c1)
c1.Update()
SetAxisTextSizes(hframe)

lat = TLatex()
lat.SetNDC()
lat.SetTextAlign(11)
lat.SetTextSize(0.03)
leg.SetTextSize(0.027)
lat.DrawLatex(0.15,0.182,'pp#rightarrowHH#rightarrow#gamma#gammab#bar{b}')

lat.SetTextSize(0.025)
lat.DrawLatex(0.50,0.182,'c_{g} = c_{2} = c_{2g} = 0')


lat.SetTextAngle(-45)
lat.SetTextSize(0.015)
lat.SetTextColor(kGray+1)
#isolines = [0.01, 0.1, 1.0, 2., 5., 10., 20., 50.]
prop = 0.12
stpx = (0.9-0.1)/(xmax - xmin)
stpy = (0.9-0.25)/(4.8)
x0 = 0.1
y0 = 0.25
lat.SetTextAlign(22)
lat.DrawLatex(x0+stpx*(2.2+20), y0 + stpy*(2.2+20)*prop, '0.01 fb')
lat.DrawLatex(x0+stpx*(4.4+20), y0 + stpy*(4.4+20)*prop, '0.1 fb')
lat.DrawLatex(x0+stpx*(8.22+20), y0 + stpy*(8.22+20)*prop, '1 fb')
lat.DrawLatex(x0+stpx*(9.82+20), y0 + stpy*(9.82+20)*prop, '2 fb')
lat.DrawLatex(x0+stpx*(12.3+20), y0 + stpy*(12.3+20)*prop, '5 fb')
lat.DrawLatex(x0+stpx*(14.7+20), y0 + stpy*(14.7+20)*prop, '10 fb')
lat.DrawLatex(x0+stpx*(17.6+20), y0 + stpy*(17.6+20)*prop, '20 fb')
lat.DrawLatex(x0+stpx*(-2.2+20), y0 + stpy*(-2.2+20)*prop, '0.01 fb')
lat.DrawLatex(x0+stpx*(-4.4+20), y0 + stpy*(-4.4+20)*prop, '0.1 fb')
lat.DrawLatex(x0+stpx*(-8.22+20), y0 + stpy*(-8.22+20)*prop, '1 fb')
lat.DrawLatex(x0+stpx*(-9.82+20), y0 + stpy*(-9.82+20)*prop, '2 fb')
lat.DrawLatex(x0+stpx*(-12.3+20), y0 + stpy*(-12.3+20)*prop, '5 fb')
lat.DrawLatex(x0+stpx*(-14.7+20), y0 + stpy*(-14.7+20)*prop, '10 fb')
lat.DrawLatex(x0+stpx*(-17.6+20), y0 + stpy*(-17.6+20)*prop, '20 fb')

#gr_th.Draw("A COLZ")
xsfunc.Draw('cont3 same')
band_2s_low.Draw('same F')
band_2s_upp.Draw('same F')
band_1s_low.Draw('same F')
band_1s_upp.Draw('same F')
for cc in contours:
  if 'cen' not in cc: continue
  for iccc,ccc in enumerate(contours[cc]):
    ccc.SetLineColor(kBlack)
    ccc.SetLineStyle(kDashed)
    ccc.SetLineWidth(3)
    ccc.Draw('same C')
    if iccc == 0: leg.AddEntry(ccc, 'Expected', 'l')
leg.AddEntry(band_1s_low, 'Expected #pm1 std. deviation', 'f')
leg.AddEntry(pt_sm, 'SM', 'p')
pt_sm.Draw('p same')
if opt.unblind:
  band_obs_upp.Draw('same F')
  band_obs_low.Draw('same F')
  for cc in contours:
    if 'obs' not in cc: continue
    for ccc in contours[cc]:
      ccc.SetLineColor(kBlack)
      ccc.SetLineWidth(2)
      ccc.Draw('same C')
leg.AddEntry(band_obs_upp, 'Observed', 'f')
leg.AddEntry(band_2s_low, 'Expected #pm2 std. deviation', 'f')


c1.Update()
c1.cd()
c1.Update()
c1.RedrawAxis()

leg.Draw('same')

DrawCMSLabels(c1, '35.9')
c1.SaveAs(outFile+"_plot.pdf")
c1.SaveAs(outFile+"_plot.png")


c0 = TCanvas('c0', 'c0', 800, 600)
set_palette()
#hframe.Draw()
#hframe.GetZaxis().SetTitle('95% C.L. Expected Limit [fb]')
#SetAxisTextSizes(hframe)
c0.SetRightMargin(0.2)
c0.Update()
for gr,gr in enumerate(myGraphs):
  if 'cen' in str(gr.GetName()):
    gr.SetTitle('')
    gr.Draw('COLZ')
    c0.Update()
    gr.GetXaxis().SetTitle('#kappa_{#lambda}')
    gr.GetYaxis().SetTitle('#kappa_{t}')
    gr.GetZaxis().SetTitle('Expected 95% upper limit [fb]')
    gr.GetZaxis().SetRangeUser(0., 5.)
    c0.Update()
    SetAxisTextSizes(gr)
    c0.Update()
    SetPadStyle(c0)
    c0.Update()
    DrawCMSLabels(c0, '35.9')
    c0.Update()
    c0.SaveAs(outFile+'_cen_lims.pdf')
    c0.SaveAs(outFile+'_cen_lims.png')
    break

if opt.unblind == True:
  cob = TCanvas('cob', 'cob', 800, 600)
  set_palette()
  #hframe.Draw()
  #hframe.GetZaxis().SetTitle('95% C.L. Expected Limit [fb]')
  #SetAxisTextSizes(hframe)
  cob.SetRightMargin(0.2)
  cob.Update()
  for gr,gr in enumerate(myGraphs):
    if 'obs' in str(gr.GetName()):
      gr.SetTitle('')
      gr.Draw('COLZ')
      cob.Update()
      gr.GetXaxis().SetTitle('#kappa_{#lambda}')
      gr.GetYaxis().SetTitle('#kappa_{t}')
      gr.GetZaxis().SetTitle('Observed 95% upper limit [fb]')
      gr.GetZaxis().SetRangeUser(0., 5.)
      cob.Update()
      SetAxisTextSizes(gr)
      cob.Update()
      SetPadStyle(cob)
      cob.Update()
      DrawCMSLabels(cob, '35.9')
      cob.Update()
      cob.SaveAs(outFile+'_obs_lims.pdf')
      cob.SaveAs(outFile+'_obs_lims.png')
      break

outf.Close()
