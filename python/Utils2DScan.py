import argparse, os, json, sys
from array import array
from ROOT import *

gROOT.SetBatch()

A13tev = [2.09078, 10.1517, 0.282307, 0.101205, 1.33191, -8.51168, -1.37309, 2.82636, 1.45767, -4.91761, -0.675197, 1.86189, 0.321422, -0.836276, -0.568156]
HHbbgg_xsec_br = 33.49*0.0026

def functionGF(kl,kt,c2,cg,c2g,A): 
  if abs(kt) < 0.00001: kt = 0.00001
  xsec =  A[0]*kt**4 + A[1]*c2**2 + (A[2]*kt**2 + A[3]*cg**2)*kl**2  + A[4]*c2g**2 + ( A[5]*c2 + A[6]*kt*kl )*kt**2  + (A[7]*kt*kl + A[8]*cg*kl )*c2 + A[9]*c2*c2g  + (A[10]*cg*kl + A[11]*c2g)*kt**2+ (A[12]*kl*cg + A[13]*c2g )*kt*kl + A[14]*cg*c2g*kl
  return xsec

def functionGF_klkt_Norm(x): 
  xsec = functionGF(x[0], x[1], 0.0, 0.0, 0.0, A13tev)
  return xsec*HHbbgg_xsec_br

def readTxt(lims):
  gr2d_cen = TGraph2D()
  gr2d_cen.SetName('gr2d_cen')
  gr2d_obs = TGraph2D()
  gr2d_obs.SetName('gr2d_obs')
  gr2d_m1s = TGraph2D()
  gr2d_m1s.SetName('gr2d_m1s')
  gr2d_p1s = TGraph2D()
  gr2d_p1s.SetName('gr2d_p1s')
  gr2d_m2s = TGraph2D()
  gr2d_m2s.SetName('gr2d_m2s')
  gr2d_p2s = TGraph2D()
  gr2d_p2s.SetName('gr2d_p2s')
  gr2d_cen_diff = TGraph2D()
  gr2d_cen_diff.SetName('gr2d_cen_diff')
  gr2d_obs_diff = TGraph2D()
  gr2d_obs_diff.SetName('gr2d_obs_diff')
  gr2d_m1s_diff = TGraph2D()
  gr2d_m1s_diff.SetName('gr2d_m1s_diff')
  gr2d_p1s_diff = TGraph2D()
  gr2d_p1s_diff.SetName('gr2d_p1s_diff')
  gr2d_m2s_diff = TGraph2D()
  gr2d_m2s_diff.SetName('gr2d_m2s_diff')
  gr2d_p2s_diff = TGraph2D()
  gr2d_p2s_diff.SetName('gr2d_p2s_diff')
  allPoints = []

  counter = 0
  for jfile in lims:
    print 'reading file', jfile
    ffile = open(jfile, 'r')
    for ll in ffile:
      ll_spl = ll.split()
      if len(ll_spl) < 9: continue
      ipt = int(ll_spl[0])
      kl = float(ll_spl[1])
      kt = float(ll_spl[2])
      exp = float(ll_spl[3])#/xsec#(HHbbgg_xsec_br*functionGF(kl,kt,0.0,0.0,0.0,A13tev))
      obs = float(ll_spl[4])#/xsec#(HHbbgg_xsec_br*functionGF(kl,kt,0.0,0.0,0.0,A13tev))
      p1s = float(ll_spl[5])#/xsec#(HHbbgg_xsec_br*functionGF(kl,kt,0.0,0.0,0.0,A13tev))
      m1s = float(ll_spl[6])#/xsec#(HHbbgg_xsec_br*functionGF(kl,kt,0.0,0.0,0.0,A13tev))
      p2s = float(ll_spl[7])#/xsec#(HHbbgg_xsec_br*functionGF(kl,kt,0.0,0.0,0.0,A13tev))
      m2s = float(ll_spl[8])#/xsec#(HHbbgg_xsec_br*functionGF(kl,kt,0.0,0.0,0.0,A13tev))
      thisPoint = [kl,kt]
      if thisPoint in allPoints:
        print 'Point already added!', thisPoint
        continue
      else:
        ipt = counter
        gr2d_cen.SetPoint(ipt, kl, kt, exp)
        gr2d_obs.SetPoint(ipt, kl, kt, obs)
        gr2d_m1s.SetPoint(ipt, kl, kt, m1s)
        gr2d_p1s.SetPoint(ipt, kl, kt, p1s)
        gr2d_m2s.SetPoint(ipt, kl, kt, m2s)
        gr2d_p2s.SetPoint(ipt, kl, kt, p2s)
        xsec = float(HHbbgg_xsec_br*functionGF(kl,kt,0.0,0.0,0.0,A13tev))
        gr2d_cen_diff.SetPoint(ipt, kl, kt, exp-xsec)
        gr2d_obs_diff.SetPoint(ipt, kl, kt, obs-xsec)
        gr2d_m1s_diff.SetPoint(ipt, kl, kt, m1s-xsec)
        gr2d_p1s_diff.SetPoint(ipt, kl, kt, p1s-xsec)
        gr2d_m2s_diff.SetPoint(ipt, kl, kt, m2s-xsec)
        gr2d_p2s_diff.SetPoint(ipt, kl, kt, p2s-xsec)
        allPoints.append(thisPoint)
        counter += 1
    ffile.close()
  return [gr2d_cen, gr2d_obs, gr2d_m1s, gr2d_p1s, gr2d_m2s, gr2d_p2s, gr2d_cen_diff, gr2d_obs_diff, gr2d_m1s_diff, gr2d_p1s_diff, gr2d_m2s_diff, gr2d_p2s_diff, allPoints]
    

def MakeThPlot(pts):
  th = TGraph2D()
  th.SetName('HHbbgg_xsec')
  for ipt,pt in enumerate(pts):
    kl = pt[0]
    kt = pt[1]
    thpt = HHbbgg_xsec_br*functionGF(kl,kt,0.0,0.0,0.0,A13tev)
    th.SetPoint(ipt, kl, kt, thpt)
  return th

def InterpolatedGraph(tgr, xmin, xmax, ymin, ymax, npx, npy):
  inter = TGraph2D()
  inter.SetName(tgr.GetName() + '_inter')
  stpx = float(xmax-xmin)/float(npx)
  stpy = float(ymax-ymin)/float(npy)
  counter = 0
  for xx in xrange(0,npx+1):
    for yy in xrange(0,npy+1):
      X = xmin + xx*stpx
      Y = ymin + yy*stpy
      Z = tgr.Interpolate(X,Y)
      inter.SetPoint(counter, X, Y, Z)
      counter += 1
  return inter

def SplitBands(name, contours):
  uhalf_m = []
  lhalf_m = []
  uhalf_p = []
  lhalf_p = []
  for c in contours:
    if name not in c: continue
    for cc in contours[c]:
      ys_org = cc.GetY()
      lencc = cc.GetN()
      ys =  [ys_org[i] for i in xrange(0, lencc)]
      ytot = 0
      for yy in ys: ytot += yy
      if ytot > 0 and str('m'+name) in c: uhalf_m.append(cc)
      elif ytot > 0 and str('p'+name) in c: uhalf_p.append(cc)
      elif ytot < 0 and str('m'+name) in c: lhalf_m.append(cc)
      elif ytot < 0 and str('p'+name) in c: lhalf_p.append(cc)
      elif ytot > 0: uhalf_m.append(cc)
      elif ytot < 0: lhalf_m.append(cc)
  return [ uhalf_m, lhalf_m, uhalf_p, lhalf_p ]

def GetContours(tgr, level=0):
  print tgr
  cc = TCanvas('cc', 'cc', 800, 600)
  tgr.Draw("ACONT5")
  hist = tgr.GetHistogram()
  gpt = hist.GetPainter()
  contours = gpt.GetContourList(level)
  myConts = []
  for ii in xrange(0, contours.GetEntries()):
#    print contours.At(ii)
    tCont = contours.At(ii)
    tCont.SetName(tgr.GetName() + '_cont_'+str(ii))
    myConts.append(tCont)
  return myConts

def MakeBands(upper, lower, pts, default=0, bname='band'):
  allkl_pts = []
  for pt in pts:
    if pt[0] not in allkl_pts: allkl_pts.append(pt[0])

  allkl_plots_upp = []
  allkl_plots_low = []
  ally_plots_upp = []
  ally_plots_low = []
  for cont in upper:
    ncont = cont.GetN()
    allkl_plots_upp += [cont.GetX()[ic] for ic in xrange(0, ncont)]
    ally_plots_upp += [cont.GetY()[ic] for ic in xrange(0, ncont)]
  for cont in lower:
    ncont = cont.GetN()
    allkl_plots_low += [cont.GetX()[ic] for ic in xrange(0, ncont)]
    ally_plots_low += [cont.GetY()[ic] for ic in xrange(0, ncont)]

  ally_upp = []
  ally_low = []
  for kl in allkl_pts:
    if kl in allkl_plots_upp:
      ind = allkl_plots_upp.index(kl)
      ally_upp.append(ally_plots_upp[ind])
    else:
      ally_upp.append(default)
    if kl in allkl_plots_low:
      ind = allkl_plots_low.index(kl)
      ally_low.append(ally_plots_low[ind])
    else:
      ally_low.append(default)

  nx = len(allkl_pts)
  band = TGraph(nx*2)
  band.SetName(bname)
  for ix,x in enumerate(allkl_pts):
    band.SetPoint(ix, x, ally_upp[ix])
    band.SetPoint(nx+ix, allkl_pts[nx - ix - 1], ally_low[nx - ix - 1])
  return band


def set_palette(name='blue', ncontours=999):
    """Set a color palette from a given RGB list
    stops, red, green and blue should all be lists of the same length
    see set_decent_colors for an example"""

    if name == "gray" or name == "grayscale":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.84, 0.61, 0.34, 0.00]
        green = [1.00, 0.84, 0.61, 0.34, 0.00]
        blue  = [1.00, 0.84, 0.61, 0.34, 0.00]
    if name == "blue":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [0.88, 0.66, 0.44, 0.22, 0.00]
        green = [1.00, 0.75, 0.51, 0.27, 0.04]
        blue  = [1.00, 0.90, 0.81, 0.71, 0.62]
    # elif name == "whatever":
        # (define more palettes)
    else:
        # default palette, looks cool
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [0.00, 0.00, 0.87, 1.00, 0.51]
        green = [0.00, 0.81, 1.00, 0.20, 0.00]
        blue  = [0.51, 1.00, 0.12, 0.00, 0.00]

    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)

    npoints = len(s)
    TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    gStyle.SetNumberContours(ncontours)
