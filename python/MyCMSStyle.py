from ROOT import *

def SetAxisTextSizes(obj, yoff=0, ysize=1, xoff=0, xsize=1):
  obj.GetYaxis().SetTitleOffset(1.1+yoff)
  obj.GetYaxis().SetTitleSize(0.0425*ysize)
  obj.GetYaxis().SetLabelSize(0.04*ysize)
  obj.GetXaxis().SetTitleOffset(1.1+xoff)
  obj.GetXaxis().SetTitleSize(0.0425*xsize)
  obj.GetXaxis().SetLabelSize(0.04*xsize)
  try:
    obj.GetZaxis().SetTitleOffset(1.1)
    obj.GetZaxis().SetTitleSize(0.0425)
    obj.GetZaxis().SetLabelSize(0.04)
  except AttributeError:
    a=1

def SetGeneralStyle():
  gStyle.SetFrameLineWidth(2)

def SetPadStyle(obj):
  obj.SetTicky()
  obj.SetTickx()

def DrawCMSLabels(obj, lumi=''):
  pad = obj.cd()
  l = pad.GetLeftMargin()
  t = pad.GetTopMargin()
  r = pad.GetRightMargin()
  b = pad.GetBottomMargin()
  lat = TLatex()
  lat.SetTextSize(0.045)
  lat.SetTextAlign(11)
  lat.SetTextFont(42)
  cmsTag = "#bf{CMS}"
  lumiTag = lumi+' fb^{-1} (13 TeV)'
  if lumi == '':
    cmsTag = "#bf{CMS} #it{Simulation}"
    lumiTag = '(13 TeV)'
  lat.DrawLatexNDC(l+0.01, 1-t+0.02, cmsTag)
  lat.SetTextAlign(31)
  lat.DrawLatexNDC(1-r-0.001, 1-t+0.02, lumiTag)

