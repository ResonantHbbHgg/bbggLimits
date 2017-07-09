import time,os,argparse
from ROOT import *

gROOT.SetBatch()
gStyle.SetOptStat(0)

parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument("-f", "--folder", dest="f", type=str)
parser.add_argument("-o", "--outfolder", dest="out", type=str)
opt = parser.parse_args()

vkl= [1.0, 0.0, 7.5,  15.0,    5.0,  10.0,    1.0,     2.4,     7.5, 10.0,    15.0,   -15.0,    2.4,     -15.0]
vkt= [1.0, 1.0, 2.5,  1.5,     2.25, 1.5,     0.5,     1.25,    2.0, 2.25,    0.5,    2.0,      2.25,    1.25]
vc2= [0.0, 0.0, -0.5, -3.0,    3.0,  -1.0,    4.0,     2.0,     0.5, 2.0,     1.0,    6.0,      2.0,     6.0]
vcg= [0.0, 0.0, 0.0,  -0.0816, 0.0,  -0.0956, -1.0,    -0.2560, 0.0, -0.2130, -0.0743, -0.1680, -0.0616, -0.0467] 
vc2g=[0.0, 0.0, 0.0,   0.3010, 0.0,  0.12,    -0.3780, -0.1480, 0.0, -0.0893, -0.0668, -0.5180, -0.1200, -0.515]

plots = [
['mgg', 'mgg', 100, 115, 135],
['mjj', 'mjj', 100, 70, 190],
['mtot_hm', 'mtot', 100, 350, 1000],
['mtot_lm', 'mtot', 100, 250, 350],
['HHTagger', 'HHTagger', 100, 0.6, 1.0]
]

counter = 0
for ii in range(0, len(vkl)):
 kl = vkl[ii]
 kt = vkt[ii]
 cg = vcg[ii]
 c2 = vc2[ii]
 c2g = vc2g[ii]
 print kl, kt, cg, c2, c2g
 #LT_NR_Nodes_All_merged_kl_10p0_kt_1p5_cg_-0p0956_c2_-1p0_c2g_0p12.root
 pointStr = 'kl_'+str(kl).replace('.', 'p')+'_kt_'+str(kt).replace('.', 'p')+'_cg_'+str(cg).replace('.', 'p')+'_c2_'+str(c2).replace('.', 'p')+'_c2g_'+str(c2g).replace('.', 'p')
 noname = str(ii)
 if ii == 1: noname = 'box'
 if ii == 0: noname = 'SM'
 rwfile = TFile(opt.f+'/LT_NR_Nodes_All_merged_'+pointStr+'.root')
 nofile = TFile(opt.f+'/LT_output_GluGluToHHTo2B2G_node_'+noname+'_13TeV-madgraph.root')
 rwtree = rwfile.Get("TCVARS")
 notree = nofile.Get("TCVARS")
 for pl in plots:
   print pl
   rwth = TH1F('rwth'+pl[0], ';'+pl[1]+';Events', pl[2], pl[3], pl[4])
   noth = TH1F('noth'+pl[0], ';'+pl[1]+';Events', pl[2], pl[3], pl[4])
   rwth.SetTitle('_Node_'+str(ii))
   rwtree.Draw(pl[1]+'>>rwth'+pl[0], 'new_evWeight*(cut_based_ct>-1)')
   notree.Draw(pl[1]+'>>noth'+pl[0], 'evWeight*(cut_based_ct>-1)')
   rwth.SetLineColor(kRed)
   noth.SetLineColor(kBlue)
   thismax = rwth.GetMaximum()
   rwth.SetMaximum(thismax*1.2)
   leg = TLegend(.7, 0.7, 0.89, 0.89)
   leg.AddEntry(rwth, 'ARW', 'l')
   leg.AddEntry(noth, 'Node', 'l')
   c = TCanvas('c', 'c', 800, 600)
   rwth.Draw()
   noth.Draw('same')
   leg.Draw('same')
   c.SaveAs(opt.out+'/'+pl[0]+'_Node_'+str(ii)+'.pdf')
   c.SaveAs(opt.out+'/'+pl[0]+'_Node_'+str(ii)+'.png')
