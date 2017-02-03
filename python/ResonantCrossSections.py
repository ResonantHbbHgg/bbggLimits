from ROOT import *
from array import array

masses = [250, 260, 300, 400, 500, 600, 700, 750, 800, 900]

#k/pl = 0.5
grav = [1e-8, 5.35e-02, 1.06e+00, 1.81e+00, 9.58e-01, 4.47e-01, 2.21e-01, 1.59e-01, 1.17e-01, 6.36e-02]
#k/pl = 1.0
grav2 = [4*x for x in grav]

gr_grav = TGraph(len(masses), array('d', masses), array('d', grav2))
grav_leg = "Bulk Graviton, #kappa/M_{Pl} = 1.0"
#LR = 3 TeV
rad = [6.7, 6.62e+00, 6.17e+00, 2.53e+00, 1.30e+00, 7.87e-01, 5.06e-01, 4.11e-01, 3.40e-01, 2.34e-01]
#LR = 1 TeV
#rad = [60, 5.96e+01, 5.56e+01, 2.28e+01, 1.17e+01, 7.08e+00, 4.55e+00, 3.70e+00, 3.06e+00, 2.10e+00]
rad_leg = "Bulk Radion, #Lambda_{R} = 3 TeV"

gr_rad = TGraph(len(masses), array('d', masses), array('d', rad))

gr_grav.SetLineColor(kRed)
gr_grav.SetLineStyle(2)
gr_grav.SetLineWidth(3)

gr_rad.SetLineColor(kRed)
gr_rad.SetLineStyle(2)
gr_rad.SetLineWidth(3)

