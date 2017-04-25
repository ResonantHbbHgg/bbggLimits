#Sample name, sum of weights, xsec*br (in fb)
hggbr = 2.27e-3
genDiPhoFilterFactor = 1./(1 - 0.06)
#factor comes from the fact that pythia includes H->gamma*gamma->llgamma, while branching fraction from YR4  does not
#for more info, see: https://indico.cern.ch/event/598436/contributions/2529023/attachments/1434414/2205057/Zenz-News-Hgg-27Mar-v2.pdf
SMHiggsNodes = [
['output_GluGluHToGG_M-125_13TeV_powheg_pythia8.root', 998200.0, 48.5*1000.*hggbr*genDiPhoFilterFactor],
['output_VBFHToGG_M-125_13TeV_powheg_pythia8.root', 3738258.272705078, 3.748*1000.*hggbr*genDiPhoFilterFactor],
['output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root', 1862202.21875, 2.257*1000*hggbr*genDiPhoFilterFactor],
['output_ttHToGG_M125_13TeV_powheg_pythia8_v2.root', 431817.3866882324, 0.5071*1000*hggbr*genDiPhoFilterFactor],
['output_bbHToGG_M-125_4FS_yb2_13TeV_amcatnlo.root', 297462.02209472656, 0.534403*1000*hggbr*genDiPhoFilterFactor],
['output_bbHToGG_M-125_4FS_ybyt_13TeV_amcatnlo.root', -24582.49609375, 0.0464028*1000*hggbr*genDiPhoFilterFactor]
]
