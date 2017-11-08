#!/usr/bin/env python

import os,sys
import argparse
parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument("--dirname", dest="dirname", default="900", type=str)
parser.add_argument("--dirloc", dest="dirloc", default="/src/HiggsAnalysis/bbggLimits/", type=str)
parser.add_argument("--mass", dest="mass", default="-10", type=str)
parser.add_argument("--resType", dest="restype", default="Radion", type=str)
parser.add_argument('--extra', dest='extra', default='', type=str)
parser.add_argument('--jsonName', dest='jsonname', default='', type=str)
opt = parser.parse_args()


tempfile= '''
{
    "LTDIR" : "DIRLOC/DIRNAME_TYPE",
    "signal" : {
	"types" : SIGTYPES,
	"signalModelCard" : "/src/HiggsAnalysis/bbggLimits/Models/models_2D_higgs_mjj70.rs"
    },
    "other" : {
        "doDoubleSidedCB": 1,
        "Combinelxbatch" : 1,
	"integratedLumi" : 36.5,
        "energy" : "13TeV",
        "higgsMass" : 125.0,
        "addHiggs" : AHIGGS,
        "doBlinding" : 0,
        "doBands" : 0,
        "ncat" : 2,
        "doSingleLimit" : 0,
        "drawSignalFit" : 0,
        "drawBackgroundFit" : 0,
        "useSigTheoryUnc" : 0,
	"doBrazilianFlag" : false,
	"runCombine" :true,
	"combineOption" : 1,
	"minMggMassFit" : 100,
	"maxMggMassFit" : 180,
	"minMjjMassFit"	: 70,
	"maxMjjMassFit"	: 190,
	"minSigFitMgg" :115,
	"maxSigFitMgg" :135,
	"minSigFitMjj" :70,
	"maxSigFitMjj" :190,
	"minHigMggFit" :115,
	"maxHigMggFit" :135,
	"minHigMjjFit":70,
	"maxHigMjjFit":190,
	"HH":false,
	"base":true,
	"low":false,
	"obs":false,
	"twotag":false,
        "doBias":false,
        "biasConfig": "BiasStudies/ForBias.json"

    },
    "data" : {
        "name" :"DoubleEGRMASS"
    },
    "higgs" : {
	"type" : {
		"vbf": "VBFHToGG_M-125_13TeV_powheg_pythia8",
		"bbh": "bbHToGG_M-125_13TeV_amcatnlo",
		"ggh": "GluGluHToGG_M-125_13TeV_powheg_pythia8",
		"tth": "ttHToGG_M125_13TeV_powheg_pythia8_v2",
		"vh" : "VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8"
		}
    }
}
'''
if opt.jsonname == '':
  jsonname = "json_"+opt.dirname+".json"
else:
  jsonname = opt.jsonname

sigTypes = '["HighMass", "LowMass"]'
toRun = ''
#toRun = '--nodes SM'
addHiggs = '1'
rMass = ''
if '-10' not in opt.mass:
  sigTypes = '["'+opt.restype+'"]'
  toRun = '--mass ' + opt.mass
  addHiggs = '0'
  rMass = '_M-'+opt.mass

outfile = open(jsonname, "w+")
towrite = tempfile.replace("DIRNAME", opt.dirname).replace("DIRLOC", opt.dirloc).replace('SIGTYPES', sigTypes).replace('AHIGGS', addHiggs).replace('RMASS', rMass)
outfile.write(towrite)
outfile.close()

if 'extraLabel' in opt.extra:
  command = "pyLimits.py -f " + jsonname + " -o LIMS_"+opt.dirname+" "+toRun+" --overwrite -j 1  -v 5 " + opt.extra
else:
  command = "pyLimits.py -f " + jsonname + " -o LIMS_"+opt.dirname+" "+toRun+" --overwrite -j 1 --extraLabel "+opt.dirname+" -v 5 " + opt.extra
print command
os.system(command)
