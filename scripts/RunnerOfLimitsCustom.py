#!/usr/bin/env python

import os,sys
import argparse
parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument("--dirname", dest="dirname", default="900", type=str)
opt = parser.parse_args()

tempfile= '''
{
    "LTDIR" : "/src/HiggsAnalysis/bbggLimits/DIRNAME_TYPE",
    "signal" : {
        "types" : ["HighMass", "LowMass"],
	"signalModelCard" : "/src/HiggsAnalysis/bbggLimits/LimitSetting/Models/models_2D_higgs.rs"
    },
    "other" : {
        "doDoubleSidedCB": 1,
        "Combinelxbatch" : 0,
        "version" : 66,
	"integratedLumi" : 36.5,
        "energy" : "13TeV",
        "higgsMass" : 125.0,
        "addHiggs" : 1,
        "doBlinding" : 0,
        "doBands" : 0,
        "ncat" : 2,
        "analysisType" : "fitTo2D_resSearch_withRegKinFit",
        "doSingleLimit" : 0,
        "drawSignalFit" : 0,
        "drawBackgroundFit" : 0,
        "useSigTheoryUnc" : 0,
	"doBrazilianFlag" : false,
	"runCombine" :true,
	"combineOption" : 2,
	"minMggMassFit" : 100,
	"maxMggMassFit" : 180,
	"minMjjMassFit"	: 60,
	"maxMjjMassFit"	: 180,
	"minSigFitMgg" :115,
	"maxSigFitMgg" :135,
	"minSigFitMjj" :60,
	"maxSigFitMjj" :180,
	"minHigMggFit" :115,
	"maxHigMggFit" :135,
	"minHigMjjFit":60,
	"maxHigMjjFit":180,
	"HH":false,
	"base":true,
	"low":false,
	"obs":false,
	"twotag":false,
        "doBias":false,
        "biasConfig": "BiasStudies/ForBias.json"

    },
    "data" : {
        "name" :"DoubleEG"
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

jsonname = "json_"+opt.dirname+".json"
outfile = open(jsonname, "w+")
towrite = tempfile.replace("DIRNAME", opt.dirname)
outfile.write(towrite)
outfile.close()
command = "pyLimits.py -f " + jsonname + " -o LIMS_"+opt.dirname+" --nodes SM --overwrite -j 1 --extraLabel "+opt.dirname+" -v 5"
print command
os.system(command)
