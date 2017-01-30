#!/usr/bin/env python

from ROOT import *
import sys, getopt, os
import argparse
import glob

parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument('-f', '--inputFolders', dest="folders", default=None, type=str, nargs='+', required=True,
                    help="input folders")
parser.add_argument('-n', '--inputNames', dest="names", default=None, type=str, nargs='+', required=True,
                    help="input folders cat names")
parser.add_argument('-l', '--lumi', dest='lumi', default='36.5', type=str, help='Integrated luminosoty')
parser.add_argument('--label', dest='label', default='', type=str, help='Label')

opt = parser.parse_args()

'''
-rw-r--r--. 1 rateixei zh 6548 Jan 30 02:26 higgsCombine_Node_SM_CatBased400CTS_qt_central.HybridNew.mH125.quant0.500.root
-rw-r--r--. 1 rateixei zh 6529 Jan 30 03:33 higgsCombine_Node_SM_CatBased400CTS_qt_m1s.HybridNew.mH125.quant0.160.root
-rw-r--r--. 1 rateixei zh 6532 Jan 30 01:36 higgsCombine_Node_SM_CatBased400CTS_qt_p1s.HybridNew.mH125.quant0.840.root
-rw-r--r--. 1 rateixei zh 6537 Jan 30 01:45 higgsCombine_Node_SM_CatBased400CTS_qt_p2s.HybridNew.mH125.quant0.975.root
'''

def main(argv):
  quantiles = ['0.025', '0.160', '0.500', '0.840', '0.975']
  qt_names = ['m2s', 'm1s', 'central', 'p1s', 'p2s']

  for ff in opt.folders:
    for iqt, qt in enumerate(quantiles):
      print glob.glob(ff+"/*.root")

if __name__ == "__main__":
  main(sys.argv[1:])

