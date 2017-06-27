#!/usr/bin/env python

from ROOT import *
import os,sys,json,time,re
import logging
from shutil import copy
from pprint import pformat
# import pebble as pb
from multiprocessing import Pool, TimeoutError, current_process
from HiggsAnalysis.bbggLimits.LimitsUtil import *

gROOT.SetBatch()

parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument('-f', '--inputFile', dest="fname", type=str, default=None, required=True,
                    help="Json config file")
parser.add_argument('-o', '--outDir', dest="outDir", type=str, default=None,
                    help="Output directory (will be created).")
parser.add_argument('--nodes', dest="nodes", default=None, type=str, nargs='+',
                    choices=['2','3','4','5','6','7','8','9','10','11','12','13','SM','box','all'],
                    help = "Choose the nodes to run")
parser.add_argument('--mass', dest="mass", default=None, type=str, nargs='+',
                    choices=['250','260','270','280','300','320','340','350','400','450','500','550','600','650','700', '750', '800', '900', 'all'],
                    help = "Choose the resonant mass to run")
parser.add_argument('--points', dest="points", default=None, type=parseNumList, nargs='+',
                    help = "Choose the points in the grid to run")
parser.add_argument('--overwrite', dest="overwrite", action="store_true", default=False,
                    help="Overwrite the results into the same directory")
parser.add_argument("-v", dest="verb", type=int, default=0,
                    help="Verbosity level: 0 - Minimal or no messages; 1 - INFO; 2 - DEBUG; 3 - Go crazy")
parser.add_argument('-j', '--ncpu',dest="ncpu", type=int, default=2,
                    help="Number of cores to run on.")
parser.add_argument('-t', '--timeout',dest="timeout", type=int, default=None,
                    help="Per job timeout (in seconds) for multiprocessing. Jobs will be killed if run longer than this.")
parser.add_argument('--extraLabel', dest='extraLabel', default='',
                    help='Extra label')
parser.add_argument('--analyticalRW', dest='analyticalRW', action='store_true', default=False)
parser.add_argument('--kl', dest='ARW_kl', type=float, default=1.0)
parser.add_argument('--kt', dest='ARW_kt', type=float, default=1.0)
parser.add_argument('--cg', dest='ARW_cg', type=float, default=0.0)
parser.add_argument('--c2', dest='ARW_c2', type=float, default=0.0)
parser.add_argument('--c2g', dest='ARW_c2g', type=float, default=0.0)

opt = parser.parse_args()


