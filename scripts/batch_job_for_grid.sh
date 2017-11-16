#!/bin/bash

workDir=$1
outDir=$2
Point=$3

echo $workDir $outDir $Point
cd $workDir

ls 

eval `scramv1 runtime -sh`

pyLimits.py -f conf_NRW_1.json -o $outDir --points $Point -j1 -v2 --overwrite
