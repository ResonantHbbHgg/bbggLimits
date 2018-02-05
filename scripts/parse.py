#!/usr/bin/env python

import numpy as np
np.set_printoptions(precision=3)

#tag = "_Node_SM"
tag = "_ARW_"

for m in ['HighMass', 'LowMass']:
    print
    print 5*"*"
    print m
    print 3*"*"
    print
    
    fileName = 'LIMS_NewBTag_central/'+m+tag+'/datacards/hhbbgg_13TeV_DataCard.txt'
    names = None
    central= None
    with open(fileName) as fp:
        line = fp.readline()
        cnt = 1
        while line:
            #if cnt in [39, 40, 41]:
            #    print cnt
            #    print line
            if cnt==39:
                names = line.split()
            if cnt==41:
                central = np.array(line.split()[1:], dtype=np.float)
                
            line = fp.readline()
            cnt += 1
    
    # print names
    # print central
    # bsys = ["lf"]
    bsys = ["jes", "lf", "hf", "hfstats1", "hfstats2", "lfstats1", "lfstats2", "cferr1", "cferr2"]

    for sy in bsys:
        yields = {}
        for var in ['up', 'down']:
            s = var+'_'+sy
            fileName = 'LIMS_NewBTag_'+s+'/'+m+tag+'/datacards/hhbbgg_13TeV_DataCard.txt'
            with open(fileName) as fp:
                line = fp.readline()
                cnt = 1
                while line:
                    if cnt==41:
                        yields[var] = np.array(line.split()[1:], dtype=np.float)
                    #if cnt in [39,41]:
                    #    print line
                
                    line = fp.readline()
                    cnt += 1
    
        
        #print '\n ** Difference for', m, sy, 'DOWN'
        #print central-yields['down']
        dn = yields['down']/central

        #print '\n ** Difference for', m, sy, 'UP'
        #print central-yields['up']
        up = yields['up']/central

        sys_max = np.maximum(abs(1-dn),abs(1-up))
        sys_avg = (abs(1-dn) + abs(1-up))/2
        # print sys_max
        outSys = 'CMS_btag_'+sy+' lnN \t'
        for i in xrange(len(central)):
            if i in [1, 8]:
                outSys += ' \t-'
            else:
                outSys += ' \t%.3f/%.3f' % (dn[i], up[i])
                #if i < 7:
                #    outSys += ' \t%.3f' % (1+sys_avg[i])
                #else:
                #    outSys += ' \t%.3f' % (1-sys_avg[i])

        #print "\n ** Final for this syst:", sy, m
        print outSys
            
