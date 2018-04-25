import os

for x in range(0,20):

    cut=-0.1+x*0.02
    print cut
    if cut < 0:
        command = "pyLimits.py -f conf_default_local.json -o outDir_m"+str(abs(cut*10))+" --analyticalRW --ttHTaggerCut "+str(cut)
    else:
        command = "pyLimits.py -f conf_default_local.json -o outDir_"+str(abs(cut*10))+" --analyticalRW --ttHTaggerCut "+str(cut) 

    print command
    os.system(command)
