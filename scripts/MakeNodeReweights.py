import time,os,argparse

parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument("-f", "--file", dest="f", type=str)
parser.add_argument("-o", "--outfolder", dest="out", type=str)
opt = parser.parse_args()

vkl= [1.0, 0.0, 7.5,  15.0,    5.0,  10.0,    1.0,     2.4,     7.5, 10.0,    15.0,   -15.0,    2.4,     -15.0]
vkt= [1.0, 1.0, 2.5,  1.5,     2.25, 1.5,     0.5,     1.25,    2.0, 2.25,    0.5,    2.0,      2.25,    1.25]
vc2= [0.0, 0.0, -0.5, -3.0,    3.0,  -1.0,    4.0,     2.0,     0.5, 2.0,     1.0,    6.0,      2.0,     6.0]
vcg= [0.0, 0.0, 0.0,  -0.0816, 0.0,  -0.0956, -1.0,    -0.2560, 0.0, -0.2130, -0.0743, -0.1680, -0.0616, -0.0467] 
vc2g=[0.0, 0.0, 0.0,   0.3010, 0.0,  0.12,    -0.3780, -0.1480, 0.0, -0.0893, -0.0668, -0.5180, -0.1200, -0.515]

counter = 0
for ii in range(0, len(vkl)):
 kl = vkl[ii]
 kt = vkt[ii]
 cg = vcg[ii]
 c2 = vc2[ii]
 c2g = vc2g[ii]
 print kl, kt, cg, c2, c2g
 command = 'python scripts/MakeARWTree.py -f ' + opt.f + ' -o ' + opt.out + ' --kl ' + str(kl) + ' --kt ' + str(kt) + ' --cg ' + str(cg) + ' --c2 ' + str(c2) + ' --c2g ' + str(c2g) + ' >/dev/null 2>&1 & '
 print command
 os.system(command)
 counter += 1
 if counter > 10:
  print 'sleeping...'
  time.sleep(60)
  print 'woke up!'
  counter = 0

