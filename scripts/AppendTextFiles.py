import argparse, os, json

parser =  argparse.ArgumentParser(description='Limit Tree maker')
parser.add_argument("-f1", dest="f1", type=str)
parser.add_argument("-f2", dest="f2", type=str)
parser.add_argument('-o', '--outFile', dest='outf', type=str)
opt = parser.parse_args()

file1 = open(opt.f1, "r")
file2 = open(opt.f2, "r")
fout = open(opt.outf, "w+")
fouttranls = open(opt.outf.replace(".txt", "_transl.txt"), "w+")
fjson = open(opt.outf.replace('.txt', '.json'), 'w+')


scan_2d = {
'kl': [float(i)*2.0 for i in range(-10,11)], #-20, -17.5, -15, -12.5, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20],
'kt': [float(i)/4.0 for i in range(-10,11)],#-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5],
'cg': [0.0],
'c2': [0.0],
'c2g': [0.0]
}

MyInputs = { }

for kl in scan_2d['kl']:
 MyInputs[str(kl)] = {}
 for kt in scan_2d['kt']:
  MyInputs[str(kl)][str(kt)] = ""

for line in file1:
  line_split = line.split()
  if len(line_split) < 9: continue
  kl = line_split[1]
  kt = line_split[2]
  l1 = line_split[3]
  l2 = line_split[4]
  l3 = line_split[5]
  l4 = line_split[6]
  l5 = line_split[7]
  l6 = line_split[8]
  MyInputs[str(kl)][str(kt)] = l1 + ' ' + l2 + ' ' + l3 + ' ' + l4 + ' ' + l5 + ' ' + l6

file1.close()

for line in file2:
  line_split = line.split()
  if len(line_split) < 9: continue
  kl = line_split[1]
  kt = line_split[2]
  if MyInputs[str(kl)][str(kt)] is not "": continue
  l1 = line_split[3]
  l2 = line_split[4]
  l3 = line_split[5]
  l4 = line_split[6]
  l5 = line_split[7]
  l6 = line_split[8]
  MyInputs[str(kl)][str(kt)] = l1 + ' ' + l2 + ' ' + l3 + ' ' + l4 + ' ' + l5 + ' ' + l6

file2.close()

jsonlimits = []

counter = 0
for kl in scan_2d['kl']:
 for kt in scan_2d['kt']:
  towrite =  str(counter) + ' '
  towrite += str(kl) + ' '
  towrite += str(kt) + ' '
  towrite += MyInputs[str(kl)][str(kt)]
  towrite += ' \n'
  fout.write(towrite)
  fouttranls.write(str(counter) + ' ' + str(kl) + ' ' + str(kt) + '\n')

  lim_split = MyInputs[str(kl)][str(kt)].split()
#  print lim_split, kl, kt
  if len(lim_split) < 6:
    lim_split = [0,0,0,0,0,0]
    print kl, kt

  someLimits = { 'one_sigma': [float(lim_split[3]), float(lim_split[2])] , 'two_sigma': [float(lim_split[5]), float(lim_split[4])], 'observed': float(lim_split[1]), 'expected': float(lim_split[0]) }

  tojson = { 'parameters': (kl,kt), 'limits': someLimits}

  jsonlimits.append(tojson)

  counter += 1

json.dump(jsonlimits, fjson)

fjson.close()
fout.close()
fouttranls.close()

