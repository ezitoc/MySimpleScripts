#!/usr/bin/env python

#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
#
#* File Name : cutmol.py
#
#* Purpose :
#
#* Creation Date : 04-06-2013
#
#* Last Modified : Tue 04 Jun 2013 04:20:03 PM ART
#
#* Created By :  Ezequiel Castillo
#
#_._._._._._._._._._._._._._._._._._._._._.

import numpy as np

FILENAME = 'outfile.xyz'
MIN = 6.0
MAX = 32.0

OUTFILE = 'newout.xyz'

infile = open(FILENAME, 'r')
data = infile.readlines()
infile.close()

head = data[:2]
pos = data[2:]

outf = open(OUTFILE, 'w')

new_atoms = []
for line in pos:
    value = float(line.split()[3])
    if value > MIN and value < MAX:
        new_atoms.append(line)

new_atoms.insert(0, str(len(new_atoms))+'\n\n')

for line in new_atoms:
    outf.write(line)

outf.close()
