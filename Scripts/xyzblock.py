#!/usr/bin/env python

#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
#
#* File Name : xyzblock.py
#
#* Purpose :
#
#* Creation Date : 04-06-2013
#
#* Last Modified : Tue 04 Jun 2013 11:42:36 AM ART
#
#* Created By :  Ezequiel Castillo
#
#_._._._._._._._._._._._._._._._._._._._._.

FILENAME = 'Las256PuntasFinales.xyz'
FRAME = 2

OUTFILE = 'outfile.xyz'

data = open(FILENAME, 'r')

natom = int(data.readline())
lines_per_block = natom + 2

ignore_lines = (FRAME-1) * lines_per_block

data.seek(0)
for i in range(ignore_lines):
    data.next()

outf = open(OUTFILE, 'w')

i = 0
for line in data:
    outf.write(line)
    i += 1
    if i >= lines_per_block:
        break

data.close()
outf.close()
