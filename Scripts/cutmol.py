#!/usr/bin/env python

#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
#
#* File Name : cutmol.py
#
#* Purpose :
#
#* Creation Date : 04-06-2013
#
#* Last Modified : Wed 05 Jun 2013 07:07:26 PM ART
#
#* Created By :  Ezequiel Castillo
#
#_._._._._._._._._._._._._._._._._._._._._.

def cut(atom_list):
    atom_list.insert(0, str(len(atom_list))+'\n\n')

    return  atom_list

def inmovilizar(steady_list, else_list):
    print 'Los atomos inmovilizados son los primeros %i atomos' % (len(steady_list))
    final = steady_list + else_list

    final.insert(0, str(len(final))+'\n\n')

    return final





FILENAME = 'newout.xyz'
MIN = 6.0
MAX = 32.0

OUTFILE = 'inmovilizar.xyz'

infile = open(FILENAME, 'r')
data = infile.readlines()
infile.close()

head = data[:2]
pos = data[2:]

outf = open(OUTFILE, 'w')

new_atoms = []
garbage = []
for line in pos:
    value = float(line.split()[3])
    if value > MIN and value < MAX:
        new_atoms.append(line)
    else:
        garbage.append(line)

# Declare what you want to do

final_list = inmovilizar(garbage, new_atoms)

# Write to file
for line in final_list:
    outf.write(line)

outf.close()
