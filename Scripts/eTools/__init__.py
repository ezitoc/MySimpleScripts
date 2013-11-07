#!/usr/bin/env python
#-*- coding: utf-8 -*-

import sys
import os
import re
import Siesta
import Selection
from itertools import count
from collections import namedtuple
from pTools import *


atno_to_symbol = \
        {1: 'H',
         2: 'He',
         6: 'C',
         7: 'N',
         8: 'O',
        16: 'S',
        79: 'Au'}

symbol_to_atno = dict([[v,k] for k,v in atno_to_symbol.items()])

def get_dirs(a_dir=None):
    if not a_dir:
        a_dir = os.curdir
    files = os.listdir(a_dir)
    return [ f for f in files if os.path.isdir(os.path.join(a_dir, f)) ]

def sort_list_of_strings(a_list):
    return sorted(a_list, key=lambda x:[int(y) for y in x.split('.')])

Point = namedtuple('Point', 'x, y, z')

class Atom(object):

    """ Atom class. Must be instantiated as:
        an_atom = Atom(symbol, x, y, z)
        where symbol is either the atomic number or the atom name.
        x, y and z are the atomic coordinates."""

    last_id = count()

    def __init__(self, symbol, x, y, z, siesta_move=(1,1,1)):
        self.__position = Point(x, y, z)
        self.siesta_move = siesta_move

        if re.match(r'[a-zA-Z]+', symbol):
            self.symbol = symbol
        else:
            self.atno = symbol
        if not self.symbol:
            self.symbol = self.atno_to_symbol()

        self.index = self.last_id.next()

    def get_position(self):
        return self.__position

    def set_position(self, x, y, z):
        self.__position = Point(x, y, z)

    def translate(self, x, y, z):
        x0, y0, z0 = self.__position
        self.set_position(x0+x, y0+y, z0+z)

    def atno_to_symbol(self):
        return atno_to_symbol[self.atno]

    def symbol_to_atno(self):
        return symbol_to_atno[self.symbol]

    def __repr__(self):
        return '{name}({sym!r}, {pos.x}, {pos.y}, {pos.z})'.format(
            name=self.__class__.__name__, sym=self.symbol,
            pos=self.get_position())

    #def __repr__(self):
        #return '%s %10.4f %10.4f %10.4f' % (self.symbol, self.__position[0],
                                            #self.__position[1],
                                            #self.__position[2])


class Molecule(object):

    last_id = count()

    def __init__(self, name='Generic', atom_list=None):

        if name == 'Generic':
            self.name = 'Generic' + str(self.last_id.next())
        else:
            self.name = name

        if atom_list:
            self.atom_list = atom_list
        else:
            self.atom_list = []

    def add_atom(self, atom):
        self.atom_list.append(atom)

    def get_position(self):
        return [atom.get_position() for atom in self.atom_list]

    def load_from_file(self, filename, filetype='xyz'):


        if filetype == 'xyz':
            with open(filename, 'r') as f:
                lines = f.readlines()[2:]

        self.load_from_list(lines)

    def load_from_list(self, a_list):

        # Erase molecule
        self.atom_list = []

        for line in a_list:
            spline = line.split()
            if re.match(r'[a-zA-Z]+', spline[0]):
                an_atom = Atom(spline[0], float(spline[1]),
                            float(spline[2]),
                            float(spline[3]))
                self.add_atom(an_atom)

    def write_to_file(self, filename, type='xyz', mode='a'):

        with open(filename, mode) as f:
            if type == 'xyz':
                f.write('%i\n\n' % self.__len__())
                for atom in self.atom_list:
                    f.write('%s %10.4f %10.4f %10.4f\n' % (atom.symbol,
                                                         atom.get_position()[0],
                                                         atom.get_position()[1],
                                                         atom.get_position()[2]))



    def __getitem__(self, index):
        return self.atom_list[index]

    def __len__(self):
        return len(self.atom_list)

    #def __repr__(self):
        #str = 'This is a molecule named %s\n' % self.name
        #str = str + 'It has %d atoms\n' % len(self.atom_list)
        #return_str = ''
        #for atom in self.atom_list:
            #return_str = return_str + `atom` + '\n'
        #return return_str

class System(object):
    def __init__(self, molecule_list=None):
        self.molecule_list = molecule_list

    def add_molecule(self, molecule):
        self.molecule_list.append(molecule)
        self.number_of_atoms = self.__len__()

    def get_position(self):
        return [pos for mol in self.molecule_list for pos in mol.get_position()]

    def get_species(self):
        return list(set([atom.symbol for mol in self.molecule_list for atom in mol]))

    #FIXME
    #def xyz_to_zmatrix(self):

        #if os.path.isfile(self.InputFile):
            #print "File \""+self.InputFile+"\" found."
            #inputin = open(self.InputFile, 'r+')
        #else:
            #print "File \""+self.InputFile+"\" not found."
            #sys.exit(1)

        #inputin.write('\n#-- Atomic coordinates\n')
        #inputin.write('%block Zmatrix\n')
        #inputin.write('cartesian\n')

        #ostr = ''
        #FIXME: self.system
        #for i, atom in enumerate(self.system):
            #ostr = ostr + '%2s%16.8f%16.8f%16.8f%3d%3d%3d%8d\n' % \
                #(self.species_in_dict()[atom.symbol], atom.get_position()[0], atom.get_position()[1],
                 #atom.get_position()[2], atom.siesta_move[0],
                 #atom.siesta_move[1], atom.siesta_move[2], i+1)

        #return ostr
        #inputin.write("""%endblock Zmatrix""")
        #inputin.close()
        #file.close()

    def species_in_dict(self):
        return {specie: i+1 for i, specie in enumerate(self.get_species())}

    def __getitem__(self, index):
        return [atom for mol in self.molecule_list for atom in mol][index]

    def __len__(self):
        total_atoms = 0
        for molecule in self.molecule_list:
            total_atoms += len(molecule)
        return total_atoms

