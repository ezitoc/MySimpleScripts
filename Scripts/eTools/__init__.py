#!/usr/bin/env python
#-*- coding: utf-8 -*-

import sys
import os
import re
import types
from itertools import count
from collections import namedtuple
from pTools import *

#__all__ = ['Atom', 'Molecule', 'Selection']

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

        if re.match(r'[a-zA-Z]+.*', symbol):
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

    def remove_atom(self, i):
        print 'Deleting {atom}'.format(atom=self.atom_list[i])
        self.atom_list.pop(i)

    def load_from_file(self, filename, filetype='xyz'):
        with open(filename, 'r') as f:
            if filetype == 'xyz':
                lines = f.readlines()[2:]
            elif filetype == 'gro':
                lines = [l.split() for l in f.readlines()[2:-1]]
                lines = map(list, zip(*lines))
                for a_list in
                coords = map(float, [ lines[3], lines[4], lines[5] ])
                lines = [lines[1], lines[3], lines[4], lines[5]]
                lines = map(list, zip(*lines))
        self.load_from_list(lines, splitted=True)

    def load_from_list(self, a_list, splitted=False):
        # Erase molecule
        self.atom_list = []
        import pdb; pdb.set_trace()  # XXX BREAKPOINT
        for line in a_list:
            if splitted:
                spline = line
            else:
                spline = line.split()
            if re.match(r'[a-zA-Z]+.*', spline[0]):
                an_atom = Atom(spline[0], float(spline[1]),
                            float(spline[2]),
                            float(spline[3]))
                self.add_atom(an_atom)

    def write_to_file(self, filename='default.xyz', fmt='xyz'):
        with open(filename, 'w') as fout:
            fout.write('{tot_at:5d}\n\n'.format(tot_at=len(self.atom_list)))
            for atom in self.atom_list:
                fout.write('{sym:3s}{x:13.8f}{y:13.8f}{z:13.8f}\n'.format(sym=atom.symbol,
                                                                        x=atom.get_position().x,
                                                                        y=atom.get_position().y,
                                                                        z=atom.get_position().z))

    def move(self, x, y, z):
        #for atom in self.atom_list:
            #atom.translate(x, y, z)
        self + [x, y, z]
        return self

    def __add__(self, an_object):

        if isinstance(an_object, (self.__class__, Selection.Selection)):
            #self.atom_list += an_object.atom_list
            for atom in an_object:
                self.add_atom(atom)

        elif isinstance(an_object, Atom):
            self.add_atom(an_object)

        elif isinstance(an_object, types.ListType):
            if len(an_object) == 3:
                x, y, z = an_object
                for atom in self.atom_list:
                    atom.translate(x, y, z)

        else:
            raise IOError('Not addition supported of {name} object with \
                            "{obj}"'.format(name=self.__class__.__name__,
                            obj=an_object))
        return self

    def __sub__(self, an_object):

        if isinstance(an_object, (self.__class__, Selection.Selection)):
            #self.atom_list += an_object.atom_list
            for i, orig_at in enumerate(self.atom_list):
                for remove_at in an_object:
                    if orig_at.index == remove_at.index:
                        self.remove_atom(i)

        elif isinstance(an_object, Atom):
            for i, orig_at in enumerate(self.atom_list):
                if orig_at.index == an_object.index:
                    self.remove_atom(i)

        elif isinstance(an_object, types.ListType):
            if len(an_object) == 3:
                x, y, z = an_object
                for atom in self.atom_list:
                    atom.translate(-1*x, -1*y, -1*z)

        else:
            raise IOError('Not addition supported of {name} object with \
                            "{obj}"'.format(name=self.__class__.__name__,
                            obj=an_object))
        return self

    def __radd__(self, an_object):
        self.__add__(an_object)

    def __rsub__(self, an_object):
        self.__sub__(an_object)

    def __getitem__(self, index):
        return self.atom_list[index]

    def __len__(self):
        return len(self.atom_list)

    def __repr__(self):
        return '{name}({atoms})'.format(name=self.__class__.__name__,
                atoms=self.atom_list)


class Selection(Molecule):

    def __init__(self, a_system):

        atom_list = []

        if isinstance(a_system, Atom):
            atom_list = a_system

        elif isinstance(a_system, Molecule):
            for atom in a_system:
                atom_list.append(atom)

        elif isinstance(a_system, System):
            for atom in a_system:
                atom_list.append(atom)

        self.atom_list = atom_list
        self.orig_alist = atom_list

    def _filter(self, func):
        return filter(func, self.atom_list)

    def remove_from_list(self, an_element):
        if an_element in self.atom_list:
            self.atom_list.remove(an_element)

    def by_symbol(self, symbol):
        self.atom_list = self._filter(lambda a: a.symbol == symbol)
        return self

    def by_xrange(self, xmin, xmax):
        self.atom_list = self._filter(lambda a: xmin <= a.get_position().x <=
                                      xmax)
        return self

    def by_yrange(self, ymin, ymax):
        self.atom_list = self._filter(lambda a: ymin <= a.get_position().y <=
                                      ymax)
        return self

    def by_zrange(self, zmin, zmax):
        self.atom_list = self._filter(lambda a: zmin <= a.get_position().z <=
                                      zmax)
        return self

    def sphere(self, radius, center_atom=None, pos=None):
        if not center_atom and not pos:
            pos = Point(0, 0, 0)
        elif pos:
            pos = Point(pos[0], pos[1], pos[2])
        elif center_atom:
            pos = center_atom.get_position()
        self.atom_list = self._filter(lambda a: (a.get_position().x-pos.x)**2 +
                                      (a.get_position().y-pos.y)**2 +
                                      (a.get_position().z-pos.z)**2 <=
                                      radius**2)
        return self

    def update_sel(self):
        for m_atom in self.atom_list:
            for i, o_atom in enumerate(self.orig_alist):
                if m_atom.index == o_atom.index:
                    self.orig_alist[i] = m_atom
        self.atom_list = self.orig_alist
        return self.atom_list
