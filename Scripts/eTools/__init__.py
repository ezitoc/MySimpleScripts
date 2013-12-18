#!/usr/bin/env python
#-*- coding: utf-8 -*-

import sys
import os
import re
import types
from itertools import count
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

DefaultSiestaKeys = \
    {'SystemName': 'siesta',
    'SystemLabel': 'siesta',
    'NumberOfSpecies': None,
    'NumberOfAtoms': None,
    'PAO_BasisType': 'split',
    'PAO_BasisSize': 'DZP',
    'PAO_EnergyShift': '0.01 eV',
    'PAO_SplitNorm': '0.15',
    'MeshCutoff': '200.0 Ry',
    'AtomicCoordinatesFormat': 'Ang',
    'AtomCoorFormatOut': 'Ang',
    'WriteCoorXmol': 'T',
    'WriteMDXmol': 'T',
    'WriteMDhistory': 'T',
    'MD_UseSaveXV': 'T',
    'DM_UseSaveDM': 'T',
    'MD_UseSaveCG': 'T',
    'XC_functional': 'GGA',
    'XC_authors': 'PBE',
    'SpinPolarized': 'T',
    'MaxSCFIterations': '500',
    'DM_MixingWeight': '0.01',
    'DM_NumberPulay': '4',
    'DM_Tolerance': '0.0001',
    'SolutionMethod': 'diagon',
    'MD_TypeOfRun': 'CG',
    'MD_MaxForceTol': '0.01 eV/Ang',
    'MD_NumCGsteps': '500',
    'MD_MaxCGDispl': '0.05  Ang',
    'LatticeConstant': '1.0 Ang',
    'ZM_UnitsLength': 'Ang',
    'ZM_UnitsAngle': 'deg'}

DefaultSiestaBlockKeys = \
    {'ChemicalSpeciesLabel': None,
    'Kgrid_Monkhorst_Pack': None,
    'LatticeVectors': None,
    'AtomicCoordinatesOrigin': [0, 0, 0],
    'Zmatrix': None}

SiestaKeysFormat = \
    {'ChemicalSpeciesLabel': '%3i%5i%10s\n',
    'Kgrid_Monkhorst_Pack': '%4i%4i%4i%8.1f\n',
    'LatticeVectors': None,
    'AtomicCoordinatesOrigin': None,
    'Zmatrix': '%4i%12.4f%12.4f%12.4f%3i%3i%3i\n'}

def get_dirs(a_dir=None):
    if not a_dir:
        a_dir = os.curdir
    files = os.listdir(a_dir)
    return [ f for f in files if os.path.isdir(os.path.join(a_dir, f)) ]

def sort_list_of_strings(a_list):
    return sorted(a_list, key=lambda x:[int(y) for y in x.split('.')])

def is_mono_dim_list(a_list):
    for elem in a_list:
        if isinstance(elem, types.ListType):
            return False
    return True

class AttrDict(dict):
    """ Class for accessing dict keys like an attribute  """
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


class Atom(object):

    """ Atom class. Must be instantiated as:
        an_atom = Atom(symbol, x, y, z)
        where symbol is either the atomic number or the atom name.
        x, y and z are the atomic coordinates."""

    last_id = count()

    def __init__(self, symbol, pos, siesta_move=(1,1,1)):

        #self.pos = Point(x, y, z)
        self.set_position(pos)
        self.siesta_move = siesta_move

        if re.match(r'[a-zA-Z]+.*', symbol):
            self.symbol = symbol
        else:
            self.atno = symbol
        if not self.symbol:
            self.symbol = self.atno_to_symbol()

        self.index = self.last_id.next()

    def set_position(self, pos):
        #if pos:
        if len(pos) == 3:
            self.pos = np.array(pos)

    def get_position(self, coord=None):
        if coord:
            if coord == 'x':
                return self.pos[0]
            if coord == 'y':
                return self.pos[1]
            if coord == 'z':
                return self.pos[2]
        else:
            return self.pos

    def translate(self, vec):
        self.set_position(self.pos + vec)

    def atno_to_symbol(self):
        return atno_to_symbol[self.atno]

    def symbol_to_atno(self):
        return symbol_to_atno[self.symbol]

    def __add__(self, vec):
        self.pos += vec

    def __sub__(self, vec):
        self.pos -= vec

    def __radd__(self, an_object):
        self.__add__(an_object)

    def __rsub____(self, an_object):
        self.__sub__(an_object)

    def __repr__(self):
        return '{name}({sym!r}, {pos})'.format(
            name=self.__class__.__name__, sym=self.symbol,
            pos=self.pos)


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

        if filetype == 'xyz':
            firstAtLine = 2
            lastAtLine = None
            nameCol = 0
            firstPosCol = 1
            lastPosCol = 4
            toAng = 1
        elif filetype == 'gro':
            firstAtLine = 2
            lastAtLine = -1
            nameCol = 1
            firstPosCol = 3
            lastPosCol = 6
            toAng = 10

        with open(filename, 'r') as f:

            lines = [l.split() for l in f.readlines()[firstAtLine:lastAtLine]]
            lines = map(list, zip(*lines))
            namelst = lines[nameCol]
            lines = [map(float, a_list) for a_list in lines[firstPosCol:lastPosCol]]
            posarr = np.array(lines)*toAng
            posarr = posarr.T

        self.atom_list = [Atom(namelst[i], pos) for i, pos in enumerate(posarr)]

    def write_to_file(self, filename='default.xyz', fmt='xyz'):
        with open(filename, 'w') as fout:
            fout.write('{tot_at:5d}\n\n'.format(tot_at=len(self.atom_list)))
            for atom in self.atom_list:
                fout.write('{sym:3s}{x:13.8f}{y:13.8f}{z:13.8f}\n'.format(sym=atom.symbol,
                                                                        x=atom.pos[0],
                                                                        y=atom.pos[1],
                                                                        z=atom.pos[2]))

    def get_geometry_center(self):
        asum = np.zeros(3)
        for atom in self.atom_list:
            asum += atom.pos
        asum /= len(self.atom_list)
        return asum

    def get_species(self):
        return list(set([atom.symbol for atom in self.atom_list]))
        #return {i+1:specie for i, specie in enumerate(species)}
        #return OrderedDict({i+1:specie for i, specie in enumerate(species)})

    def __add__(self, an_object):

        if isinstance(an_object, (self.__class__, Selection)):
            #self.atom_list += an_object.atom_list
            for atom in an_object:
                self.add_atom(atom)

        elif isinstance(an_object, Atom):
            self.add_atom(an_object)

        elif isinstance(an_object, (types.ListType, np.array)):
            for atom in self.atom_list:
                atom + an_object

        else:
            raise IOError('Not addition supported of {name} object with \
                            "{obj}"'.format(name=self.__class__.__name__,
                            obj=an_object))
        return self

    def __sub__(self, an_object):

        if isinstance(an_object, (self.__class__, Selection)):
            #self.atom_list += an_object.atom_list
            for i, orig_at in enumerate(self.atom_list):
                for remove_at in an_object:
                    if orig_at.index == remove_at.index:
                        self.remove_atom(i)

        elif isinstance(an_object, Atom):
            for i, orig_at in enumerate(self.atom_list):
                if orig_at.index == an_object.index:
                    self.remove_atom(i)

        elif isinstance(an_object, (types.ListType, np.array)):
            for atom in self.atom_list:
                atom - an_object

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
        #return self.atom_list[index-1]

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
        self.atom_list = self._filter(lambda a: xmin <= a.pos[0] <=
                                      xmax)
        return self

    def by_yrange(self, ymin, ymax):
        self.atom_list = self._filter(lambda a: ymin <= a.pos[1] <=
                                      ymax)
        return self

    def by_zrange(self, zmin, zmax):
        self.atom_list = self._filter(lambda a: zmin <= a.pos[2] <=
                                      zmax)
        return self

    def sphere(self, radius, pos=None, center_atom=None):
        if pos.any():
            center = pos
        elif center_atom:
            center = center_atom.pos
        else:
            center = np.zeros(3)
        self.atom_list = self._filter(lambda a: np.linalg.norm(a.pos - center) <=
                                      radius**2)
        return self

    def update_sel(self):
        for m_atom in self.atom_list:
            for i, o_atom in enumerate(self.orig_alist):
                if m_atom.index == o_atom.index:
                    self.orig_alist[i] = m_atom
        self.atom_list = self.orig_alist
        return self.atom_list


class Siesta(object):

    def __init__(self, sys=None, inFile='input.fdf', outFile='output.out',
            errFile='err.log', siesta='siesta', suffix='-PBE'):
        self.system = sys
        self.inFile = inFile
        self.outFile = outFile
        self.errFile = errFile
        self.siesta = siesta
        self.pseudoSuffix = suffix
        self.keys = AttrDict(DefaultSiestaKeys)
        self.blockKeys = AttrDict(DefaultSiestaBlockKeys)

        self.species = AttrDict({specie: i+1 for i, specie in
                                 enumerate(self.system.get_species())})
        self.siesta_system_attr()
        self.setSysDepKeys()
        #self.xyzFile = SystemLabel+'.xyz'
        #self.aniFile = SystemLabel+'.ANI'

    def siesta_system_attr(self):
        for atom in self.system:
            atom.zmatMove = [1, 1, 1]
            atom.specieNo = self.species[atom.symbol]

    def constrain(self, atomNumber, vec=[0, 0, 0]):
        self.system[atomNumber-1].zmatMove = vec

    def setSysDepKeys(self):
        self.keys.NumberOfSpecies = len(self.system.get_species())
        self.keys.NumberOfAtoms = len(self.system)
        self.blockKeys.ChemicalSpeciesLabel = \
        [[self.species[specie], symbol_to_atno[specie],
          specie+self.pseudoSuffix] for specie in self.species]

        self.blockKeys.Zmatrix = \
        [[atom.specieNo, atom.pos[0], atom.pos[1], atom.pos[2],
          atom.zmatMove[0], atom.zmatMove[1], atom.zmatMove[2]] for atom in
         self.system]

        self.setKgrid()
        self.setLatticeVectors()

    def setKgrid(self, kpoints=1):
        self.blockKeys.Kgrid_Monkhorst_Pack = np.hstack((np.identity(3),np.zeros((3,1))))

    def setLatticeVectors(self, latvec=None):
        if latvec:
            self.blockKeys.LatticeVectors = np.array(latvec)
        else:
            self.blockKeys.LatticeVectors = np.identity(3)*10

    #def zmaMolecule(self, mol):

    def createInput(self):
        inF = open(self.inFile, 'w')

        for a_key in self.keys:
            inF.write('{keyword:<30}{value:>15}\n'.format(keyword=a_key,
                                                          value=self.keys[a_key]))

        for a_key in self.blockKeys:
            inF.write('%block {}\n'.format(a_key))

            if a_key is 'ChemicalSpeciesLabel':
                for line in self.blockKeys[a_key]:
                    inF.write(SiestaKeysFormat['ChemicalSpeciesLabel'] % tuple(line))

            elif a_key is 'Zmatrix':
                for line in self.blockKeys[a_key]:
                    inF.write(SiestaKeysFormat['Zmatrix'] % tuple(line))

            elif a_key is 'Kgrid_Monkhorst_Pack':
                for line in self.blockKeys[a_key]:
                    inF.write(SiestaKeysFormat['Kgrid_Monkhorst_Pack'] % tuple(line))

            elif isinstance(self.blockKeys[a_key], np.ndarray):
                np.savetxt(inF, self.blockKeys[a_key], fmt='%6.3f')

            elif isinstance(self.blockKeys[a_key], types.ListType):
                if is_mono_dim_list(self.blockKeys[a_key]):
                    for elem in self.blockKeys[a_key]:
                        inF.write('%6.3f' % elem)
                    inF.write('\n')

            inF.write('%endblock {}\n'.format(a_key))
        inF.close()

    def run(self):
        inF = open(self.inFile, 'r')
        outF = open(self.outFile, 'w')
        errF = open(self.errFile, 'w')

        p = subprocess.Popen([self.siesta], stdin=inF, stdout=outF,
                             stderr=errF)
        p.wait()
        inF.close()
        outF.flush()
        errF.flush()
        outF.close()
        errF.close()

    def Input_load(self):
        if os.path.isfile(self.inFile):
            print 'File "{}" found.'.format(self.inFile)
            inF = open(inFile, 'r+')
        else:
            raise IOError('File "{}" not found.'.format(self.inFile))

        #for line in inF:


