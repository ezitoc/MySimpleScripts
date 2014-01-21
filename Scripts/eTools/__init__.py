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
        26: 'Fe',
        79: 'Au'}

symbol_to_atno = dict([[v,k] for k,v in atno_to_symbol.items()])

VdW_Radii = \
        {1: 0.37,
         2: 0.32,
         6: 0.77,
         7: 0.75,
         8: 0.73,
        16: 1.02,
        26: 0.75,
        79: 1.44
         }

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


def get_distance(obj1, obj2):
    if isinstance(obj1, Atom):
        pos1 = obj1.pos
    if isinstance(obj2, Atom):
        pos2 = obj2.pos
    return np.linalg.norm(pos1 - pos2)


def find_shortest_path(graph, start, end, path=[]):
    path = path + [start]
    if start == end:
        return path
    if not graph.has_key(start):
        return None
    shortest = None
    for node in graph[start]:
        if node not in path:
            newpath = find_shortest_path(graph, node, end, path)
            if newpath:
                if not shortest or len(newpath) < len(shortest):
                    shortest = newpath
    return shortest


def get_angle(vec1, vec2):
    angle = np.arccos(np.inner(vec1, vec2)/
                      (np.linalg.norm(vec1)*np.linalg.norm(vec2)))
    return np.rad2deg(angle)


#def find_all_paths(graph, start, end, path=[]):
    #path = path + [start]
    #if start == end:
        #return [path]
    #if not graph.has_key(start):
        #return []
    #paths = []
    #for node in graph[start]:
        #if node not in path:
            #newpaths = find_all_paths(graph, node, end, path)
            #for newpath in newpaths:
                #paths.append(newpath)
    #return paths


class AttrDict(dict):
    """ Class for accessing dict keys like an attribute  """
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


class Graph(object):
    def __init__(self, name='', adjList=None):
        self.name = name
        #self.vertices = {} # vertex list
        self.vertices = [] # vertex list
        self.edges = [] # edge list
        if adjList:
            self.origAdjList = adjList
            self.adjacencyList = self.load_graph(adjList)

    def load_graph(self, graph):
        for v in graph:
            self.add_vertex(v)

        adjacencyList = {}
        for u in graph:
            for vertex1 in self.vertices:
                if vertex1.vertex is u:
                    adjacencyList[vertex1] = []
                    for i, v in enumerate(graph[u]):
                        for vertex2 in self.vertices:
                            if vertex2.vertex is v:
                                adjacencyList[vertex1].append(vertex2)
                                self.add_edge(vertex1, vertex2)
        return adjacencyList

    def add_vertex(self, v):
        if isinstance(v, Vertex):
            vertex = v
        else:
            vertex = Vertex(v)

        self.vertices.append(vertex)
        #self.vertices[vertex] = vertex.vertex
        #return vertex

    def add_edge(self, u, v):
        self.edges.append(Edge(u, v))

    #def breadth_first_search(self, start):
        #Q = [start]
        #D = {node:'inf' for node in self.adjacencyList}
        #D[start] = 0
        #T = []
        #while len(Q) > 0:
            #u = Q.pop(0)
            #for v in graph[u]:
                #if D[v] == 'inf':
                    #D[v] = D[u]+1
                    #Q.append(v)
                    #T.append([u, v])
        #return (D, T)

    def breadth_first_search(self, start):
        T = Tree(adjList=self.origAdjList)
        #vertices = self.adjacencyList.copy()
        #T.add_vertex(v=vertices.pop(start), lvl=0, parent=None)
        for vertex in T.vertices:
            if vertex.vertex is start.vertex:
                T.root = vertex
                vertex.level = 0
                vertex.parent = None
                start = vertex
            else:
                vertex.level = 'inf'
                vertex.parent = None
        #D = {node:'inf' for node in self.adjacencyList}
        #D[start] = 0
        Q = [start]
        while len(Q) > 0:
            u = Q.pop(0)
            for v in T.adjacencyList[u]:
                if v.level == 'inf':
                    v.level = u.level + 1
                    v.parent = u
                    #T.add_vertex(v)
                    T.add_edge(u, v)
                    Q.append(v)
        return T

    #def depth_first_search(self):
        #for vertex in self.vertices:
            #vertex.level = 'inf'
            #vertex.parent = None
        #time = 0
        #for vertex in self.vertices:
            #if vertex.color  is 'inf':
                #DFS_visit()

    #def DFS_visit(self, v):
        #pass

    def find_path(graph, start, end, path=[]):
        path = path + [start]
        if start == end:
            return path
        if not graph.has_key(start):
            return None
        for node in graph[start]:
            if node not in path:
                newpath = find_path(graph, node, end, path)
                if newpath: return newpath
        return None

    def find_all_paths(self, start, end, path=[]):
        path = path + [start]
        if start == end:
            return [path]
        if not self.adjacencyList.has_key(start):
            return []
        paths = []
        for node in self.adjacencyList[start]:
            if node not in path:
                newpaths = self.find_all_paths(node, end, path)
                for newpath in newpaths:
                    paths.append(newpath)
        return paths


class Vertex(object):
    def __init__(self, v=None):
        self.vertex = v

    def __repr__(self):
        return '%s' % (self.vertex,)


class Edge(object):
    def __init__(self, u=None, v=None):
        self.edge = (u, v)

    def __repr__(self):
        return '%s--%s' % (self.edge[0], self.edge[1])


class Tree(Graph):
    def __init__(self, adjList=None):
        Graph.__init__(self, name=None, adjList=adjList)
        self.root = None

    def add_vertex(self, v):
        if isinstance(v, Vertex):
            vertex = TreeVertex(v.vertex)
        else:
            vertex = TreeVertex(v)

        self.vertices[vertex] = vertex.vertex

    def get_length_paths(self, n):
        result = []
        for vertex in self.vertices:
            if vertex.level is n-1:
                a_path = []
                in_vertex = vertex
                #while in_vertex.parent:
                for i in range(n):
                    a_path.append(in_vertex.vertex)
                    in_vertex = in_vertex.parent
                result.append(a_path)
        return result


class TreeVertex(Vertex):
    def __init__(self, v, level=None, parent=None):
        self.vertex = v
        self.level = level
        self.parent = parent


class Atom(object):

    """ Atom class. Must be instantiated as:
        an_atom = Atom(symbol, x, y, z)
        where symbol is either the atomic number or the atom name.
        x, y and z are the atomic coordinates."""

    last_id = count()

    def __init__(self, symbol, pos, label=None, siesta_move=(1,1,1)):

        #self.pos = Point(x, y, z)
        self.set_position(pos)
        self.siesta_move = siesta_move

        if re.match(r'[a-zA-Z]+.*', symbol):
            self.symbol = symbol
            atno = None
        else:
            atno = symbol
        if not self.symbol:
            self.symbol = self.atno_to_symbol()
        if not atno:
            atno = self.symbol_to_atno()
        self.atno = atno

        self.index = self.last_id.next()
        if label:
            self.label = label
        else:
            self.label = str(self.index)

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
        return '{label}-{sym}'.format(sym=self.symbol, label=self.label)
        #return '{name}({sym!r}, {pos}, {index})'.format(
            #name=self.__class__.__name__, sym=self.symbol,
            #pos=self.pos, index=self.index)


class Molecule(object):

    last_id = count()

    def __init__(self, name='Generic', atoms=None, fromFile=None):

        if name == 'Generic':
            self.name = 'Generic' + str(self.last_id.next())
        else:
            self.name = name

        if atoms:
            self.atoms = atoms
        else:
            self.atoms = []

        if fromFile:
            self.load_from_file(fromFile)

        self.get_conections()
        #self.get_structure()

    def update_mol(self):
        for i, atom in enumerate(self.atoms):
            atom.label = str(i)

    def add_atom(self, atom):
        self.atoms.append(atom)
        self.update_mol()

    def remove_atom(self, i):
        print 'Deleting {atom}'.format(atom=self.atoms[i])
        self.atoms.pop(i)
        self.update_mol()

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

        self.atoms = [Atom(namelst[i], pos) for i, pos in enumerate(posarr)]

        self.update_mol()

    def write_to_file(self, filename='default.xyz', fmt='xyz'):
        with open(filename, 'w') as fout:
            fout.write('{tot_at:5d}\n\n'.format(tot_at=len(self.atoms)))
            for atom in self.atoms:
                fout.write('{sym:3s}{x:13.8f}{y:13.8f}{z:13.8f}\n'.format(sym=atom.symbol,
                                                                        x=atom.pos[0],
                                                                        y=atom.pos[1],
                                                                        z=atom.pos[2]))

    def get_geometry_center(self):
        asum = np.zeros(3)
        for atom in self.atoms:
            asum += atom.pos
        asum /= len(self.atoms)
        return asum

    def get_species(self):
        return list(set([atom.symbol for atom in self.atoms]))
        #return {i+1:specie for i, specie in enumerate(species)}
        #return OrderedDict({i+1:specie for i, specie in enumerate(species)})

    def get_conections(self):
        conections = {}
        for atom in self.atoms:
            lconections = []
            ddist = self.get_closest_atom(atom, dist=True)
            for close_atom in ddist:
                pair_dist = VdW_Radii[atom.atno] + VdW_Radii[close_atom.atno]
                if ddist[close_atom] < pair_dist:
                    lconections.append(close_atom)
            conections[atom] = lconections
        self.conections = conections
        return conections

    def get_structure(self):
        G = Graph(adjList=self.conections)
        bonds = []
        angles = []
        dihedrals = []
        not_bonded = []
        for vertex1 in G.vertices:
            for vertex2 in G.vertices:
                conections = G.find_all_paths(vertex1, vertex2)
                if len(conections[0]) > 1:
                    for conection in conections:
                        if len(conection) == 2:
                            bonds.append(Bond(conection))
                        elif len(conection) == 3:
                            angles.append(Angle(conection))
                        elif len(conection) == 4:
                            dihedrals.append(Dihedral(conection))
                        else:
                            not_bonded.append(conection)

        #for vertex in G.vertices:
            #T = G.breadth_first_search(vertex)
            #import pdb; pdb.set_trace()  # XXX BREAKPOINT
            #bonds += T.get_length_paths(2)
            #angles += T.get_length_paths(3)
            #dihedrals += T.get_length_paths(4)
        #for i, f_atom in enumerate(self.atoms):
            #for s_atom in self.atoms[i+1:]:
                #conect = find_shortest_path(self.conections, f_atom, s_atom)
                #if conect:
                    #if len(conect) == 2:
                        #bonds.append(conect)
                    #elif len(conect) == 3:
                        #angles.append(Angle(conect))
                    #elif len(conect) == 4:
                        #dihedrals.append(Dihedral(conect))
                    #else:
                        #not_bonded.append(conect)
        for type_ in [bonds, angles, dihedrals, not_bonded]:
            for conection in type_:
                import pdb; pdb.set_trace()  # XXX BREAKPOINT
                if conection[::-1] in type_:
                    type_.remove(conection)
        self.bonds = bonds
        self.angles = angles
        self.dihedrals = dihedrals
        self.not_bonded = not_bonded

    def get_closest_atom(self, atom, dist=False):
        rest = self.atoms[:]  # Make a copy.
        rest.remove(atom)
        closest_list = {an_atom: get_distance(atom, an_atom) for an_atom in rest}
        if dist:
            return closest_list
        else:
            return sorted(closest_list)

    def __add__(self, an_object):

        if isinstance(an_object, (self.__class__, Selection)):
            #self.atoms += an_object.atoms
            for atom in an_object:
                self.add_atom(atom)

        elif isinstance(an_object, Atom):
            self.add_atom(an_object)

        elif isinstance(an_object, (types.ListType, np.array)):
            for atom in self.atoms:
                atom + an_object

        else:
            raise IOError('Not addition supported of {name} object with \
                            "{obj}"'.format(name=self.__class__.__name__,
                            obj=an_object))
        return self

    def __sub__(self, an_object):

        if isinstance(an_object, (self.__class__, Selection)):
            #self.atoms += an_object.atoms
            for i, orig_at in enumerate(self.atoms):
                for remove_at in an_object:
                    if orig_at.index == remove_at.index:
                        self.remove_atom(i)

        elif isinstance(an_object, Atom):
            for i, orig_at in enumerate(self.atoms):
                if orig_at.index == an_object.index:
                    self.remove_atom(i)

        elif isinstance(an_object, (types.ListType, np.array)):
            for atom in self.atoms:
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
        return self.atoms[index]
        #return self.atoms[index-1]

    def __len__(self):
        return len(self.atoms)

    def __repr__(self):
        return '{name}({atoms})'.format(name=self.__class__.__name__,
                atoms=self.atoms)


class Conection(object):

    def get_atoms_from_list(self, list_):

        self.atoms = []

        if len(list_) is 1:
            list_ = list_[0]

        #if type_ is 'Bond':
            #l = 2
        #elif type_ is 'Angle':
            #l = 3
        #elif type_ is 'Dihedral':
            #l = 4

        for elem in list_:
            if isinstance(elem, Atom):
                self.atoms.append(elem)
            elif isinstance(elem, Vertex):
                self.atoms.append(elem.vertex)

    def __getitem__(self, index):
        return self.atoms[index]

    def __len__(self):
        return len(self.atoms)

    def __contains__(self, x):
        return x in self.atoms

class Bond(Conection):
    def __init__(self, *args):

        self.get_atoms_from_list(args)

        vec1 = self.atoms[0].pos - self.atoms[1].pos

        self.distance = get_distance(self.atoms[0], self.atoms[1])

    def set_angle(self):
        pass

    def __add__(self):
        pass

    def __radd__(self):
        pass

    def __sub__(self):
        pass

    def __rsub__(self):
        pass

    def __repr__(self):
        return '{name} {at1} {at2} = {distance} Angstroms'.format(name=self.__class__.__name__,
                at1=self.atoms[0], at2=self.atoms[1], distance=self.distance)


class Angle(Conection):
    def __init__(self, *args):

        self.get_atoms_from_list(args)

        #if len(args) is 1:
            #self.get_atoms_from_list(args[0])
        #elif len(args) is 3:
            #self.get_atoms_from_list(args)
        #else:
            #raise IOError('A list of three atom or three atoms needed')

        vec1 = self.atoms[0].pos - self.atoms[1].pos
        vec2 = self.atoms[2].pos - self.atoms[1].pos

        self.angle = get_angle(vec1, vec2)

    def set_angle(self):
        pass

    def __add__(self):
        pass

    def __radd__(self):
        pass

    def __sub__(self):
        pass

    def __rsub__(self):
        pass

    def __repr__(self):
        return '{name} {at1} {at2} {at3} = {angle} º'.format(name=self.__class__.__name__,
                at1=self.atoms[0], at2=self.atoms[1], at3=self.atoms[2], angle=self.angle)


class Dihedral(Conection):
    def __init__(self,*args):

        self.get_atoms_from_list(args)

        #if len(args) is 1:
            #if isinstance(args[0], types.ListType):
                #self.atom1, self.atom2, self.atom3, self.atom4 = args[0]
        #elif len(args) is 4:
            #for elem in args:
                #if not isinstance(elem, Atom):
                    #raise TypeError('Not Atom Instance')
            #self.atom1, self.atom2, self.atom3, self.atom4 = args
        #else:
            #raise IOError('A list of four atom or four atoms needed')

        vec1 = self.atoms[0].pos - self.atoms[1].pos
        vec2 = self.atoms[2].pos - self.atoms[1].pos
        vec3 = - vec2
        vec4 = self.atoms[3].pos - self.atoms[2].pos

        n1 = np.cross(vec1, vec2)
        n2 = np.cross(vec3, vec4)

        self.dihedral = get_angle(n1, n2)

    def __add__(self):
        pass

    def __radd__(self):
        pass

    def __sub__(self):
        pass

    def __rsub__(self):
        pass

    def __repr__(self):
        return '{name} {at1} {at2} {at3} {at4} = {dihedral} º'.format(
            name=self.__class__.__name__, at1=self.atoms[0], at2=self.atoms[1],
            at3=self.atoms[2], at4=self.atoms[3], dihedral=self.dihedral)


class Selection(Molecule):

    def __init__(self, a_system):

        atoms = []

        if isinstance(a_system, Atom):
            atoms = a_system

        elif isinstance(a_system, Molecule):
            for atom in a_system:
                atoms.append(atom)

        elif isinstance(a_system, System):
            for atom in a_system:
                atoms.append(atom)

        self.atoms = atoms
        self.orig_alist = atoms

    def _filter(self, func):
        return filter(func, self.atoms)

    def remove_from_list(self, an_element):
        if an_element in self.atoms:
            self.atoms.remove(an_element)

    def by_symbol(self, symbol):
        self.atoms = self._filter(lambda a: a.symbol == symbol)
        return self

    def by_xrange(self, xmin, xmax):
        self.atoms = self._filter(lambda a: xmin <= a.pos[0] <=
                                      xmax)
        return self

    def by_yrange(self, ymin, ymax):
        self.atoms = self._filter(lambda a: ymin <= a.pos[1] <=
                                      ymax)
        return self

    def by_zrange(self, zmin, zmax):
        self.atoms = self._filter(lambda a: zmin <= a.pos[2] <=
                                      zmax)
        return self

    def sphere(self, radius, pos=None, center_atom=None):
        if pos.any():
            center = pos
        elif center_atom:
            center = center_atom.pos
        else:
            center = np.zeros(3)
        self.atoms = self._filter(lambda a: np.linalg.norm(a.pos - center) <=
                                      radius**2)
        return self

    def update_sel(self):
        for m_atom in self.atoms:
            for i, o_atom in enumerate(self.orig_alist):
                if m_atom.index == o_atom.index:
                    self.orig_alist[i] = m_atom
        self.atoms = self.orig_alist
        return self.atoms


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
        #zmaMolecule = []
        #for i, atom in enumerate(mol):
            #if i is 0:
                #zmaMolecule.append([atom.specieNo, 0, 0, 0, atom.pos[0],
                                    #atom.pos[1], atom.pos[2],
                                    #atom.siesta_move[0], atom.siesta_move[1],
                                    #atom.siesta_move[2]])
            #if i is 1:
                #zmaMolecule.append([atom.specieNo, 1, 0, 0]<++>)<++>

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


