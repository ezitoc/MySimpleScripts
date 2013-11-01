#!/usr/bin/env python
#-*- coding: utf-8 -*-

import eTools


class Selection():

    def __init__(self, a_system, atom_list=[]):

        if isinstance(a_system, eTools.Atom):
            atom_list = a_system

        elif isinstance(a_system, eTools.Molecule):
            for atom in a_system:
                atom_list.append(atom)

        elif isinstance(a_system, eTools.System):
            for atom in a_system:
                atom_list.append(atom)

        self.atom_list = atom_list

    def remove_from_list(self, an_element):
        if an_element in self.atom_list:
            self.atom_list.remove(an_element)

    def range(self, axmin, axmax, coord='z'):
        if coord == 'x':
            ax = 0
        elif coord == 'y':
            ax = 1
        elif coord == 'z':
            ax = 2
        else:
            raise IOError
        r_list = []
        for atom in self.atom_list:
            pos = atom.get_position()[ax]
            if pos > axmin and pos < axmax:
                r_list.append(atom)
        self.atom_list = r_list

    def sel_by_type(self, at_type):
        r_list = []
        for atom in self.atom_list:
            if atom.symbol == at_type:
                r_list.append(atom)
        self.atom_list = r_list

    def __getitem__(self, index):
        return self.atom_list[index]


if __name__ == "__main__":
    main()
