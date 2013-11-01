#!/usr/bin/env python
#-*- coding: utf-8 -*-

import eTools


class Selection(eTools.Molecule):

    def __init__(self, a_system, atom_list=[]):
        if isinstance(a_system, eTools.Atom):
            self.atom_list = a_system

        if isinstance(a_system, eTools.Molecue):
            for atom in a_system:
                Molecule.add_atom(atom)

        if isinstance(a_system, eTools.System):
            for atom in a_system:
                Molecule.add_atom(atom)

    #def zrange(self):
        #for atom in atom_list:


if __name__ == "__main__":
    main()
