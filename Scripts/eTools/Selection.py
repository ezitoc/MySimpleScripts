#!/usr/bin/env python
#-*- coding: utf-8 -*-

import eTools


class Selection():

    def __init__(self, a_system):

        atom_list = []

        if isinstance(a_system, eTools.Atom):
            atom_list = a_system

        elif isinstance(a_system, eTools.Molecule):
            for atom in a_system:
                atom_list.append(atom)

        elif isinstance(a_system, eTools.System):
            for atom in a_system:
                atom_list.append(atom)

        self.atom_list = atom_list
        self.orig_alist = atom_list

    def __repr__(self):
        return '{name}({atoms})'.format(name=self.__class__.__name__,
                atoms=self.atom_list)

    def __add__(self, an_object):
        if isinstance(an_object, list):
            if len(an_object) == 3:
                x, y, z = an_object
                for atom in self.atom_list:
                    atom.translate(x, y, z)
            else:
                raise IOError('Not addition supported of {name} object with \
                              "{obj}"'.format(name=self.__class__.name,
                              obj=an_object))
        return self

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

    def update_sel(self):
        for m_atom in self.atom_list:
            for i, o_atom in enumerate(self.orig_alist):
                if m_atom.index == o_atom.index:
                    self.orig_alist[i] = m_atom
        self.atom_list = self.orig_alist
        return self.atom_list

    def move(self, x, y, z):
        for atom in self.atom_list:
            atom.translate(x, y, z)
        return self

    def __len__(self):
        return len(self.atom_list)

    def __getitem__(self, index):
        return self.atom_list[index]


if __name__ == "__main__":
    main()
