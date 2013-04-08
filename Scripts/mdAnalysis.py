#!/usr/bin/env python

#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
#
#* File Name : mdAnalysis.py
#
#* Purpose :
#
#* Creation Date : 06-03-2013
#
#* Last Modified : Mon 11 Mar 2013 02:30:21 PM ART
#
#* Created By :  Ezequiel Castillo
#
#_._._._._._._._._._._._._._._._._._._._._.

import os


class Base(object):
    """General functions definitios"""

    def __init__(self, coordFile):
        self.coordFile = os.path.basename(coordFile)
        self.baseName, self.ext = os.path.splitext(self.coordFile)
        print "%s %s" % (self.baseName, self.ext)
        exit(0)



    def getCoords(xyzFileName=None):
        """Get coordinates as a list of list, i.e as [[x1,y1,z1],[x2,y2,z2],etc]"""
        if self.ext == '.xyz':
            this
        elif self.ext:
            pass
    pass

class Atom(object):
    """Class containing atom atributes"""
    pass

if __name__ == '__main__':
    systemA = Base('initial.xyz')
    xyzCoord = systemA.getCoords()
    print xyzCoords
