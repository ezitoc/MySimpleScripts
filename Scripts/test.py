#!/usr/bin/env python

#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
#
#* File Name : test.py
#
#* Purpose :
#
#* Creation Date : 25-04-2013
#
#* Last Modified : Thu 25 Apr 2013 01:23:33 PM ART
#
#* Created By :  Ezequiel Castillo
#
#_._._._._._._._._._._._._._._._._._._._._.

import itertools
import numpy as np
from matplotlib import pyplot as plt

class Base(object):
    def __init__(self, a, multiPlot=True, numColGraph=None, numRowGraph=None,
                 figSize=None, DPI=None, num=None):

        self.a = a
        self.x = np.linspace(0, 5)

        if multiPlot:
        self.nCG = numColGraph
            self.nRG = numRowGraph
        else:
            self.nCG = 1
            self.nRG = 1

        if figSize and DPI:
            self.thePlot = plt.figure(figsize=figSize, dpi=DPI)

        if num == 0:
            self.plotId = itertools.count(1)

    def createPlot1(self):
        y = self.x**(a/2)
        self.thePlot.add_subplot(self.nRG, self.nCG, next(self.plotId))
        plt.plot(self.x, y, label=str(self.a)+'/2')

    def createPlot2(self):
        y = self.x**a
        self.thePlot.add_subplot(self.nRG, self.nCG, next(self.plotId))
        plt.plot(self.x, y, label=self.a)

    def createPlot3(self):
        y = self.x**(2*a)
        self.thePlot.add_subplot(self.nRG, self.nCG, next(self.plotId))
        plt.plot(self.x, y, label=str(self.a)+'*2')


if __name__ == "__main__":

    A = np.linspace(0, 2, 5)

    for i, a in enumerate(A):
        Instance = Base(a, numColGraph=3, numRowGraph=len(A),
                 figSize=(12,10), DPI=100, num=i)
        Instance.createPlot1()
        Instance.createPlot2()
        Instance.createPlot3()

    plt.show()
