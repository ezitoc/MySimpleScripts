#!/usr/bin/env python

#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
#
#* File Name : calcE.py
#
#* Purpose :
#
#* Creation Date : 11-03-2013
#
#* Last Modified : Fri 22 Mar 2013 06:46:09 PM ART
#
#* Created By :  Ezequiel Castillo
#
#_._._._._._._._._._._._._._._._._._._._._.

import sys
import math
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
from scipy.optimize import curve_fit

rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


def checkCorr(Xvalues, Yvalues):

    def func(t, a, b, c):
        """Fitting exponential function"""
        return a*np.exp(-t/b)+np.exp(-t/c)

    # Ignore overflow error (for curve_fit)
    np.seterr(over='ignore')

    Aprom = np.average(Yvalues)
    A = Yvalues-Aprom
    corr = np.correlate(A, A, "full")

    # Auto Correlation Function
    ACF = corr[corr.size/2:]/(np.var(A)*len(A))
    print len(Xvalues)


    # Exponential fit with Scipy's curve_fit
    popt, pcov = curve_fit(func, Xvalues, ACF)
    print popt, pcov

    # Tau, integration of correlation function
    tau = popt[0]*popt[1]+popt[2]
    print tau

    #print popt[0], popt[1], popt[2], popt[3]
    plt.plot(Xvalues, ACF, 'o', label="ACF")
    plt.plot(Xvalues, func(Xvalues, *popt), label="Fitted Curve")
    plt.legend(loc='upper right')
    plt.xlabel('Tiempo [ps]')
    plt.ylabel('$C(t)$')
    plt.axis([0,5*tau,-0.05,1])
    plt.show()
    plt.savefig('corr.pdf', format='pdf', bbox_inches='tight')



class calcE(object):
    """ Super class"""

    def __init__(self, vprom=None, alpha=None, vmin=None, Ecut=None, freq=None,
                 T=None, biasFile=None, enerFile=None, tempFile=None,
                 smoothFactor=None, tau=None):

        if Ecut and freq:
            raise NameError('Energy cut and scape frequence not at the same time!')
        if smoothFactor and tau:
            raise NameError('Smooth factor and tau not at the same time!')

        self.kB = 8.6173324*10**-5 # eV/K
        self.T = T # K
        self.beta = 1./(self.kB*self.T)
        self.vprom = vprom
        self.alpha = alpha
        self.vmin = vmin
        self.freq = freq
        self.biasFile = biasFile
        self.VSinBias = np.loadtxt(self.biasFile, usecols=[0])
        self.deltaBias = np.loadtxt(self.biasFile, usecols=[1])
        #self.VBias = np.loadtxt(self.biasFile, usecols=[2])
        self.tempFile = tempFile
        self.enerFile = enerFile
        self.enerPot = np.loadtxt(self.enerFile, usecols=[0])
        self.timeNoBoost = np.loadtxt(self.tempFile, usecols=[0])
        self.timeBoost = np.loadtxt(self.tempFile, usecols=[0])
        self.smoothFactor = smoothFactor
        self.tau = tau

    def evalE(self):
        self.ener = (self.vprom + self.vmin * (self.alpha - 1) - self.kB *\
                    self.temp * math.log(self.freq)) / self.alpha
        return self.ener

    def diffE(self):
        self.diffE = self.ener - self.vmin
        return self.diffE

    def evalExVb(self):
        self.ExVb = self.vprom + (self.vmin - self.ener) * (self.alpha - 1)
        return self.ExVb

    def evalVbProm(self):
        self.VbProm = np.average(self.VBias)
        return self.VbProm

    def evalVrepesado(self):
        """ We want to get <Vb>b"""
        #VbSN = self.VSinBias * math.exp(self.beta*self.deltaBias)
        VbSN = np.array([ self.VSinBias[i]*math.exp(self.beta*self.deltaBias[i]) for i in range(len(self.VSinBias))  ])
        Norm = np.array([math.exp(self.beta*self.deltaBias[i]) for i in range(len(self.VSinBias))])
        Vbb = np.sum(VbSN)/np.sum(Norm)
        checkCorr(self.timeNoBoost, self.VSinBias)
        print Vbb
        sys.exit()
        plt.plot(self.timeNoBoost, self.enerPot, 'o-g', label='Potencial no boosteado', markersize=1)
        plt.legend()
        plt.show()

    def time(self):
        A = np.vstack([self.timeNoBoost, np.ones(len(self.timeNoBoost))]).T
        S = np.linalg.lstsq(A, self.timeBoost)
        a, b = S[0]
        res = S[1][0]
        plt.plot(self.timeNoBoost, self.timeBoost, 'o-g', label='alpha='+str(self.alpha), markersize=1)
        plt.plot(timeBoostx, a*timeBoostx + b, 'r')
        plt.legend(loc='upper left')
        plt.show()
        plt.savefig(str(self.alpha)+'.pdf', format='pdf', bbox_inches='tight')

    def plotPot(self):
        smoothX, smoothY = self.smoothData(self.timeBoost, self.enerPot)
        plt.plot(self.timeBoost, self.enerPot)
        plt.plot(smoothX, smoothY)
        plt.show()

    def smoothData(self, x, y):
        # Check
        if len(x)!=len(y):
            print 'ERROR: Length of Xdata not equal to Ydata'
        elif self.smoothFactor:
            #Check if the problem can be divided in equal spaces
            if len(x)*self.smoothFactor % 1/self.smoothFactor != 0:
                print 'ERROR: No divisibility between Xdata and Ydata'
            segNo = int(len(x) * self.smoothFactor) # Numero de segmentos
            segQ = int(1 / self.smoothFactor) # Numero de valores por segmento
            x.shape=(segNo, segQ)
            y.shape=(segNo, segQ)

        elif self.tau:
            # Tau in picosecond
            size = int(self.tau/(x[1]-x[0]))
            xa = [x[i:i+size] for i in range(0, len(x), size)]
            ya = [y[i:i+size] for i in range(0, len(y), size)]

        newX = []
        for list in xa:
            newX.append(np.average(list))

        newY = []
        for list in ya:
            newY.append(np.average(list))

        return np.array(newX), np.array(newY)


if __name__ == "__main__":
    """ A continuacion se definen las variables a utilizar por el programa. """


    #f = calcE(vprom=-1540.45, alpha=0.8, vmin=-1554.63, freq=1*10-3,
              #temp=300, biasFile='bias.dat', enerFile='ener.dat',
              #tempFile='temp.dat', smoothFactor=0.001)

    g = calcE(vprom=-1540.45, alpha=0.8, vmin=-1554.63, freq=1*10-3,
              T=300, biasFile='bias.dat', enerFile='ener.dat',
              tempFile='temp.dat', tau=3.52)

    #energy = f.evalE()
    #diffener = f.diffE()
    #ExVb = f.evalExVb()
    #VbProm = f.evalVbProm()
    ##f.plotData()
    V = g.evalVrepesado()
    #V = f.plotPot()
    #E = g.plotPot()

    #checkCorr('temp.dat', 'ener.dat', 0, 0)

    #print 'Energy = %g eV' % (energy)
    #print 'DiffEnergy = %g eV' % (diffener)
    #print '<Vb>b = %g eV' % (VbProm)
    #print 'Expected <Vb>b = %g eV' % (ExVb)
