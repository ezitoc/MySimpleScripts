#!/usr/bin/env python

#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
#
#* File Name : calcE.py
#
#* Purpose :
#
#* Creation Date : 11-03-2013
#
#* Last Modified : Tue 23 Apr 2013 04:19:12 AM ART
#
#* Created By :  Ezequiel Castillo
#
#_._._._._._._._._._._._._._._._._._._._._.

import sys
import os
import re
import math
import pdb
import shutil
import itertools
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
from scipy.optimize import curve_fit

rc('font',**{'family':'serif','serif':['Palatino'],'size':10})
rc('text', usetex=True)
#rc('lines', linewidth='1.0')
#rc('axes', linewidth='0.5')
#rc('axes', style='plain')
#rc('figure', **{'subplot.bottom':0.1})

np.seterr(over='print')

dictAtomMass = \
{
'Au':196.967,
'C':12.011,
'H':1.008,
'S':32.066}

def count(firstval=0, step=1):
    x = firstval
    while True:
        yield x
        x += step

def pairwise(iterable):
    itnext = iter(iterable).next
    while True:
        yield itnext( ), itnext( )

def dictFromSequence(seq):
    return dict(pairwise(seq))

def replace_words(text, word_dic):
    """
    take a text and replace words that match a key in a dictionary with
    the associated value, return the changed text
    """
    rc = re.compile('|'.join(map(re.escape, word_dic)))
    def translate(match):
        return word_dic[match.group(0)]
    return rc.sub(translate, text)

def multiple_subst(infile, replace_pattern):
    # read the file
    match = re.search(r'.*/(.*)$', infile)
    filename = match.group(1)
    fin = open(infile, "r")
    str1 = fin.read()
    fin.close()
    # add 'replace pattern' : 'replacement' to dictionary
    dict_make = dictFromSequence(replace_pattern)
    # call the function and get the changed text
    str2 = replace_words(str1, dict_make)
    # write changed text back out
    fout = open(os.path.join(os.getcwd(),filename), "w")
    fout.write(str2)
    fout.close()


def checkCorr(Xvalues, Yvalues, showPlot=False):

    def func(t, a, b, c, d, e):
    #def func(t, a, b, c):
        """Fitting exponential function"""
        #return a*np.exp(-t/b)+np.exp(-t/c)
        return np.exp(-t/a)+np.exp(-t/b)+np.exp(-t/c)+np.exp(-t/d)+np.exp(-t/e)

    # Ignore overflow error (for curve_fit)

    Aprom = np.average(Yvalues)
    A = Yvalues-Aprom
    corr = np.correlate(A, A, "full")

    # Auto Correlation Function
    ACF = corr[corr.size/2:]/(np.var(A)*len(A))


    # Exponential fit with Scipy's curve_fit
    popt, pcov = curve_fit(func, Xvalues, ACF, maxfev = 10000)

    # Tau, integration of correlation function
    tau = popt[0]+popt[1]+popt[2]+popt[3]+popt[4]

    if showPlot:
        plt.plot(Xvalues, ACF, 'o', label="ACF")
        plt.plot(Xvalues, func(Xvalues, *popt), label="Fitted Curve")
        plt.legend(loc='upper right')
        plt.xlabel('Tiempo [ps]')
        plt.ylabel('$C(t)$')
        plt.axis([0,30*tau,-0.05,1])
        #plt.show()
        #plt.savefig('corr.pdf', format='pdf', bbox_inches='tight')

    return tau

def makeBarsFromHisto(bins, widthFactor):
    center = (bins[:-1]+bins[1:])/2
    width = (bins[1]-bins[0])*widthFactor
    return center, width


class createInputs(object):

    def __init__(self, projectName='project', subFolders=[], filesToCopy=[],
                 modFilesAndPattern=None, pattern=None, ignoreFiles=None):
        self.rootDir = os.path.abspath(os.curdir)
        self.projectName = projectName
        self.subFolders = subFolders
        #self.modFiles = modFiles
        self.modFiles = [modFile for modFile in modFilesAndPattern]
        self.filesToCopy = filesToCopy
        self.modFilesAndPattern = modFilesAndPattern
        self.ignoreFiles = ignoreFiles
        #self.pattern = pattern

    def selectFiles(self):
        """Return a list of files that won't be modified"""
        if not self.filesToCopy:
            #Copy all files by default.
            filesList = [f for f in os.listdir(self.rootDir) if
                         os.path.isfile(os.path.join(self.rootDir, f))]
        else:
            filesList = self.filesToCopy
        if self.ignoreFiles:
            filesListMod = [f for f in filesList if f not in self.ignoreFiles]
            filesList = filesListMod
        elif self.modFiles:
            filesListMod = [f for f in filesList if f not in self.modFiles]
            filesList = filesListMod
        return filesList

    def createFolders(self, foldersSuffix=None):

        """Create subfolders containing all files declared in the filesToCopy
        option (all files included, otherwise) and ignoring those declared at the
        modFiles optioni.
        At the moment you can only include one modFile.
        """

        #if suffixNo and suffixName:
            #print "ERROR: Two suffix types declared. Declare only one."
            #sys.exit(1)

        #if suffixNo or suffixName:
            #if suffixNo:
                #self.suffixList = range(suffixNo)
            #else:
                #self.suffixList = suffixName
        #else:
            #print "ERROR: No suffix declared."
            #sys.exit(1)

        if not os.path.isdir(self.projectName):
            os.mkdir(self.projectName)
        else:
            print 'WARN: folder "%s" already exists.' % (self.projectName)
        os.chdir(self.projectName)

        selectedFiles = self.selectFiles()

        # TODO: Check incopatibility between modFiles and selectedFiles

        for i, subFolder in enumerate(self.subFolders):
            replacePatt = subFolder
            if not os.path.isdir(subFolder):
                os.mkdir(subFolder)
            else:
                print 'WARN: folder "%s" already exists.' % (subFolder)

            os.chdir(subFolder)

            for file in selectedFiles:
                shutil.copy(os.path.join(self.rootDir, file), os.curdir)

            for modFile in self.modFilesAndPattern:
                pattern = self.modFilesAndPattern[modFile][i]
                multiple_subst(os.path.join(self.rootDir, modFile), pattern)

            os.chdir(os.pardir)








        #for suffix in self.suffixList:
            #self.folderName = '%s_%s' % (self.foldersName, suffix)
        ## CWD = CWD/projectName/folderName
            #os.mkdir(self.folderName)
            #os.chdir(self.folderName)
            #copySelectedFiles = self.selectFiles()
            #for file in copySelectedFiles:
                #shutil.copy(os.path.join(self.rootDir, file), os.curdir)
        ## CWD = CWD/projectName
            #os.chdir(os.pardir)

        #os.chdir(self.rootDir)


class calcE(object):
    """ Super class"""

    def __init__(self, vprom=None, alpha=None, vmin=None, Ecut=None, freq=None,
                 T=None, biasFile=None, enerFile=None, tempFile=None,
                 smoothFactor=None, tau=None, AMD=False, tMax=None, dt=None,
                 dFrame=None, xyzFile=None, MD_enerFile=None, basedir=os.getcwd(),
                 dFrameXYZ=None, generatePlot=False, filePlotName=None,
                 multiPlot=False, figSize=None, DPI=100, stepId=0,
                 numColGraph=None, numRowGraph=None):

        self.basedir = basedir

        if Ecut and freq:
            raise NameError('Energy cut and scape frequence not at the same time!')
        if smoothFactor and tau:
            raise NameError('Smooth factor and tau not at the same time!')
        if AMD:
            if not tMax and dt and dFrame:
                raise NameError('Acelerated Molecular Dynamics, declare tMax and dt to compute ACF')
            else:
                self.dtSteps = dt * dFrame
                self.timeNoBoost = np.arange(self.dtSteps, tMax+self.dtSteps, self.dtSteps)
        if not vprom and not MD_enerFile:
            raise NameError('Must declare at least one of the following: vprom or MD_enerFile')
        elif MD_enerFile:
            self.MD_enerFile = MD_enerFile
            self.MD_enerPot = np.loadtxt(os.path.join(self.basedir,self.MD_enerFile), usecols=[0])
        if generatePlot:
            self.generatePlot = generatePlot
            if figSize and dpi:
                #NOTE: La misma duda de plotId aplica aca tambien.
                self.thePlot = plt.figure(figsize=figSize, dpi=DPI)
            else:
                NameError('generatePlot option set to true, must declare figSize and dpi')
            if multiPlot:
                if numColGraph and numRowGraph:
                    self.nCG = numColGraph
                    self.nRG = numRowGraph
                    self.numTotGraph = numRowGraph*numColGraph
                else:
                    NameError('MultiPlot option set to True, must declare both numColGraph and numRowGraph')
            else:
                self.nCG = 1
                self.nRG = 1
                self.numTotGraph = self.nRG*self.nCG

            if stepId == 0:
                # NOTE: No se si debe ir plotId solo o con el self.
                self.plotId = count(1)


        self.kB = 8.6173324*10**-5 # eV/K
        self.T = T # K
        self.beta = 1./(self.kB*self.T)
        self.alpha = alpha
        self.vmin = vmin
        self.freq = freq

        if biasFile:
            self.biasFile = biasFile
            self.VSinBias = np.loadtxt(self.biasFile, usecols=[0])
            self.deltaBias = np.loadtxt(self.biasFile, usecols=[1])
            self.VBias = np.loadtxt(self.biasFile, usecols=[2])
        if tempFile:
            self.tempFile = tempFile
            self.timeBoost = np.loadtxt(self.tempFile, usecols=[0])
        if enerFile:
            self.enerFile = enerFile
            self.enerPot = np.loadtxt(self.enerFile, usecols=[0])
            self.enerTot = np.loadtxt(self.enerFile, usecols=[2])
        if dFrameXYZ:
            self.dtStepsXYZ = dt * dFrameXYZ
            self.timeXYZ = np.arange(self.dtStepsXYZ, tMax+self.dtStepsXYZ, self.dtStepsXYZ)
        self.smoothFactor = smoothFactor
        self.tau = tau
        self.xyzFile = xyzFile

        if vprom:
            self.vprom = vprom
        else:
            self.vprom = np.average(self.MD_enerPot)

        # Valores importantes

        self.Ecut = (self.vprom + self.vmin * (float(self.alpha) - 1) - self.kB *\
                    self.T * math.log(self.freq)) / float(self.alpha)
        self.diffE = self.Ecut - self.vmin
        self.ExVb = self.vprom + (self.vmin - self.Ecut) * (self.alpha - 1)


    def histo(self):
        #plt.figure(1, figsize=(8,6))
        self.thePlot.add_subplot(self.nRG, self.nCG, next(self.plotId))
        heights1, bins1, patches = plt.hist(self.enerPot, 30)
        heights2, bins2, patches2 = plt.hist(self.MD_enerPot, 30)
        plt.close()
        plt.title('Alpha = %.2f' % (self.alpha))
        center = (bins1[:-1]+bins1[1:])/2
        # Repesado del histograma
        heightsR = []
        for E, height in zip(center, heights1):
            if E < self.Ecut:
                heightsR.append(height*np.exp((E-self.Ecut)*self.beta*(self.alpha-1)))
            else:
                heightsR.append(height)
        center, width = makeBarsFromHisto(bins=bins2, widthFactor=1)
        plt.bar(center, map(float, heights2)/(np.sum(heights2)*width), width=width, label='DM Comun', color='orange')

        center, width = makeBarsFromHisto(bins=bins1, widthFactor=0.8)
        plt.bar(center, map(float, heights1)/(np.sum(heights1)*(bins1[1]-bins1[0])), width=width, label='Acelerada No Repesada', color='blue')

        center, width = makeBarsFromHisto(bins=bins1, widthFactor=0.6)
        plt.bar(center, heightsR/(np.sum(heightsR)*(bins1[1]-bins1[0])), width=width, label='Acelerada Repesada', color='red')

        locs,labels = plt.xticks()
        plt.xticks(locs, map(lambda x: "%g" % x, locs-min(locs)))
        plt.text(0.92, -0.07, "+%g" % min(locs), fontsize=10, transform = plt.gca().transAxes)
        plt.legend(loc='upper left')
        plt.axvline(x=self.Ecut, linewidth=2, color='g', ls='--')
        #plt.savefig('Histogramas.pdf', format='pdf', bbox_inches='tight', dpi=50)
        #plt.locator_params(nbins=4)
        #plt.ticklabel_format(style='plain')
        #plt.show()
        #pdb.set_trace()

        #plt.show()

    #def evalEcut(self):
        #self.Ecut = (self.vprom + self.vmin * (float(self.alpha) - 1) - self.kB *\
                    #self.T * math.log(self.freq)) / float(self.alpha)
        #return self.Ecut

    #def diffE(self):
        #self.diffE = self.Ecut - self.vmin
        #return self.diffE

    #def evalExVb(self):
        #self.ExVb = self.vprom + (self.vmin - self.Ecut) * (self.alpha - 1)
        #return self.ExVb

    def evalVbProm(self):
        self.VbProm = np.average(self.VBias)
        return self.VbProm

    def plotVOverlaped(self, plot=True):
        """ We want to overlap the V curve of a normal MD with the acelerated one """

        # Arreglo con VbSN[i]=V[i]*exp(beta*DeltaV[i]), i.e., el
        # numerador del valor de expectacion
        VbSN = np.array([ self.VSinBias[i]*np.exp(self.beta*self.deltaBias[i]) for i in range(len(self.VSinBias)) ])

        # Arreglo con normFac[i]=exp(beta*DeltaV[i]), i.e., el
        # denominador del valor de expectacion
        normFac = np.array( [np.exp(self.beta*self.deltaBias[i]) for i in range(len(self.deltaBias))]  )

        # Calculo de correlacion
        tau = checkCorr(self.timeNoBoost, self.VSinBias, showPlot=False)
        decor = 2*tau

        # Cortamos VbSN y normFac de acuerdo a la cantidad de elementos
        # del potencial que entran en un tiempo igual a decor
        size = int(decor/self.dtSteps)
        VBlock = [VbSN[i:i+size] for i in range(0, len(VbSN), size)]
        deltaVBlock = [normFac[i:i+size] for i in range(0, len(normFac), size)]

        pesos = []
        for block in deltaVBlock:
            pesos.append(np.sum(block))

        Vrepes = []
        for i, data in enumerate(VBlock):
            Vrepes.append(np.sum(data)/pesos[i])

        # Tenemos que armar ahora un arreglo con los tiempos de la
        # dinamica acelerada, que lo hacemos igual a:
        #
        # sum_i^N (deltaT*exp(beta*<deltaV(r_i)>))
        #
        # donde <deltaV> es el promedio en el intervalo de tiempo

        deltaVBlock = [self.deltaBias[i:i+size] for i in range(0, len(self.deltaBias), size)]

        thyper = []
        value = 0
        for block in deltaVBlock:
            value = value + decor*np.exp(self.beta*np.average(block))
            thyper.append(value)

        if plot:
            plt.subplot(224)
            self.thePlot.add_subplot(self.nRG, self.nCG, next(self.plotId))
            plt.plot(self.timeNoBoost, self.MD_enerPot, 'o-r', label='Potencial DM comun', markersize=3)
            plt.plot(thyper, Vrepes, 'o-g', label='Potencial DM acelerada', markersize=1)
            plt.xlabel('Tiempo [ps]')
            plt.ylabel('$V$ [eV]')
            plt.legend(loc='lower center')
            #plt.show()

    def compareTimes(self, plot=True):
        """ At this point we want to check the linearity of time vs hyperdynamic time  """
        A = np.vstack([self.timeNoBoost, np.ones(len(self.timeNoBoost))]).T
        S = np.linalg.lstsq(A, self.timeBoost)
        a, b = S[0]
        res = S[1][0]
        if plot:
            plt.subplot(222)
            self.thePlot.add_subplot(self.nRG, self.nCG, next(self.plotId))
            plt.plot(self.timeNoBoost, self.timeBoost, 'o-g', label='alpha='+str(self.alpha), markersize=1)
            plt.plot(self.timeNoBoost, a*self.timeNoBoost + b, 'r')
            plt.legend(loc='upper left')
            #plt.savefig('Times.pdf', format='pdf', bbox_inches='tight', dpi=100)
            #plt.show()
        return a

    def plotPot(self):
        smoothX, smoothY = self.smoothData(self.timeBoost, self.enerPot)
        plt.plot(self.timeBoost, self.enerPot)
        plt.plot(smoothX, smoothY)
        #plt.show()

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

    def inertia(self, plot=True):
        f = open(self.xyzFile, 'r')
        lines = f.readlines()
        f.close()

        totalAtoms = int(lines[0])
        size = totalAtoms+2

        totalFrames = len(lines)/(totalAtoms+2)

        index = []
        for frame in range(totalFrames):
            index.append(frame * (totalAtoms+2))
            index.append(frame * (totalAtoms+2) + 1)

        coordsRaw = np.delete(lines, index)

        size = totalAtoms
        coords = [coordsRaw[i:i+size] for i in range(0, len(coordsRaw), size)]

        I=[]
        for i in range(totalFrames):
            qI = []
            for line in coords[i]:
                element = line.split()[0]
                x, y, z = map(float, line.split()[1:])
                qI.append(dictAtomMass[element]*(x**2+y**2+z**2))
            I.append(np.sum(qI))


        if generatePlot:
            plt.subplot(223)
            self.thePlot.add_subplot(self.nRG, self.nCG, next(self.plotId))
            plt.plot(self.timeXYZ, I, '-' )
            plt.xlabel('Time [ps]')
            plt.ylabel('$I$ [$\\textrm{uma}\\cdot\\AA ^2$]')
            plt.legend()
            #plt.show()

        return np.amax(I)-np.amin(I)

    def escapeFactor(self):
        escFac = 1. - float(np.count_nonzero(self.deltaBias))/float(len(self.deltaBias))
        return escFac

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------


if __name__ == "__main__":
    """ A continuacion se definen las variables a utilizar por el programa. """

    projectName = 'AlphaTEST'

    #vmin = -1554.63
    #alfas = map(str, np.linspace(0.1, 0.5, 6))
    #patternGems = ['$ALPHA$', '$ECUT$'] # El primer elemento de la lista tiene que corresponder con el parametro crudo a variar
    #patternModelo = ['$ALPHA$']

    #Ecut = [calcE(alpha=ALPHA, MD_enerFile='DM_ener.dat', vmin=vmin, T=300,
        #freq=0.001).evalEcut() - vmin for ALPHA in alfas]
    #Ecut = map(str, Ecut)

    #fullPatternGems = []
    #for alfa, ecut in zip(alfas, Ecut):
        #fullPatternGems.append([patternGems[0], alfa, patternGems[1], ecut])

    #fullPatternModelo = [[patternModelo[0], alpha] for alpha in alfas]

    #filesAndPattern = {'gems.gms' : fullPatternGems,
            #'SIDRA.sge' : fullPatternModelo}

    #A = createInputs(projectName=projectName, subFolders=alfas,
                     #modFilesAndPattern=filesAndPattern,
                     #ignoreFiles=['DM_ener.dat', 'calcE.py'])
                     ##filesToCopy=['ener.dat', 'DM_ener.dat'])
    #A.createFolders()
    #pdb.set_trace()

#==========================================================================

    # The current working directory
    #thedir = os.path.abspath(os.curdir)
    rootdir = os.getcwd()


    os.chdir(projectName)

    thedir = os.getcwd()
    dirs = [ float(name) for name in os.listdir(thedir) if
            os.path.isdir(os.path.join(thedir, name)) ]
    dirs.sort()
    dirs.reverse()


    dirs = [element for element in dirs if float(element) > 0.5]

    cmPerGraph = np.array([5.2, 5.2]) # (Ancho, Alto) en centimetros
    inPerGraph = cmPerGraph/2.54
    numColGraph = 5
    numRowGraph = len(dirs)
    figSize = (inPerGraph[0]*numColGraph, inPerGraph[1]*numRowGraph)
    pdb.set_trace()
    #dirs = [0.1, 0.18, 0.26, 0.34, 0.42, 0.5]

    if os.path.isfile('factores.dat'):
        os.remove('factores.dat')

    fout = open('factores.dat', 'a')

    fa = []
    fev = []
    fes = []
    for i,Alpha in enumerate(dirs):
        os.chdir(str(Alpha))
        Instance = calcE(basedir=rootdir, alpha=Alpha, vmin=-1554.63, freq=1*10**(-3),
                    T=300, biasFile='bias.dat', enerFile='ener.dat',
                    tempFile='temp.dat', AMD=True, tMax=50000, dt=0.2, dFrame=10,
                    xyzFile='traj.xyz', MD_enerFile='DM_ener.dat', dFrameXYZ=50, 
                    generatePlot=True, filePlotName='alphaTest.pdf', multiPlot=True,
                    numColGraph=numColGraph, numRowGraph=numRowGraph, stepId=i)
        #Ecut = Instance.evalEcut()
        #diffE = Instance.diffE()
        #ExVb = Instance.evalExVb()
        #Instance.histo(i)
        #VbProm = Instance.evalVbProm()
        #Instance.plotVOverlaped()
        fa.append(Instance.compareTimes())  # Factor de aceleracion
        fev.append(Instance.escapeFactor()) # Factor de escape verdadero
        fes.append(Instance.freq)           # Factor de escape supuesto
        #fout.write("Alpha = %s\n" % Alpha)
        #fout.write("Factor aceleracion: %.2f\n" % fa)
        #fout.write("Factor escape verdadero: %.2f\n" % fev)
        #fout.write("Factor de escape supuesto: %.3f\n\n" % fes)
        #inertiaDif = Instance.inertia()
        #plt.savefig(str(Alpha)+'.pdf', format='pdf', bbox_inches='tight', dpi=50)
        #plt.close()

        # Factor de escape verdadero
        os.chdir(os.pardir)

    fout.close()
    pdb.set_trace()
    figure1 = plt.figure(title='AlphaTEST')
    figure1.plot(dirs, fa, 'ro-', dirs, fev, 'g^-')

    #for i,Alpha in enumerate(dirs):
        #os.chdir(str(Alpha))
        #Instance = calcE(basedir=rootdir, alpha=Alpha, vmin=-1554.63, freq=1*10**(-3),
                    #T=300, biasFile='bias.dat', enerFile='ener.dat',
                    #tempFile='temp.dat', AMD=True, tMax=50000, dt=0.2, dFrame=10,
                    #xyzFile='traj.xyz', MD_enerFile='DM_ener.dat', dFrameXYZ=50)
        #Instance.histo(i+1)

        ## Factor de escape verdadero
        #os.chdir(os.pardir)

    ##plt.subplot(bottom=0.3)
    #plt.savefig('Histograms.pdf', format='pdf', bbox_inches='tight', dpi=100)
    #plt.close()



    #print 'Energy = %g eV' % (energy)
    #print 'DiffEnergy = %g eV' % (diffener)
    #print '<Vb>b = %g eV' % (VbProm)
    #print 'Expected <Vb>b = %g eV' % (ExVb)
