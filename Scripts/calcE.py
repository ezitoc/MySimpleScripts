#!/usr/bin/env python

#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
#
#* File Name : calcE.py
#
#* Purpose :
#
#* Creation Date : 11-03-2013
#
#* Last Modified : Mon 08 Apr 2013 06:20:03 PM ART
#
#* Created By :  Ezequiel Castillo
#
#_._._._._._._._._._._._._._._._._._._._._.

import sys
import os
import math
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
from scipy.optimize import curve_fit
import pdb

rc('font',**{'family':'serif','serif':['Palatino'],'size':8})
rc('text', usetex=True)

dictAtomMass = \
{
'Au':196.967,
'C':12.011,
'H':1.008,
'S':32.066}

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
    fin = open(infile, "r")
    str1 = fin.read()
    fin.close()
    # add 'replace pattern' : 'replacement' to dictionary
    dict_make = dictFromSequence(replace_pattern)
    # call the function and get the changed text
    str2 = replace_words(str1, dict_make)
    # write changed text back out
    fout = open("run-"+infile, "w")
    fout.write(str2)
    fout.close()


def checkCorr(Xvalues, Yvalues, showPlot=False):

    def func(t, a, b, c, d, e):
    #def func(t, a, b, c):
        """Fitting exponential function"""
        #return a*np.exp(-t/b)+np.exp(-t/c)
        return np.exp(-t/a)+np.exp(-t/b)+np.exp(-t/c)+np.exp(-t/d)+np.exp(-t/e)

    # Ignore overflow error (for curve_fit)
    np.seterr(over='ignore')

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

class createInputs(object):

    def __init__(self, projectName='project', modFiles=None, filesToCopy=[]):
        self.projectName = projectName
        self.rootDir = os.path.abspath('.')
        self.ignoreFiles = modFiles
        self.filesToCopy = filesToCopy

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
            return filesListMod
        else:
            return filesList

    def createFolders(self, suffixNo=None, suffixName=None, foldersName="run"):
        """Create subfolders containing all files declared in the filesToCopy
        option (all files included, otherwise) and ignoring those declared at the
        modFiles option"""

        if suffixNo and suffixName:
            print "ERROR: Two suffix types declared. Declare only one."
            sys.exit(1)

        if suffixNo or suffixName:
            if suffixNo:
                self.suffixList = range(suffixNo)
            else:
                self.suffixList = suffixName
        else:
            print "ERROR: No suffix declared."
            sys.exit(1)

        self.foldersName = foldersName
        # CWD = CWD/projectName
        if not os.path.isdir(self.projectName):
            os.mkdir(self.projectName)
        else:
            print 'WARN: folder "%s" already exists.' % (self.projectName)
        os.chdir(self.projectName)
        for suffix in self.suffixList:
            self.folderName = '%s_%s' % (self.foldersName, suffix)
        # CWD = CWD/projectName/folderName
            os.mkdir(self.folderName)
            os.chdir(self.folderName)
            copySelectedFiles = self.selectFiles()
            for file in copySelectedFiles:
                shutil.copy(os.path.join(self.rootDir, file), os.curdir)
        # CWD = CWD/projectName
            os.chdir(os.pardir)

        os.chdir(self.rootDir)

class createProject(createInputs):
    """  Ready to create inputs """

    def __init__(self, templateFile=None, replacePatternAsList=None):
        self.templateFile = templateFile
        self.replacePatternAsList = replacePatternAsList

    def modify(self):
        multiple_subst(self.templateFile, self.replacePatternAsList)




class calcE(object):
    """ Super class"""

    def __init__(self, vprom=None, alpha=None, vmin=None, Ecut=None, freq=None,
                 T=None, biasFile=None, enerFile=None, tempFile=None,
                 smoothFactor=None, tau=None, AMD=False, tMax=None, dt=None,
                 dFrame=None, xyzFile=None, MD_enerFile=None, basedir=None,
                 dFrameXYZ=None):

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


        self.kB = 8.6173324*10**-5 # eV/K
        self.T = T # K
        self.beta = 1./(self.kB*self.T)
        self.alpha = alpha
        self.vmin = vmin
        self.freq = freq
        self.biasFile = biasFile
        self.VSinBias = np.loadtxt(self.biasFile, usecols=[0])
        self.deltaBias = np.loadtxt(self.biasFile, usecols=[1])
        self.VBias = np.loadtxt(self.biasFile, usecols=[2])
        self.tempFile = tempFile
        self.enerFile = enerFile
        self.enerPot = np.loadtxt(self.enerFile, usecols=[0])
        self.enerTot = np.loadtxt(self.enerFile, usecols=[2])
        self.timeBoost = np.loadtxt(self.tempFile, usecols=[0])
        self.smoothFactor = smoothFactor
        self.tau = tau
        self.xyzFile = xyzFile
        self.dtStepsXYZ = dt * dFrameXYZ
        self.timeXYZ = np.arange(self.dtStepsXYZ, tMax+self.dtStepsXYZ, self.dtStepsXYZ)

        if vprom:
            self.vprom = vprom
        else:
            self.vprom = np.average(self.MD_enerPot)

    def histo(self):
        plt.figure(1)
        plt.subplot(221)
        plt.hist(self.MD_enerPot, 50, color='g')
        plt.axvline(x=self.evalEcut(), linewidth=2, color='r')
        #plt.show()

    def evalEcut(self):
        self.Ecut = (self.vprom + self.vmin * (self.alpha - 1) - self.kB *\
                    self.T * math.log(self.freq)) / self.alpha
        return self.Ecut

    def diffE(self):
        self.diffE = self.Ecut - self.vmin
        return self.diffE

    def evalExVb(self):
        self.ExVb = self.vprom + (self.vmin - self.Ecut) * (self.alpha - 1)
        return self.ExVb

    def evalVbProm(self):
        self.VbProm = np.average(self.VBias)
        return self.VbProm

    def plotVOverlaped(self, plot=True):
        """ We want to overlap the V curve of a normal MD with the acelerated one """

        # Arreglo con VbSN[i]=V[i]*exp(beta*DeltaV[i]), i.e., el
        # numerador del valor de expectacion
        VbSN = np.array([ self.VSinBias[i]*math.exp(self.beta*self.deltaBias[i]) for i in range(len(self.VSinBias)) ])

        # Arreglo con normFac[i]=exp(beta*DeltaV[i]), i.e., el
        # denominador del valor de expectacion
        normFac = np.array( [math.exp(self.beta*self.deltaBias[i]) for i in range(len(self.deltaBias))]  )

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
            value = value + decor*math.exp(self.beta*np.average(block))
            thyper.append(value)

        if plot:
            plt.subplot(224)
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
            plt.plot(self.timeNoBoost, self.timeBoost, 'o-g', label='alpha='+str(self.alpha), markersize=1)
            plt.plot(self.timeNoBoost, a*self.timeNoBoost + b, 'r')
            plt.legend(loc='upper left')
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


        if plot:
            plt.subplot(223)
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

    # The current working directory
    #thedir = os.path.abspath(os.curdir)
    thedir = os.getcwd()

    dirs = [ float(name) for name in os.listdir(thedir) if os.path.isdir(os.path.join(thedir, name)) ]

    pdb.set_trace()



    #for Alpha in dirs:
        #os.chdir(os.path.join(thedir, str(Alpha)))
        #Instance = calcE(basedir=thedir, alpha=Alpha, vmin=-1554.63, freq=1*10**(-3),
                    #T=300, biasFile='bias.dat', enerFile='ener.dat',
                    #tempFile='temp.dat', AMD=True, tMax=50000, dt=0.2, dFrame=10,
                    #xyzFile='traj.xyz', MD_enerFile='DM_ener.dat', dFrameXYZ=50)
        #Ecut = Instance.evalEcut()
        #diffE = Instance.diffE()
        #ExVb = Instance.evalExVb()
        #VbProm = Instance.evalVbProm()
        #Instance.plotVOverlaped()
        #fa = Instance.compareTimes()  # Factor de aceleracion
        #fev = Instance.escapeFactor() # Factor de escape verdadero
        #fes = Instance.freq           # Factor de escape supuesto
        #inertiaDif = Instance.inertia()
        #Instance.histo()
        #plt.savefig(str(Alpha)+'.pdf', format='pdf', bbox_inches='tight')
        #plt.close()

        ## Factor de escape verdadero
        #os.chdir(thedir)

    #print 'Energy = %g eV' % (energy)
    #print 'DiffEnergy = %g eV' % (diffener)
    #print '<Vb>b = %g eV' % (VbProm)
    #print 'Expected <Vb>b = %g eV' % (ExVb)
