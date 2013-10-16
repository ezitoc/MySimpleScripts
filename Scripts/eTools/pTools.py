#!/usr/bin/env python
#-*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import numpy as np
import os
import re


def DiffData(x, y):
    dif_x = x[1:] - x[:-1]
    x_diff = x[:-1] + dif_x/2
    dif_y = y[1:] - y[:-1]
    y_diff = dif_y/dif_x
    return x_diff, y_diff

def FitBlock(data=None, x=None, y=None, sysName=None, poli='Polynomial'):
    """
    Linear regression of square fit of some x,y data.
    type, can be 'Polynomial' or 'Chevbyshev'.
    interv, boolean. Weather it creates a total or partial fit.
    """

    def plot_partial_fit(x_to_fit, y_to_fit, grad=None, poli=poli):
        """
        Add partial fits to axes.
        """
        xfit, yfit, dyfit = FitData(x_to_fit, y_to_fit,
                                    force=True, grad=grad, type=poli)
        ax1.plot(xfit, yfit, color='r', linestyle='-',
                    linewidth='1.5')
        ax2.plot(xfit, dyfit, color='blue', linewidth='1.5')

    if data:
        x = [a_data[0] for a_data in data]
        y = [a_data[1] for a_data in data]
    elif not x and not y:
        print 'Must declare either "x" and "y" or two-element set'

    if not sysName:
        sysName = os.path.basename(os.getcwd())


    real_p = r'-?(\d+(\.\d*)?|\d*\.d+)'

    labels = ['P{0}'.format(i) for i in range(len(x))]
    f0 = plt.figure()
    f0.suptitle(sysName)
    ax0 = f0.add_subplot(111)
    ax0.set_ylabel('y')
    ax0.set_xlabel('x')
    ax0.scatter(x, y, s=50, facecolors='none', edgecolors='red')
    for label, ax, ay in zip(labels, x, y):
        plt.annotate(label, xy=(ax, ay), xytext=(0,20),
                        textcoords='offset points')
    plt.savefig('fit_'+sysName+'.png', format='png', bbox_inches='tight')
    plt.show()

    f2, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=False)
    ax1.set_ylabel('Energy [eV]')
    ax2.set_ylabel('Force [nN]')
    ax2.set_xlabel('Stretching distance [\AA]')
    ax2.axhline(linestyle='--', color='gray', linewidth='2')
    f2.subplots_adjust(hspace=0)
    f2.suptitle(sysName)
    ax1.scatter(x, y, s=30, facecolors='none', edgecolors='blue')

    left_x = x
    left_y = y

    while True:
        print 'Range of point numbers? (comma separated):'
        print 'Enter "r"  to fit all remaining points'
        print 'Enter "s"  to stop fitting'

        sel = raw_input('--> ')
        selist = re.split(',|:', sel)

        if selist[0] == 'r':
            if len(selist) == 2:
                grado = int(selist[1])
            else:
                grado = None
            x_to_fit = []
            y_to_fit = []
            for i in range(len(left_x)):
                if re.match(real_p, str(left_x[i])):
                    x_to_fit.append(left_x[i])
                    y_to_fit.append(left_y[i])
                    if i == len(left_x) - 1 and len(x_to_fit) > 1:
                        plot_partial_fit(x_to_fit, y_to_fit, grad=grado)
                    elif left_x[i+1] == 'x' and len(x_to_fit) > 1:
                        plot_partial_fit(x_to_fit, y_to_fit, grad=grado)
                        x_to_fit = []
                        y_to_fit = []
            break

        elif selist[0] == 's':
            break
        else:
            a = int(selist[0])
            b = int(selist[1])
            if len(selist) == 3:
                grado = int(selist[2])
            else:
                grado = None
            x_to_fit = x[a:b+1]
            y_to_fit = y[a:b+1]
            left_x = ['x' if i in range(a,b+1) else e for i, e in
                      enumerate(left_x)]
            plot_partial_fit(x_to_fit, y_to_fit, grad=grado)

    plt.savefig('energy-force_'+sysName+'.png', format='png',
                bbox_inches='tight')
    plt.show()

def FitData(x, y, type='Polynomial', force=False, grad=None):
    """
    Return x and y values of best fit for x y data set by square fit of desired
    type.

    type: at this moment 'Polynomial' or 'Chebyshev' are the posible options
    """

    x_new = np.linspace(x[0], x[-1], num=len(x)*20)

    if not grad:
        grad = int(len(x)/2)+1
        if grad > 6:
            grad = 6

    if type == 'Polynomial':
        coefs = np.polynomial.polynomial.polyfit(x, y, grad)
        ffit = np.poly1d(coefs[::-1])
        ffitd = np.polyder(ffit)

    elif type == 'Chebyshev':
        ffit = np.polynomial.Chebyshev.fit(x, y, grad)
        ffitd = np.polynomial.Chebyshev.deriv(ffit)

    if force:
        return x_new, ffit(x_new), ffitd(x_new)
    else:
        return x_new, ffit(x_new)


def main():
    pass

if __name__ == "__main__":
    main()
