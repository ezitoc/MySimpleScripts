#!/usr/bin/env python

#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
#
#* File Name : pyGeMS.py
#
#* Purpose :
#
#* Creation Date : 22-05-2013
#
#* Last Modified : Wed 22 May 2013 07:38:56 PM ART
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

    def first_run():
        FOLDERNAME = 'firstRun'
        FILES = 'gems.py initial.xyz'

        # Chequeamos directorio
        if not os.path.isdir(FOLDERNAME):
            os.mkdir(FOLDERNAME)

        # Nos vamos al directorio creado
        os.chdir(FOLDERNAME)

        # Copiamos archivos pertinentes

