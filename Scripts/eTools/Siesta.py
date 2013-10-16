#!/usr/bin/env python
#-*- coding: utf-8 -*-

import os
import re
import eTools

class Run(object):

    def __init__(self, sysName=None, SystemLabel='siesta', InputFile='input.fdf',
                 outFile='output.out', final_energy=None):
        self.sysName = sysName
        self.systemLabel = SystemLabel
        self.InputFile = InputFile
        self.xyzFile = SystemLabel+'.xyz'
        self.aniFile = SystemLabel+'.ANI'
        self.outFile = outFile
        self.final_energy = final_energy
        if os.path.isfile(self.xyzFile):
            self.FinalMol = eTools.Molecule(name=sysName)
            self.FinalMol.load_from_file(self.xyzFile)
        self.process_input()
        self.process_output()
        self.frames = self.get_frames()


    def process_input(self):

        total_atoms_p = re.compile(r'NumberOfAtoms\s+(\d+)')

        if os.path.isfile(os.path.join(self.sysName, self.InputFile)):
            with open(os.path.join(self.sysName, self.InputFile)) as f:
                for line in f:
                    total_atoms_m = total_atoms_p.match(line)
                    if total_atoms_m:
                        self.total_atoms = int(total_atoms_m.group(1))


    def process_output(self):
        real_p = r'-?(\d+(\.\d*)?|\d*\.d+)'
        final_ener_p = re.compile(r'.*Total =\s+(' + real_p + r')')
        minim_ener_p = re.compile(r'siesta: E_KS\(eV\) =\s+('+ real_p +r')')
        minim_forc_p = re.compile(r'\s+Max\s+(' + real_p + r')\s+constrained')

        minim_ener = []
        minim_forc = []

        with open(os.path.join(self.sysName, self.outFile), 'r') as f:
            for line in f:
                minim_ener_m = minim_ener_p.match(line)
                minim_forc_m = minim_forc_p.match(line)
                final_ener_m = final_ener_p.match(line)

                if minim_ener_m:
                    ener = float(minim_ener_m.group(1))
                    minim_ener.append(ener)
                elif minim_forc_m:
                    force = float(minim_forc_m.group(1))
                    minim_forc.append(force)
                elif final_ener_m:
                    final_ener = float(final_ener_m.group(1))
                    self.final_energy = final_ener

        self.min_forces = minim_forc
        self.min_energies = minim_ener


    def get_frames(self):

        frames = []
        #a_mol = eTools.Molecule()

        if os.path.isfile(os.path.join(self.sysName, self.aniFile)):
            with open(os.path.join(self.sysName, self.aniFile), 'r') as f:
                lines = f.readlines()
            while len(lines) > 0:
                a_mol = eTools.Molecule()
                #a_mol = eTools.Molecule()
                mol_as_list = lines[2:self.total_atoms+2]
                a_mol.load_from_list(mol_as_list)
                frames.append(a_mol)
                lines = lines[self.total_atoms+2:]

        else:
            a_mol = eTools.Molecule()
            a_mol.load_from_file(os.path.join(self.sysName, self.xyzFile))
            frames.append(a_mol)

        return frames

    def xyz_movie(self, fileName=None):

        if not fileName:
            fileName = self.sysName + '-movie.xyz'

        for mol in self.frames:
            mol.write_to_file(fileName)

    #def job_complete(self):
        #complete_pattern = r'.*End of run.*'


class BigRun(object):

    def __init__(self, sysName=None):
        self.sysName = sysName
        os.chdir(sysName)
        ordered_dirs = eTools.sort_list_of_strings(eTools.get_dirs())
        self.Runs = [Run(sysName=run) for run in ordered_dirs]
        os.chdir(os.pardir)

    def xyz_movie(self, fileName=None):

        if not fileName:
            fileName = self.sysName + '-movie.xyz'

        if os.path.isfile(fileName):
            os.remove(fileName)

        for run in self.Runs:
            run.xyz_movie(fileName)

    def get_runs(self):

        return [float(run.sysName) for run in self.Runs]


    def get_energies(self):

        return [(float(run.sysName), run.final_energy) for run in self.Runs
                     if run.final_energy]

    #def extract_energy(self):
        #"""
        #For a given stretching job, return as a numpy array, the stretching
        #distance, energy_value.
        #"""

        #real_p = r'-?(\d+(\.\d*)?|\d*\.d+)'
        #patt = r'.*Total =\s+(' + real_p + r')'
        #p = re.compile(patt)

        #dist_list = []
        #ener_list = []
        #for d in self.get_dirs():
            #os.chdir(d)
            #if os.path.isfile(self.outFile):
                #with open(self.outFile, 'r') as fopen:
                    #for line in fopen:
                        #match = p.search(line)
                        #if match:
                            #ener = float(match.group(1))
                            #dist_list.append(float(d))
                            #ener_list.append(ener)
            #else:
                #print('WARNING: file %s not found. Continuing with next'
                        #' directory...' % (self.outFile))
            #os.chdir(os.pardir)

        #return np.array(dist_list), np.array(ener_list)



if __name__ == "__main__":
    pass
