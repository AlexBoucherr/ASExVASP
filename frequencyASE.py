import os
from ase.calculators.vasp import Vasp2
from ase.io import read, Trajectory
from ase.constraints import FixAtoms
import numpy as np

traj = Trajectory('nebOutput.traj', 'r')
atoms = traj[1]

dft = Vasp2(atoms=atoms, ncore=15, directory='./', prec='normal', xc='rpbe', ivdw=11, icharg=0, lmaxmix=2, ispin=2,
            kspacing=0.2, gamma=True, charge=0, encut=600, ediff=1E-5, lreal='auto', nelm=500, algo='fast', lasph=True,
            laechg=False, lelf=False, ldipol=False, idipol=3, ismear=0, sigma=0.1, ediffg=1e-3, ibrion=-1, nsw=0,
            potim=0.2, lwave=True, lcharg=False)

# First calculation.
atoms.calc = dft
atoms.get_potential_energy()

'''
Now, set-up the frequency calculation.
'''

# first freeze the max. atoms possible to make calculation lighter.
atoms = read('./CONTCAR')
atoms.constraint = False
frozenAtoms = []
for atom in atoms:
    if atom.index != 96:
        frozenAtoms.append(atom.index)

atoms.set_constraint(FixAtoms(indices=frozenAtoms))

# Reset the first DFT calculation and set-up frequency calculation.
dft.reset()
dft.set(atoms=atoms,
        directory='./',
        ncore=1,      # /!\ TEST FOR TS RESEARCH.
        ediff=1E-7,
        ediffg=1E-7,
        ibrion=5,     # Request frequency calculation.
        nsw=500,      # 500 iteration max.
        potim=0.03)   # really small POTIM recommended for freq. research.

atoms.calc = dft
atoms.get_potential_energy()
