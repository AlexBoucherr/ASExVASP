from ase.optimize.bfgs import BFGS
from ase.io import read
from ase.calculators.vasp import Vasp2


atoms = read('POSCAR')

dft = Vasp2(atoms=atoms, directory='./', ncore=10, prec='normal', xc='rpbe', lwave=False, lcharg=False, ivdw=11,
            icharg=0, lmaxmix=2, ispin=2, kpts=(1, 1, 1), gamma=True, charge=0, encut=600, ediff=1E-6, lreal='auto',
            nelm=700, algo='conjugate', lasph=True, laechg=False, lelf=False, ldipol=False, idipol=3, ismear=0,
            sigma=0.1, ediffg=-0.01, ibrion=1, nsw=500, potim=0.4)

atoms.calc = dft
atoms.get_potential_energy()
