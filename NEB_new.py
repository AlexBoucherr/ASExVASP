"""
Runs NEB calculation on a trajectory .traj file build by NEB_input_gen.py script.
"""

from ase.calculators.vasp import Vasp2
from ase.io import Trajectory, write
from ase.neb import NEB
from ase.optimize import BFGS
import numpy as np
from ase.neb import NEBTools

# Read .traj file.
neb_traj = Trajectory('neb_input.traj', 'r')

# Build the serie of images, and asosciate them with a calculator.
images = [atoms for atoms in neb_traj]

for image in images:
    dft = Vasp2(atoms=image, directory='./', ncore=10, prec='normal', xc='rpbe', lwave=True, lcharg=True,
                ivdw=11, icharg=0, lmaxmix=4, nupdown=0, ispin=2, kspacing=0.2, gamma=True, charge=0, encut=600,
                ediff=1E-5, lreal='auto', nelm=500, algo='conjugate', lasph=True, laechg=False, lelf=True,
                ldipol=False, idipol=3, ismear=2, sigma=0.2, ediffg=-0.01, ibrion=-1, nsw=0, potim=0.4)
    image.calc = dft

# Set-up cNEB.
neb = NEB(images=images, climb=True)

# Linear interpolation between the images in .traj file.
neb.interpolate()

# Associate an optimizer to the system (BFGS or FIRE)
optimizer = BFGS(neb, trajectory='nebOutput.traj', logfile='optimizer.txt')
optimizer.run(fmax=0.01)

# Get the output trajectory.
write('nebOutput.traj', images)

barrier = NEBTools.get_barrier(neb)
fmax = NEBTools.get_fmax(neb)

outputNEB = open('NEBOUT', 'w')
outputNEB.write('# Energy barrier (eV) of the transformation and dE as TS:\n' + str(barrier) + '\n')
outputNEB.write('\n# Max forces applied on the system:\n' + str(round(fmax, 4)) +
                ' eV/A. Threshold: ' + str(optimizer.fmax) + ' eV/A.\n')

for i in range(0, len(images)):
    outputNEB.write(' \n# Image no. ' + str(i + 1) + '\n')
    outputNEB.write('Total energy: ' + str(round(images[i].get_potential_energy(), 4)) + ' eV\n')
    forces = round(np.linalg.norm(images[i].get_forces()[6]), 4)
    outputNEB.write('Total forces norm: ' + str(forces) + ' eV/A-1\n')

outputNEB.close()
