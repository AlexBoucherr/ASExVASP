from ase.constraints import FixAtoms
from ase.calculators.vasp import Vasp2
from ase.io import Trajectory, write, read
from ase.neb import NEB
from ase.optimize import BFGS
import numpy as np
from ase.neb import NEBTools

"""
Set necessary parameters in the dft calculator l. 76-110.
If required, change / delete the section setting the constraints.
n_images = number of images requested during the extrapolation process. Extrapolation is performed linearly in this
version of the script. Next version will include verification of the images generated.
atoms_displaced = index of the atoms to displace; [atom1, atom2, ..., atomN].
displacement = displacement vector associated with each displaced atom; [[dx1, dy1, dz1], ..., [dxN, dyN, dzN]].
climb = True/False; enables c-NEB or disables it.

Optimizer is BFGS, uses a force convergence as thresehold to 0.01 eV/A. tuneable, l. 121.
"""

n_images = 3
atoms_displaced = [8, 97]
displacement = [[0.000, 0.000, 0.000], [0.000, 0.000, 0.000]]
climb = True

# Read initial positions.
atoms_initial = read('POSCAR_init')

# Change constraints.
frozenAtoms = []
atoms_initial.constraints = False
for atom in atoms_initial:
    if atom.position[2] < 14.200:
        frozenAtoms.append(atom.index)

atoms_initial.set_constraint(FixAtoms(frozenAtoms))


# Create a function to linearly extrapolate reactant to product.
def NEBextrapolate(initialPosition, n, displacedIndex, finalPositions):
    displacement = [np.array(i) / (n - 1) for i in finalPositions]
    trajectory = Trajectory('neb.traj', 'w', atoms=initialPosition)
    motions = []

    for i in range(0, len(displacement)):
        motions.append([displacedIndex[i], displacement[i]])

    structures = initialPosition
    trajectory.write(atoms=structures)
    for i in range(1, n):
        atoms = structures
        for delta in motions:
            for atom in atoms:
                if atom.index == delta[0]:
                    atom.position[0] = atom.position[0] + (delta[1][0])
                    atom.position[1] = atom.position[1] + (delta[1][1])
                    atom.position[2] = atom.position[2] + (delta[1][2])
        trajectory.write(atoms=atoms)
        structures = atoms

    trajectory.close()
    return


# Create a trajectory file containing images through our path.
NEBextrapolate(initialPosition=atoms_initial,
               n=n_images,
               displacedIndex=atoms_displaced,
               finalPositions=displacement)

# Reads the trajectory file created by NEBextrapolate.
trajectory = Trajectory('neb.traj', 'r')

# Build the serie of images, and asosciate them with a calculator.
images = [atoms for atoms in trajectory]
for image in images:
    dft = Vasp2(atoms=image,  # working structure.
                directory='./',  # Working dir.
                ncore=10,  # number of cores.
                prec='normal',  # Default accuracy.
                xc='rpbe',  # xc functionnal.
                lwave=True,
                lcharg=True,
                ivdw=11,  # VdW correction. 11 = DFT-D3.
                icharg=0,  # Initial guess on charge density. Default: 0.
                # Electronic relaxation.
                lmaxmix=4,  # Used d-orbitals in initial guess for WF. Default: 2 (p-orbitals).
                nupdown=0,
                ispin=2,  # Spin-polarized calculation.
                kspacing=0.2,  # Min distance (rec. space) between 2 k-points.
                gamma=True,  # Gamma-centred mesh.
                charge=0,  # Charge of the system.
                encut=600,  # Kinetic cut-off energy.
                ediff=1E-5,  # Electronic relaxation cut-off.
                lreal='auto',
                nelm=500,  # Max iteration electronic relaxation.
                algo='conjugate',  # algo used for electronic relaxation.
                lasph=True,
                laechg=False,
                lelf=False,
                ldipol=False,  # Dipole correction.
                idipol=3,  # 1: x-axis, 2: y-axis, 3: z-axis, 4: all axis.
                # Density of states.
                ismear=2,  # Smearing method. n>0 = Methfessel-Paxton order n.
                sigma=0.2,  # Smearing parameter (eV).
                # Ionic relaxation.
                ediffg=-0.01,  # Energy threshold for ionic relaxation (IR).
                ibrion=1,  # IR algo. 2 = CG algo, 1 = QN algo.
                nsw=100,  # Max iterations for IR.
                potim=0.4)
    image.calc = dft

# Set-up cNEB.
neb = NEB(images=images, climb=climb)

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
