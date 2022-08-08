"""
Use this script to generate input files in the form of .traj file for NEB calculations using VASP.
First a linear extrapolatio is performed. If man_check (manual check, l.83) is enabled, a loop runs over each structure
in the linear extrapolation asking the used to move an atom if needed.

To leave the image as it is, answer none to the line asking 'Select the atom index you need to move:'.

Alex, v.1.0, 03/08/2022.
"""

from ase.io import read, Trajectory
from ase.visualize import view
from ase.constraints import FixAtoms
import numpy as np
import os


# Following function returns the initial, non-corrected path.
def get_neb_extrapolation(initialPosition, n, displacedIndex, finalPositions):
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


# The following function moves atoms in an image in the initial .traj file.
def move_atom(file):
    modif = 'yes'
    view(file)
    while modif == 'yes':
        atomToMove = input('Select the atom index you need to move:\n')
        if atomToMove == 'none':
            return file
        else:
            atomToMove = atomToMove.split()
            print('Select the delta you need to mode the atom along x, y and z-axis (A).')
            deltaX, deltaY, deltaZ = float(input('delta X: ')), float(input('delta Y: ')), float(input('delta Z: '))

            for element in atomToMove:
                element = int(element)
                for atom in file:
                    if atom.index == element:
                        atom.position[0] += deltaX
                        atom.position[1] += deltaY
                        atom.position[2] += deltaZ

            view(file)
            accept = input('Accept changes? (yes/no)\n')
            if accept == 'no':
                for element in atomToMove:
                    element = int(element)
                    for atom in file:
                        if atom.index == element:
                            atom.position[0] -= deltaX
                            atom.position[1] -= deltaY
                            atom.position[2] -= deltaZ
            modif = input('Perform further modifications? (yes/no)\n')
    return file


initial = read('/home/alex/Downloads/coordination_sites/hollow_site/frequencies/CONTCAR')
final = None
n_images = 5
displaced_atoms = [96]
displacement_vector = [[-2.4, -4.3, 0.0]]
man_check = True

"""
If needed, set up new constraints to your system.
"""
frozenAtoms = []
initial.constraints = False
for atom in initial:
    if round(atom.position[2], 3) < 13.870:
        frozenAtoms.append(atom.index)

initial.set_constraint(FixAtoms(frozenAtoms))

# First situation; Final = None, i.e., a final position is not provided.
if final is None:
    """
    Generates a first trajectory file containing images of the TS path.
    """
    get_neb_extrapolation(initial, n=n_images, displacedIndex=displaced_atoms, finalPositions=displacement_vector)

    # Reads the trajectory file created by get_neb_extrapolation.
    trajectory = Trajectory('neb.traj', 'r')

    """
    If man_check (manual check), then view each image one by one and adjust positions if needed.
    If you don't want to move an atom, answer none when asked an atom index.
    """
    if man_check:
        os.remove('./neb.traj')
        trajectory_final = Trajectory('neb_final.traj', 'w')
        for i in range(0, len(trajectory)):
            atoms = move_atom(file=trajectory[i])
            trajectory_final.write(atoms)

    trajectory = Trajectory('neb_final.traj', 'r')
    print(trajectory[1])
    view(trajectory[1])

# If a final position is provided:
if final is not None:
    displaced_atoms = []
    displacement_vector = []
    for i in range(0, len(initial)):
        delta = final[i].position - initial[i].position
        if round(np.linalg.norm(delta), 3) > 0.000:
            displaced_atoms.append(i)
            displacement_vector.append([round(delta[0], 3), round(delta[1], 3), round(delta[2], 3)])

