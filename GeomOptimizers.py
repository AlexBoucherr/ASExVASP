########################################################################################################################
# This script contains a simple way to run a calculation using different optimising algorithm implemented in Python,   #
# (DIIS, CG, DMD) or ASE (FIRE, BFGS, LBFGS).                                                                          #
# Version 1.0 created on 11/11/2021 by Alex. Current version: 1.0, edited: 11/11/2021.                                 #
########################################################################################################################
"""
Set the variable algo to the algorithm needed. If DMD, CG or DIIS, set manually ibrion, nsw, potim...
In the DFT calculator. Some parameters (addgrid, isym, nsim, ldipol, nfree) are set to run calculations
on Pd nanoparticles so far. Proceed with caution working with other systems, and possibly disable these parameters.
When using BFGS or FIRE, set ibrion to -1, and nsw to 0.
"""
from ase.calculators.vasp import Vasp2
from ase.optimize import BFGS, FIRE
from ase.io import read
import ase.io.vasp
from datetime import datetime

algo = 'DIIS'  # BFGS, FIRE (from ASE) or DIIS, DMD or CG (from VASP).

atoms = read('./POSCAR')

# DFT parameters
dft = Vasp2(directory='./',  # Working dir.
            ncore=5,  # number of cores.
            prec='normal',  # Default accuracy.
            xc='rpbe',  # xc functionnal.
            lwave=True,
            lcharg=True,
            ivdw=11,  # VdW correction. 11 = DFT-D3.
            icharg=0,  # Initial guess on charge density. Default: 0.
            # Electronic relaxation.
            lmaxmix=4,  # Used d-orbitals in initial guess for WF. Default: 2 (p-orbitals).
            ispin=2,  # Spin-polarized calculation.
            kspacing=0.15,  # Min distance (rec. space) between 2 k-points.
            gamma=True,  # Gamma-centred mesh.
            charge=0,  # Charge of the system.
            encut=500,  # Kinetic cut-off energy.
            ediff=1E-5,  # Electronic relaxation cut-off.
            lreal='auto',
            nelm=500,  # Max iteration electronic relaxation.
            algo='conjugate',  # algo used for electronic relaxation.
            lasph=True,
            laechg=False,
            lelf=False,
            ldipol=True,  # Dipole correction.
            idipol=4,  # 1: x-axis, 2: y-axis, 3: z-axis, 4: all axis.
            # Density of states.
            ismear=2,  # Smearing method. n>0 = Methfessel-Paxton order n.
            sigma=0.2,  # Smearing parameter (eV).
            # Ionic relaxation.
            ediffg=-0.01,  # Energy threshold for ionic relaxation (IR).
            ibrion=1,  # IR algo. 2 = CG algo, 1 = QN algo.
            nsw=500,  # Max iterations for IR.
            potim=0.8,  # Default: 0.4
            # Optimised parameters.
            nfree=0,
            isym=2,
            addgrid=False,
            nsim=10)

atoms.calc = dft

if algo == 'BFGS':
    dynamics = BFGS(atoms,
                    trajectory='TRAJECTORY.txt',  # Create a file containing trajectories of each atom.
                    logfile='logfile',  # logfile containing data on successive iterations.
                    alpha=58)  # Initial guess for Hessian matrix. Default: 70.

    startTime = datetime.now()  # Starting time of the calculation.
    dynamics.run(fmax=0.01)  # fmax equivalent to ediffg in DFT.

    # Write a file containing time required to run a calculation.
    deltaTime = datetime.now() - startTime  # Time needed to perform the calculation.
    time = open('TIME', 'w')
    time.write('Time required to run calculation:\n' + str(deltaTime))
    time.close()

    # Write a contcar file containing final position.
    ase.io.vasp.write_vasp('CONTCAR2',  # The name of the file we create.
                           atoms,  # The object we use to create the file (our surface).
                           direct=True,
                           vasp5=True,  # VASP5 compatible file.
                           sort=True,
                           ignore_constraints=True)  # If True, the 'FixAtoms set-up' will be ignored.

classicAlgo = ['DIIS', 'DMD', 'CJ']
if algo in classicAlgo:
    atoms.get_potential_energy()
