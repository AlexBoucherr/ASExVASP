########################################################################################################################
# Performs a classic calculation. Change calculator settings if needed.                                                #
########################################################################################################################
from ase.build import bulk
from ase.calculators.vasp import Vasp2
from ase.io import read

fermi = True  # Get fermi energy, True/False.

'''
Set up system and parameter DFT. Bulk Pd is test system. Probably better to read input structure from a file or just 
change the script.
'''
# The system.
atoms = read('./POSCAR')

# DFT parameters
dft = Vasp2(atoms=atoms,      # working structure.
            directory='./',   # Working dir.
            ncore=5,          # number of cores.
            prec='accurate',  # Default accuracy.
            nupdown=0,        # Force electrons to be paired. Default: no set.
            xc='rpbe',        # xc functionnal.
            lwave=True,
            lcharg=True,
            ivdw=11,          # VdW correction. 11 = DFT-D3.
            icharg=0,         # Initial guess on charge density. Default: 0.
            # Electronic relaxation.
            lmaxmix=4,        # Used d-orbitals in initial guess for WF. Default: 2 (p-orbitals).
            ispin=2,          # Spin-polarized calculation.
            kspacing=0.15,    # Min distance (rec. space) between 2 k-points.
            gamma=True,       # Gamma-centred mesh.
            charge=0,         # Charge of the system.
            encut=500,        # Kinetic cut-off energy.
            ediff=1E-5,       # Electronic relaxation cut-off.
            lreal='auto',
            nelm=500,         # Max iteration electronic relaxation.
            algo='normal',    # algo used for electronic relaxation.
            lasph=True,
            laechg=False,
            lelf=False,
            ldipol=False,     # Dipole correction.
            idipol=3,         # 1: x-axis, 2: y-axis, 3: z-axis, 4: all axis.
            # Density of states.
            ismear=2,         # Smearing method. n>0 = Methfessel-Paxton order n.
            sigma=0.2,        # Smearing parameter (eV).
            # Ionic relaxation.
            ediffg=-0.01,     # Energy threshold for ionic relaxation (IR).
            ibrion=2,         # IR algo. 2 = CG algo, 1 = QN algo.
            nsw=100,          # Max iterations for IR.
            potim=0.2,
            smass=0.5)

'''
Performs calculation.
'''
atoms.calc = dft
atoms.get_potential_energy()

if fermi:
    # Get Fermi level of the system in a file called FERMI.
    E_Fermi = dft.get_fermi_level()
    EFERMI = open('FERMI', 'w')
    EFERMI.write('Fermi level eV): ' + '\n' + str(round(E_Fermi, 4)) + '\n')
    EFERMI.close()
