from ase.build import bulk
from ase.calculators.vasp import Vasp2
import matplotlib.pyplot as plt
from ase.dft import DOS
import numpy as np

computeDOS = True  # Compute DOS or no.
splitSpin = True  # Split spin up and spin down.

'''
Set up system and parameter DFT.
'''
# The system.
bulk = bulk('Pd', a=3.929)

# DFT parameters
dft = Vasp2(atoms=bulk,  # working structure.
            directory='./',  # Working dir.
            ncore=5,  # number of cores.
            prec='accurate',  # Default accuracy.
            nupdown=0,  # Force electrons to be paired.
            xc='rpbe',  # xc functionnal.
            ivdw=11,  # VdW correction. 11 = DFT-D3.
            icharg=0,
            lmaxmix=4,  # Used d-orbitals in initial guess for WF.
            ispin=2,  # Spin-polarized calculation.
            kspacing=0.15,
            gamma=True,  # Gamma-centred mesh.
            charge=0,  # Charge of the system.
            encut=500,  # Kinetic cut-off energy.
            ediff=1E-5,  # Electronic relaxation cut-off.
            lreal='auto',
            nelm=500,  # Max iteration electronic relaxation.
            algo='normal',  # algo used for electronic relaxation.
            lasph=True,
            laechg=False,
            lelf=False,
            ismear=2,  # Smearing method. n>0 = Methfessel-Paxton order n.
            sigma=0.2,  # Smearing parameter (eV).
            ediffg=-0.01,  # Energy threshold for ionic relaxation (IR).
            ibrion=2,  # IR algo. 2 = CG algo.
            nsw=100,  # Max iterations for IR.
            potim=0.2,
            smass=0.5,
            lwave=True,
            lcharg=True)

'''
Performs classic energy calculation.
'''
bulk.calc = dft
bulk.get_potential_energy()

# Get Fermi level of the system in a file called FERMI.
E_Fermi = dft.get_fermi_level()
EFERMI = open('FERMI', 'w')
EFERMI.write('Fermi level eV): ' + '\n' + str(round(E_Fermi, 4)) + '\n')
EFERMI.close()

'''
Set up Desity of States (DOS) calculation.
'''
if computeDOS:
    if not splitSpin:
        # parameters for DOS calcualtion.
        energies, weights = dft.get_dos(spin=None,  # Both spin states accounted for (0, spin up, 1 spin down).
                                        npts=1001,  # Number of points.
                                        width=0.2)  # Smearing parameter (eV).

        # Write data in output file.
        EnergyFile = open('DOSenergies.txt', 'w')
        for element in energies:
            EnergyFile.write(str(element) + '\n')
        EnergyFile.close()
        DosFile = open('DOSweights.txt', 'w')
        for element in weights:
            DosFile.write(str(element) + '\n')
        DosFile.close()

        # Plot and save DOS figure.
        plt.plot(energies, weights, color='#CC0000', linewidth=0.7)
        plt.xlabel('Energy, eV')
        plt.ylabel('Density of State (DOS)')
        plt.xlim(min(energies), max(energies))
        plt.ylim(0, max(weights) + 0.1)
        plt.savefig('DOS.png')
        plt.show()
    if splitSpin:
        dos = DOS(dft, npts=1001, width=0.2)
        # Calculation spin up:
        dUP = dos.get_dos(spin=0)  # spin=0; spin up.
        eUP = dos.get_energies()

        # Calculation spin down:
        dDOWN = dos.get_dos(spin=1)
        eDOWN = dos.get_energies()

        energyUP = open('ENERGUP.txt', 'w')
        energyDOWN = open('ENERGDOWN.txt', 'w')
        dosUP = open('DOSUP.txt', 'w')
        dosDOWN = open('DOSDOWN.txt', 'w')

        for element in dUP:
            dosUP.write(str(element) + '\n')
        for element in dDOWN:
            dosDOWN.write(str(element) + '\n')
        for element in eUP:
            energyUP.write(str(element) + '\n')
        for element in eDOWN:
            energyDOWN.write(str(element) + '\n')

        energyUP.close()
        energyDOWN.close()
        dosUP.close()
        dosDOWN.close()

        filePath = './'
        dosUPpath = filePath + 'DOSUP.txt'
        dosDOWNpath = filePath + 'DOSDOWN.txt'
        energyPath = filePath + 'ENERGUP.txt'


        def ExtractData(path):
            list = []
            # Read file.
            with open(path) as f:
                line = f.readline()
                while line:
                    line = f.readline()
                    list.append(line)

            # Remove last element which is ''.
            for element in list:
                if element == '':
                    list.remove(element)

            # Remove the '\n' from all elements except the last one, turn str into float.
            for i in range(0, len(list)):  # The last element doesn't have the '\n'. So len - 1.
                list[i] = list[i].strip('\n')
                list[i] = float(list[i])
            return list


        # The data
        energies = ExtractData(path=energyPath)
        dosUP = ExtractData(path=dosUPpath)
        dosDOWN = ExtractData(path=dosDOWNpath)
        dosDOWNneg = [(-1 * i) for i in dosDOWN]  # Let's get the downspin on the bottom of the figure.

        # Equation for axis y = 0.
        X = np.linspace(min(energies), max(energies), 1000)
        Y = [0 for i in X]

        plt.plot(energies, dosUP, color='r', linewidth=0.7, label='DOS, spin up.')
        plt.plot(energies, dosDOWNneg, color='g', linewidth=0.7, label='DOS, spin down.')
        plt.plot(X, Y, '-', color='k', linewidth=0.7)
        plt.axvline(0, color='k', linestyle='-', linewidth=0.8, label='E$_F$$_e$$_r$$_m$$_i$')
        plt.legend()
        plt.xlim(min(energies), max(energies))
        plt.ylim(min(dosDOWNneg) - 0.1, max(dosUP) + 0.1)
        plt.xlabel('E - E$_F$$_e$$_r$$_m$$_i$, eV')
        plt.ylabel('Density of States')
        plt.savefig('DOS.png')
        plt.show()
