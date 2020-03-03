
import os
import shutil

from ase.lattice.hexagonal import *
from ase.calculators.dftb import Dftb
from ase.io import write

def run_dftb(system, c, records):
	calc = Dftb(label='graphene', atoms=graphene,kpts=(1,1,1),
                    run_manyDftb_steps=False,
                    # Specific DFTB+ keywords
                    Hamiltonian_SCC='Yes',
                    Hamiltonian_SCCTolerance=1E-5,
                    Hamiltonian_Charge = -1.0,
                    Hamiltonian_Filling_ = 'Fermi',
                    Hamiltonian_Filling_Temperature='1E-6',
                    Hamiltonian_SpinPolarisation_ = 'Colinear',
                    Hamiltonian_SpinPolarisation_UnpairedElectrons = 1.0,
                    Hamiltonian_SpinConstants_ = '',
                    Hamiltonian_SpinConstants_C = '{ -0.023}',
                    Hamiltonian_ElectricField_ = '',
                    Hamiltonian_ElectricField_PointCharges_ = '',
                    Hamiltonian_ElectricField_PointCharges_CoordsAndCharges_ = 'DirectRead',
                    Hamiltonian_ElectricField_PointCharges_CoordsAndCharges_Records = records,
                    Hamiltonian_ElectricField_PointCharges_CoordsAndCharges_File = 'charges.dat'
                   )


	graphene.set_calculator(calc)
	calc.calculate(system)

	return()

def clean_up_folder(resultsdir, resultsfile, tidyfile):
	shutil.move(resultsfile,os.path.join(resultsdir, tidyfile+'_'+str(int(c)))) 

	return()

# Set up an infinite graphene sheet
simulation_name = 'graphene_9x9_1e_spinp_ft_1e-6_pointcharges'
#simulation_name = 'graphene_9x9_ft_1e-5'
dftb_files = ['detailed.out', 'band.out', 'dftb_in.hsd', 'graphene.out']
n    = 9
alat = 2.472
clat = [10., 20., 30., 40., 50., 60., 70., 80., 90., 100.]

# Make a dedicated folder for the simulation results
resultsdir = os.path.join(os.getcwd(),simulation_name)
if not os.path.exists(resultsdir):
	os.mkdir(resultsdir)

# Loop over the layer spacing
for c in clat:
	graphene = Graphene(symbol='C', latticeconstant={'a':alat, 'c':c}, size=(n,n, 1))

	# Figure out how many charges in the cell
	charges = open('charges.dat', 'r').readlines()
	records = 0
	for line in charges:
		if float(line.split()[2]) < (c - 0.5):
			records += 1

	run_dftb(graphene, int(c), records)

	# Move the results to the new folder
	for f in dftb_files:
		output=os.path.join(os.getcwd(),f)
		clean_up_folder(resultsdir, output, f)
