import os, sys, time
from pathlib import Path
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
cwd = os.getcwd()
cosmo_dir = cwd[: cwd.find('simulation_analysis')] + 'simulation_analysis/'
tools_dir = cosmo_dir + 'tools'
sys.path.append( tools_dir )
from tools import *
#Append analysis directories to path
extend_path()
from constants_cosmo import G_COSMO

use_mpi = False
if use_mpi:
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  n_procs = comm.Get_size()
else:
  rank = 0
  n_procs = 1

# data_dir = '/home/bruno/Desktop/data/'
# data_dir = '/home/bruno/Desktop/ssd_0/data/'
data_dir = '/raid/bruno/data/'
input_dir = data_dir + 'cosmo_sims/1024_hydro_50Mpc/output_files_pchw18/'
output_dir = input_dir + 'neutral_fraction/'
if rank == 0: create_directory( output_dir )

H0 = 67.66
Omega_b =  0.0497 
h = H0 / 100

rho_crit =  3*(H0*1e-3)**2/(8*np.pi*G_COSMO)/ h**2
rho_gas_mean = rho_crit * Omega_b 

dens_max = 10 * rho_gas_mean



files_per_snapshot = 16

n_snap = 169

chem_type = 'HI'
chem_dens_name = f'{chem_type}_density'

X = 0.75984603480
if chem_type == 'HI': chem_fraction = X 
else: chem_fraction = 1 - X


n_file = 0
in_file_name = input_dir + f'{n_snap}.h5.{n_file}'
in_file = h5.File( in_file_name, 'r' )
current_z = in_file.attrs['Current_z'][0]

density = in_file['density'][...]
chem_density = in_file[chem_dens_name][...]
indices = density <= dens_max
# n_samples_local =  
density = density[indices]
chemn_density = chem_density[indices]
 




