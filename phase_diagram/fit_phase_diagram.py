import os, sys
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from load_data import load_analysis_data
from phase_diagram_functions import fit_thermal_parameters_mcmc, get_density_tyemperature_values_to_fit

#Parse Command Parameters
args = sys.argv[1:]
n_args = len(args)

use_mpi = True
if use_mpi:
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  n_procs = comm.Get_size()
else:
  rank = 0
  n_procs = 1

data_dir = '/raid/bruno/data/'
input_dir = data_dir + 'cosmo_sims/sim_grid/512_50Mpc_pchw19/analysis_files/'
output_dir = data_dir + 'cosmo_sims/sim_grid/512_50Mpc_pchw19/figures/phase_diagram/'
fit_dir = input_dir + 'fit_mcmc/'
if rank == 0:
  create_directory( fit_dir )
  create_directory( output_dir )
if use_mpi: comm.Barrier()



indices = range(150)
indices_to_generate = split_indices( indices, rank,  n_procs )
if len(indices_to_generate) == 0: exit()
# print(f'Generating: {rank} {indices_to_generate}\n' ) 

for n_file in indices_to_generate:
  data = load_analysis_data( n_file, input_dir )
  values_to_fit = get_density_tyemperature_values_to_fit( data['phase_diagram'], delta_min=-1, delta_max=1, n_samples_line=50, fraction_enclosed=0.70 )
  fit_values = fit_thermal_parameters_mcmc( n_file, values_to_fit, fit_dir )


