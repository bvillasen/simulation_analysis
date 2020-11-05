import os, sys
from pathlib import Path
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

input_dir = args[0]
fit_dir = input_dir + 'fit_mcmc/'

use_mpi = True
if use_mpi:
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  n_procs = comm.Get_size()
else:
  rank = 0
  n_procs = 1


if rank == 0: create_directory( fit_dir )
if use_mpi: comm.Barrier()


files = [f for f in listdir(input_dir) if (isfile(join(input_dir, f)) and ( f.find('_analysis') > 0) ) ]
indices = [ '{0:03}'.format( int(file.split('_')[0]) ) for file in files ]
indices.sort()
n_files = len( files )
if rank == 0: print( f' N_Analysis_Files: {n_files}' )


indices_to_generate = split_indices( indices, rank,  n_procs )
if len(indices_to_generate) == 0: exit()
# print(f'Generating: {rank} {indices_to_generate}\n' ) 

for n_file in indices_to_generate:
  fit_file = fit_dir + f'fit_{n_file}.pkl'
  file_path = Path(fit_file)
  if file_path.is_file():
    print( f' Skiping File: {n_file} ') 
    continue
  data = load_analysis_data( n_file, input_dir )
  values_to_fit = get_density_tyemperature_values_to_fit( data['phase_diagram'], delta_min=-1, delta_max=1, n_samples_line=50, fraction_enclosed=0.70 )
  fit_values = fit_thermal_parameters_mcmc( n_file, values_to_fit, fit_dir )


if use_mpi: comm.Barrier()
if rank == 0:
  files_fit = [f for f in listdir(fit_dir) if ( f.find('fit_') >= 0) ) ]
  n_files_fit = len( files_fit )
  if rank == 0: print( f' N_Fit_Files: {n_files_fit}' )
  if n_files_fit != n_files:
    print( f'ERROR: Fit files not match N_Files: {input_dir}' )
    exit(-1) 
  
  print(f'Successfully fit: {input_dir}')
