import os, sys
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from load_data import load_analysis_data
from spectra_functions import compute_optical_depth

#Parse Command Parameters
args = sys.argv[1:]
n_args = len(args)

use_mpi = False
if use_mpi:
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  n_procs = comm.Get_size()
else:
  rank = 0
  n_procs = 1

data_dir = '/raid/bruno/data/'
input_dir = data_dir + 'cosmo_sims/256_hydro_50Mpc/analysis_files/'
output_dir = data_dir + 'cosmo_sims/256_hydro_50Mpc/analysis_files/figures/'
if rank == 0: create_directory( output_dir )
if use_mpi: comm.Barrier()



n_snapshot = 0
data = load_analysis_data( n_snapshot, input_dir, phase_diagram=False, lya_statistics=True )
cosmology = data['cosmology']
box = data['box'] 
lya_statistics = data['lya_statistics']
n_skewers = lya_statistics['n_skewers']
Flux_mean = lya_statistics['Flux_mean']



