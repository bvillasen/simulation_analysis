import sys, os, time
import numpy as np
import h5py as h5
from tools import *
from load_data import load_snapshot_data_distributed

use_mpi = False
if use_mpi:
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  n_procs = comm.Get_size()
else:
  rank = 0
  n_procs = 1

data_dir = '/gpfs/alpine/csc434/scratch/bvilasen/'
input_dir = data_dir + 'cosmo_sims/rescaled_P19/2048_50Mpc/snapshot_files/'
output_dir = data_dir + 'cosmo_sims/rescaled_P19/2048_50Mpc/reduced_snapshots/'
if rank == 0: create_directory( output_dir )

n_snap = 0


type = 'hydro'
# type = 'particles'

if rank == 0: print(f'Reducing Snapshot: {n_snap }' )

fields_hydro = [ 'density', 'temperature' ]
fields_particles = [ 'density' ]

if type == 'hydro': base_file_name = '.h5.'
if type == 'particles': base_file_name = '_particles.h5.'



files_snapshot = [f for f in listdir(input_dir) if f.find(f'{n_snap}{base_file_name}') == 0 ]
n_files = len( files_snapshot )

if rank == 0: print(f'N files per snapshot: {n_files}')



