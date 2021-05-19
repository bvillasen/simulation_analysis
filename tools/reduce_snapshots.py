import sys, os, time
import numpy as np
import time
import h5py as h5
from tools import *

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

n_snapshots = 1

precision = np.float32


file_counter = 0
time_start = time.time()


n_snap = 0

type = 'hydro'
# type = 'particles'

if rank == 0: print(f'Reducing Snapshot: {n_snap }' )

fields_hydro = [ 'density', 'temperature' ]
fields_particles = [ 'density' ]

if type == 'hydro': 
  base_file_name = '.h5.'
  fields_list = fields_hydro
  
if type == 'particles': 
  base_file_name = '_particles.h5.'
fields_list = fields_particles




files_snapshot = [f for f in listdir(input_dir) if f.find(f'{n_snap}{base_file_name}') == 0 ]
n_files_per_snap = len( files_snapshot )

if rank == 0: print(f'N files per snapshot: {n_files_per_snap}')
if rank == 0: print( f'Splitting over {n_procs} processes ' )

indices_local = split_indices( range(n_files_per_snap), rank, n_procs )
n_files_local = len(indices_local)
n_total_local = n_snapshots * n_files_local


for file_id in indices_local:
  file_name = files_snapshot[file_id]
  in_file  = h5.File( input_dir + file_name, 'r' )
  out_file = h5.File( output_dir + file_name, 'w' )
    
  for key in in_file.attrs.keys():
    out_file.attrs[key] = in_file.attrs[key]
  
  for field in fields_list:
    data = in_file[field][...].astype( precision )
    out_file.create_dataset( field, data=data )
  
  
  
  
  
  in_file.close()
  out_file.close() 
  
  file_counter += 1
  if rank == 0: print_progress( file_counter, n_total_local, time_start )
  
  break
  







