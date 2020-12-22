import os, sys
import numpy as np
import h5py as h5
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from bootstrap_functions import bootstrap_sample_mean
from tools import *


use_mpi = False
if use_mpi :
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  nprocs = comm.Get_size()
else:
  rank = 0
  nprocs = 1
  
print_out = False
if rank == 0: print_out = True 

n_snap = 169

parameters = sys.argv
if print_out: print( parameters )
for option in parameters:
  if option.find("n_snap=") != -1: n_snap = int(option[option.find('=')+1:])


# if print_out: print( f'Snapshot: {n_snap}' )


# uvb = 'pchw18'
suvb = 'hm12'
# dataDir = '/home/bruno/Desktop/ssd_0/data/'
# dataDir = '/raid/bruno/data/'
dataDir = '/data/groups/comp-astro/bruno/'
simulation_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/'
input_dir = simulation_dir + 'transmited_flux_{0}_review/flux_power_spectrum/'.format(uvb)
output_dir = simulation_dir + 'transmited_flux_{0}_review/bootstraped_power_spectrum/'.format(uvb)
if rank == 0: create_directory( output_dir )

snaps = [ 83, 90,  96, 102,  119, 124, 130, 136, 143, 151, 159, 169, ]
snaps_boss = [  96,  102, 106, 110,  114, 119, 124, 130, 136, 143, 151, 159 ]
snapshots = list( set( snaps_boss ).union(set(snaps)))
snapshots.sort()
print(snapshots)


for n_snap in snapshots:

  print( f'Bootstraping Snapshot: {n_snap}' )

  # Load Power Spectrum Samples
  file_name = input_dir + f'flux_ps_{n_snap}.h5'
  if print_out: print( f'Loading  File: {file_name}' )
  file = h5.File( file_name, 'r')
  current_z = file.attrs['current_z'] 
  k_vals = file['k_vals'][...] 
  power_spectrum_all = file['flux_power_spectrum'][...]
  file.close()
  mean_all = power_spectrum_all.mean(axis=0)

  data = power_spectrum_all
  n_total = data.shape[0]
  print( f'N total in sampe: {n_total}' )

  n_iterations = 10000


  file_name = output_dir + f'bootstrap_{n_snap}.h5'
  print( f'Writing to File: {file_name}' )
  file = h5.File( file_name, 'w')
  file.attrs['current_z'] = current_z
  file.attrs['n_iterations'] = n_iterations
  file.create_dataset( 'k_vals', data=k_vals )
  file.create_dataset( 'mean',   data=mean_all )

  n_in_sample_list = [ 100, 500, 1000, 2500, 5000, 10000, 25000, 50000, 60000 ]
  for n_in_sample in n_in_sample_list:
    print(f'\nN in sample: {n_in_sample} ')
    sample_mean_all = bootstrap_sample_mean( n_iterations, n_in_sample, data, print_out )
    group = file.create_group( str(n_in_sample) )
    group.create_dataset( 'sample_mean', data=sample_mean_all )

  file.close()