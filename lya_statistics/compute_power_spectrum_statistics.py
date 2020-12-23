import os, sys
import numpy as np
import h5py as h5
import pickle
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

parameters = sys.argv
if print_out: print( parameters )
for option in parameters:
  if option.find("uvb")    != -1: uvb = option[option.find('=')+1:]


# if print_out: print( f'Snapshot: {n_snap}' )


# uvb = 'pchw18'
# uvb = 'hm12'
# dataDir = '/home/bruno/Desktop/ssd_0/data/'
# dataDir = '/raid/bruno/data/'
dataDir = '/data/groups/comp-astro/bruno/'
simulation_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/'
input_dir = simulation_dir + 'transmited_flux_{0}_review/bootstraped_power_spectrum/'.format(uvb)
output_dir = simulation_dir + 'transmited_flux_{0}_review/bootstraped_power_spectrum/statistics/'.format(uvb)
if rank == 0: create_directory( output_dir )

snaps = [ 83, 90,  96, 102,  119, 124, 130, 136, 143, 151, 159, 169, ]
snaps_boss = [  96,  102, 106, 110,  114, 119, 124, 130, 136, 143, 151, 159 ]
snapshots = list( set( snaps_boss ).union(set(snaps)))
snapshots.sort()
print(snapshots)

# n_snap = 169
for n_snap in snapshots:

  print( f'Covariance Matrix Snapshot: {n_snap}' )

  file_name = input_dir + f'bootstrap_{n_snap}.h5'
  print( f'Writing to File: {file_name}' )
  file = h5.File( file_name, 'r')
  current_z = file.attrs['current_z']
  n_iterations = file.attrs['n_iterations']
  k_vals = file['k_vals'][...]
  ps_mean = file['mean'][...]


  ps_statistics_all = {}
  ps_statistics_all['current_z'] = current_z
  ps_statistics_all['n_iterations'] = n_iterations
  ps_statistics_all['mean'] = ps_mean
  ps_statistics_all['k_vals'] = k_vals

  ps_statistics_all['bootstrap'] = {}

  n_in_sample_list = [ 100, 500, 1000, 2500, 5000, 10000, 25000, 50000, 60000 ]
  for n_in_sample in n_in_sample_list:
    print(f'N in sample: {n_in_sample}')
    ps_statistics_all['bootstrap'][n_in_sample] = {}
    ps_mean_samples = file[str(n_in_sample)]['sample_mean'][...]
    n_k_samples = ps_mean_samples.shape[1]

    # Obtain distribution of P(k)
    n_bins_for_dist = 30
    ps_statistics = {}
    for i in range(n_k_samples):
      ps_statistics[i] = {}
      ps_statistics[i]['k_val'] = k_vals[i]
      p_vals = ps_mean_samples[:,i]
      delta_power = p_vals * k_vals[i]  / np.pi 
      ps_statistics[i]['mean']  = delta_power.mean() 
      ps_statistics[i]['sigma'] = delta_power.std() 
      bin_edges = np.linspace( delta_power.min()*0.99, delta_power.max()*1.01, n_bins_for_dist )
      power_hist, bin_edges = np.histogram( delta_power, bins=bin_edges )
      bin_centers = ( bin_edges[1:] + bin_edges[:-1] ) / 2
      power_hist = power_hist.astype(np.float) / power_hist.sum()
      ps_statistics[i]['distribution'] = {} 
      ps_statistics[i]['distribution']['bin_centers'] = bin_centers
      ps_statistics[i]['distribution']['histogram']   = power_hist
      

    matrix = np.zeros( [n_k_samples, n_k_samples] )
    for i in range(n_k_samples):
      for j in range(n_k_samples):
        vals_i = ps_mean_samples[:,i] * k_vals[i]  / np.pi
        vals_j = ps_mean_samples[:,j] * k_vals[j]  / np.pi
        mean_i = vals_i.mean()
        mean_j = vals_j.mean()
        n_samples_i = vals_i.shape[0]
        n_samples_j = vals_j.shape[0]
        if n_samples_i != n_samples_j: 
          print('ERROR: Number od samples mismatch')
          exit
        n_samples = n_samples_i
        matrix[i,j] = np.array([ (vals_i[k] - mean_i)*(vals_j[k] - mean_j) for k in range(n_samples)  ]).mean()
        
    ps_statistics_all['bootstrap'][n_in_sample]['statistics'] = ps_statistics
    ps_statistics_all['bootstrap'][n_in_sample]['cov_matrix'] = matrix


  file_name = output_dir + f'stats_{n_snap}.pkl'
  f = open( file_name, 'wb' )
  pickle.dump( ps_statistics_all, f)
  f.close()
  print ( f'Saved File: {file_name }' )


