import os, sys
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
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


uvb = 'pchw18'
# uvb = 'hm12'
# dataDir = '/home/bruno/Desktop/ssd_0/data/'
dataDir = '/raid/bruno/data/'
# dataDir = '/data/groups/comp-astro/bruno/'
simulation_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/'
input_dir = simulation_dir + 'transmited_flux_{0}_review/bootstraped_power_spectrum/'.format(uvb)
output_dir = simulation_dir + 'transmited_flux_{0}_review/bootstraped_power_spectrum/mean_ps_distribution/'.format(uvb)
if rank == 0: create_directory( output_dir )

snaps = [ 83, 90,  96, 102,  119, 124, 130, 136, 143, 151, 159, 169, ]
snaps_boss = [  96,  102, 106, 110,  114, 119, 124, 130, 136, 143, 151, 159 ]
snapshots = list( set( snaps_boss ).union(set(snaps)))
snapshots.sort()
print(snapshots)

n_snap = 169
# for n_snap in snapshots:

print( f'Covariance Matrix Snapshot: {n_snap}' )

file_name = input_dir + f'bootstrap_{n_snap}.h5'
print( f'Writing to File: {file_name}' )
file = h5.File( file_name, 'r')
current_z = file.attrs['current_z']
k_vals = file['k_vals'][...]
ps_mean = file['mean'][...]

n_in_sample = 60000
ps_mean_samples = file[str(n_in_sample)]['sample_mean'][...]
n_k_samples = ps_mean_samples.shape[1]


n_bins_for_dist = 30
ps_data = {}
for i in range(n_k_samples):
  ps_data[i] = {}
  ps_data[i]['k_val'] = k_vals[i]
  p_vals = ps_mean_samples[:,i]
  delta_power = p_vals * k_vals[i]  / np.pi 
  ps_data[i]['mean']  = delta_power.mean() 
  ps_data[i]['sigma'] = delta_power.std() 
  bin_edges = np.linspace( delta_power.min()*0.99, delta_power.max()*1.01, n_bins_for_dist )
  power_hist, bin_edges = np.histogram( delta_power, bins=bin_edges )
  bin_centers = ( bin_edges[1:] + bin_edges[:-1] ) / 2
  power_hist = power_hist.astype(np.float) / power_hist.sum()
  ps_data[i]['distribution'] = {} 
  ps_data[i]['distribution']['bin_centers'] = bin_centers
  ps_data[i]['distribution']['histogram']   = power_hist
  

for i in range(n_k_samples):
  
  fig_dir = output_dir + f'snap_{n_snap}/'
  create_directory(fig_dir)
  
  fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,5))
  
  bin_centers = ps_data[i]['distribution']['bin_centers']
  hist = ps_data[i]['distribution']['histogram']
  ax.plot( bin_centers, hist )
  
  file_name = fig_dir + f'distribution_k{i}.png'
  fig.savefig( file_name,  bbox_inches='tight', dpi=200 )
  print( f'Saved Figure: {file_name}' )
  
  

     





