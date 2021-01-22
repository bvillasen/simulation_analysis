import os, sys
import numpy as np
import h5py as h5
import pickle
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from tools import *
from stats_functions import compute_distribution, get_highest_probability_interval

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
# 
# parameters = sys.argv
# if print_out: print( parameters )
# for option in parameters:
#   if option.find("uvb")    != -1: uvb = option[option.find('=')+1:]


# if print_out: print( f'Snapshot: {n_snap}' )

n_points = 2048

# uvb = 'pchw18'
uvb = 'hm12'

# dataDir = '/home/bruno/Desktop/ssd_0/data/'
dataDir = '/raid/bruno/data/'
# dataDir = '/data/groups/comp-astro/bruno/'
simulation_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/'
input_dir  = simulation_dir + 'transmited_flux_{0}_review/flux_power_spectrum_new/'.format(uvb)
output_dir = simulation_dir + 'transmited_flux_{0}_review/flux_power_spectrum_new/'.format(uvb)
if rank == 0: create_directory( output_dir )

snaps = [ 83, 90,  96, 102,  119, 124, 130, 136, 143, 151, 159, 169, ]
snaps_boss = [  96,  102, 106, 110,  114, 119, 124, 130, 136, 143, 151, 159 ]
snapshots = list( set( snaps_boss ).union(set(snaps)))
snapshots.sort()
print(snapshots)

# n_snap = 159
for n_snap in snapshots:



  file_name = input_dir + f'flux_ps_{n_snap}.h5'
  print( f'Loading File: {file_name}' )
  file = h5.File( file_name, 'r')
  current_z = file.attrs['current_z']
  vel_Hubble = file['vel_Hubble'][...]
  k_vals = file['k_vals'][...]
  ps_vals = file['flux_power_spectrum'][...]
  file.close()

  n_skewers, n_bins = ps_vals.shape
  vel_max = vel_Hubble.max()

  print(f'N Skewers: {n_skewers}    n_bins:{n_bins} ' )


  n_bins_for_distribution = 100
  fill_sum = 0.70

  ps_stats = {}

  # index = 6
  for index in range( 25 ):
    k_val = k_vals[index]
    vel = 2*np.pi / k_val
    stride = n_points * ( vel / vel_max ) 
    n_steps = int( 2048 / stride )
    stride = int( stride )
    ids_1d = ( np.arange( 0, n_steps, 1 ) * stride ).astype( np.int )
    n_1d = len( ids_1d )
    n_independent = n_1d**2
    print ( f' id: {index},  val: {k_val:.1e}    n_independent: {n_independent}'   )
    delta_vals = ps_vals[:, index] * k_val / np.pi 
    delta_mean = delta_vals.mean()
    delta_sigma = delta_vals.std()
    distribution, bin_centers = compute_distribution( delta_vals, n_bins_for_distribution, log=True )
    v_l, v_r, v_max, sum = get_highest_probability_interval( bin_centers, distribution, fill_sum, log=True, n_interpolate=1000)
    ps_stats[index] = {}
    ps_stats[index]['k_val'] = k_val
    ps_stats[index]['bin_centers'] = bin_centers
    ps_stats[index]['distribution'] = distribution
    ps_stats[index]['delta_mean'] = delta_mean
    ps_stats[index]['delta_sigma'] = delta_sigma
    ps_stats[index]['sigma_l'] = v_l
    ps_stats[index]['sigma_r'] = v_r
    ps_stats[index]['sigma_max'] = v_max
    ps_stats[index]['n_independent'] = n_independent
    
  
  n_indp_list = []
  k_list = []
  mean_list, sigma_list = [], []
  sigma_asim_l, sigma_asim_r = [], []
  for index in range( 25 ):
    n_indp_list.append( ps_stats[index]['n_independent'] )
    k_list.append( ps_stats[index]['k_val'] )
    mean_list.append( ps_stats[index]['delta_mean'] )
    sigma_list.append( ps_stats[index]['delta_sigma'] )
    sigma_asim_l.append( ps_stats[index]['sigma_l']  )
    sigma_asim_r.append( ps_stats[index]['sigma_r']  )

  n_independent = np.array( n_indp_list )
  k_array = np.array( k_list )
  mean_array = np.array( mean_list )
  sigma_array = np.array( sigma_list )
  sigma_l_array = np.array( sigma_asim_l )
  sigma_r_array = np.array( sigma_asim_r )

  ps_stats['current_z'] = current_z
  ps_stats['k_vals'] = k_array
  ps_stats['n_independent'] = n_independent
  ps_stats['delta_mean'] = mean_array
  ps_stats['delta_sigma'] = sigma_array
  ps_stats['delta_sigma_l'] = sigma_l_array
  ps_stats['delta_sigma_r'] = sigma_r_array

  file_name = output_dir + f'stats_{n_snap}.pkl'
  f = open( file_name, 'wb' )
  pickle.dump( ps_stats, f)
  f.close()
  print ( f'Saved File: {file_name }' )


