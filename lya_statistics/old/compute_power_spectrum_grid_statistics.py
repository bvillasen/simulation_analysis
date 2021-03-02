import sys, os, time
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import pickle
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from flux_power_spectrum import get_skewer_flux_power_spectrum
from tools import *
from stats_functions import compute_distribution, get_highest_probability_interval


# dataDir = '/gpfs/alpine/proj-shared/ast149/'
# dataDir = '/data/groups/comp-astro/bruno/'
dataDir = '/raid/bruno/data/'

n_points = 2048
nx = n_points
ny = n_points
nz = n_points
ncells = nx * ny * nz

uvb = 'pchw18'
# uvb = 'hm12'

# cosmo_name = 'cosmo_3'
cosmo_name = ''


input_dir   = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/flux_power_spectrum_grid_{1}/'.format(n_points, uvb, cosmo_name )
figures_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/figures/'.format( n_points )
output_dir  = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/flux_power_spectrum_grid_{1}/'.format(n_points, uvb, cosmo_name )
create_directory( output_dir )




snaps = [ 83, 90,  96, 102,  119, 124, 130, 136, 143, 151, 159, 169, ]
snaps_boss = [  96,  102, 106, 110,  114, 119, 124, 130, 136, 143, 151, 159 ]
snapshots = list( set( snaps_boss ).union(set(snaps)))
snapshots.sort()
# print(snapshots)

n_snapshot = snapshots[-2]
snapshot_dir = output_dir + f'snap_{n_snapshot:03}/'
create_directory( snapshot_dir )

axis = 'x'
  
file_name = input_dir + f'ps_grid_{axis}_{n_snapshot}.h5'
print ( f'Loading File: {file_name}' )
file = h5.File( file_name, 'r' )
current_z = file.attrs['current_z']
F_mean_global = file.attrs['F_mean_global']
grid_size = file.attrs['grid_shape']
k_vals = file['k_vals'][...]
vel_Hubble = file['vel_Hubble'][...]
PS_grid = file['power_spectrum_grid'][...]
file.close()

vel_max = vel_Hubble[-1]
dv_Hubble = vel_Hubble[1] - vel_Hubble[0] 

ni, nj, n_bins = grid_size
n_skewers = ni * nj



delta_all = []
print( '\nRearranging Skewers')
for i in range( ni ):
  for j in range( nj ):
    delta_all.append( PS_grid[ i, j, :]* k_vals / np.pi )
delta_all = np.array( delta_all )




n_bins_for_distribution = 100
fill_sum = 0.70

ps_stats = {}

# index = 6
for index in range( 25 ):
  k_val = k_vals[index]
  delta_vals = delta_all[:, index] 
  delta_mean = delta_vals.mean()
  delta_sigma = delta_vals.std()
  distribution, bin_centers = compute_distribution( delta_vals, n_bins_for_distribution, log=True )
  v_l, v_r, v_max, sum = get_highest_probability_interval( bin_centers, distribution, fill_sum, log=True, n_interpolate=1000)

  vel = 2*np.pi / k_val
  stride = n_points * ( vel / vel_max ) 
  n_steps = int( 2048 / stride )
  stride = int( stride )
  ids_1d = ( np.arange( 0, n_steps, 1 ) * stride ).astype( np.int )
  n_1d = len( ids_1d )
  n_independent = n_1d**2

  n_groups = stride**2
  print ( f'\nk_index: {index},  stride: {stride}  N_independent: {n_independent},  N_groups: {n_groups}'  )
  n_processed = 0
  mean_samples = []
  start = time.time()
  for shift_x in range(stride):
    for shift_y in range(stride):

      n_processed += 1

      if n_processed % ( n_groups// min( 10, n_groups ) ) == 0:
        now = time.time()
        etr = ( n_groups - n_processed ) / n_processed * ( now - start ) 
        text = ' N processed: {0} / {1}    {2:.1f}%    ETR= {3:.2f} secs  '.format(n_processed, n_groups, n_processed/n_groups*100, etr)
        print_line_flush( text )

      ids_2d = [ ( id_y + shift_y, id_x + shift_x) for id_y in ids_1d for id_x in ids_1d ]
      ids_2d_linear = [ id_y*n_points + id_x for ( id_y, id_x ) in ids_2d ]

      group_vals = delta_vals[ids_2d_linear]
      group_mean = group_vals.mean()
      mean_samples.append( group_mean )
  mean_samples = np.array( mean_samples )


  ps_stats[index] = {}
  ps_stats[index]['k_val'] = k_val
  ps_stats[index]['delta_mean'] = delta_mean
  ps_stats[index]['delta_sigma'] = delta_sigma
  ps_stats[index]['bin_centers'] = bin_centers
  ps_stats[index]['distribution'] = distribution
  ps_stats[index]['sigma_left'] = v_l
  ps_stats[index]['sigma_max'] = v_max
  ps_stats[index]['sigma_right'] = v_r
  ps_stats[index]['mean_samples'] = mean_samples
  ps_stats[index]['n_independent'] = n_independent


print('\n')

file_name = output_dir + f'ps_grid_stats_{axis}_{n_snapshot}.pkl'
f = open( file_name, 'wb' )
pickle.dump( ps_stats, f)
f.close()
print ( f'Saved File: {file_name }' )





