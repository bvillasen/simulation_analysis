import sys, os, time
import numpy as np
import h5py as h5
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from flux_power_spectrum import get_skewer_flux_power_spectrum
from tools import *


dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/gpfs/alpine/proj-shared/ast149/'

n_points = 2048
nx = n_points
ny = n_points
nz = n_points
ncells = nx * ny * nz

uvb = 'pchw18'
# uvb = 'hm12'

# cosmo_name = 'cosmo_3'
cosmo_name = ''


input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/skewers_grid_{1}/'.format(n_points, uvb, cosmo_name )


snaps = [ 83, 90,  96, 102,  119, 124, 130, 136, 143, 151, 159, 169, ]
snaps_boss = [  96,  102, 106, 110,  114, 119, 124, 130, 136, 143, 151, 159 ]
snapshots = list( set( snaps_boss ).union(set(snaps)))
snapshots.sort()
# print(snapshots)

n_snapshot = snapshots[-2]
snapshot_dir = input_dir + 'snapshot_{0:03}/'.format( n_snapshot ) 
create_directory( snapshot_dir )

axis = 'x'

file_base_name = f'power_spectrum_subgrid_{axis}'
file_list, n_files = get_files_names( file_base_name, snapshot_dir )

allocated_memory = False

print('\nLoading Files...')
for index, file_name in enumerate( file_list ):
  
  text = f' Loading File: {index+1} / {n_files} '
  print_line_flush( text )
  
  file = h5.File( snapshot_dir + file_name, 'r' )
  current_z = file.attrs['current_z']
  F_mean_global = file.attrs['F_mean_global']
  subgrid_shape = file.attrs['subgrid_shape']
  k_vals = file['k_vals'][...]
  PS_subgrid = file['power_spectrum_subgrid']
  file.close()
  
  if not allocated_memory:
    n_i, n_j, n_bins = subgrid_shape
    full_size = ( n_i, n_j* n_files, n_bins )
    PS_grid  = np.ones( full_size ) * -1
    allocated_memory = True 

    


