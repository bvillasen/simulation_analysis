import sys, os, time
import numpy as np
import h5py as h5
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from flux_power_spectrum import get_skewer_flux_power_spectrum
from tools import *


use_mpi = True

if use_mpi :
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  nprocs = comm.Get_size()
else:
  rank = 0
  nprocs = 1

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

show_progess = False
if rank == 0: show_progess = True


input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/skewers_grid_{1}/'.format(n_points, uvb, cosmo_name )


snaps = [ 83, 90,  96, 102,  119, 124, 130, 136, 143, 151, 159, 169, ]
snaps_boss = [  96,  102, 106, 110,  114, 119, 124, 130, 136, 143, 151, 159 ]
snapshots = list( set( snaps_boss ).union(set(snaps)))
snapshots.sort()
# print(snapshots)

n_snapshot = snapshots[-2]
snapshot_dir = input_dir + 'snapshot_{0:03}/'.format( n_snapshot ) 
if rank == 0: create_directory( snapshot_dir )

axis = 'x'
    
file_name = snapshot_dir + f'skewers_subgrid_{axis}_{rank}.h5'
file = h5.File( file_name, 'r' )
current_z     = file.attrs['current_z']
subgrid_shape = file.attrs['subgrid_shape']
vel_Hubble = file['vel_hubble'][...]
F_subgrid = file['F_subgrid'][...]
file.close()

F_mean_local = F_subgrid.mean()

if use_mpi:
  F_mean_global = comm.allreduce( F_mean_local, MPI.SUM ) / nprocs
else:
  F_mean_global = F_mean_local
  
tau = - np.log( F_mean_global  )
if rank == 0: 
  print( f'tau:    {tau}')
  print( f'F_mean: {F_mean_global}')
  print( f'current_z: {current_z}')
  print( f'subgrid_shape: {subgrid_shape}')
  
if axis == 'x': n_los, n_i, n_j = subgrid_shape
if axis == 'y': n_i, n_los, n_j = subgrid_shape
if axis == 'z': n_i, n_j, n_los = subgrid_shape


n_skewers = n_i * n_j
n_processed = 0
text = f'N processed: {n_processed} / {n_skewers}    {n_processed/n_skewers*100}%'
if rank == 0: print_line_flush( text )
start = time.time() 
allocated_memory = False
for i in range( n_i ):
  for j in range( n_j ):
    
    n_processed += 1
    
    if n_processed % ( n_skewers//128 ) == 0:
      now = time.time()
      etr = ( n_skewers - n_processed ) / n_processed * ( now - start ) / 60
      text = 'N processed: {0} / {1}    {2:.1f}%    ETR= {3:.2f} min  '.format(n_processed, n_skewers, n_processed/n_skewers*100, etr)
      if rank == 0: print_line_flush( text )
      
      if axis == 'x': los_F = F_subgrid[:, i, j]
      if axis == 'y': los_F = F_subgrid[i, :, j]
      if axis == 'z': los_F = F_subgrid[i, j, :]
      
      delta_F = ( los_F - F_mean_global ) / F_mean_global

      d_log_k = 0.1
      bin_centers, skewer_power_spectrum = get_skewer_flux_power_spectrum(vel_Hubble, delta_F, d_log_k=d_log_k )

      if not allocated_memory:
        n_bins = len(bin_centers)
        PS_subgrid = np.ones( [n_i, n_j, n_bins ] ) * -1
      
      PS_subgrid[i, j, :] = skewer_power_spectrum


comm.Barrier()

neg_indices = PS_subgrid < 0
if neg_indices.sum() > 0: print ('ERROR: Negative Values in F_subgrid')
      
    
file_name = snapshot_dir + f'power_spectrum_subgrid_{axis}_{rank:03}.h5'
file = h5.File( file_name, 'w' )
file.attrs['current_z'] = current_z
file.attrs['F_mean_global'] = F_mean_global
file.attrs['subgrid_shape'] = subgrid_shape


file.create_dataset( 'k_vals',  data=bin_centers )
file.create_dataset( 'power_spectrum_subgrid',  data=PS_subgrid )
file.close()

print ( f'Saved File: {file_name}' )



comm.Barrier()
if rank == 0: 
  print( f'N Bins: {n_bins}' )
  print( 'Finished Succesfully' )


