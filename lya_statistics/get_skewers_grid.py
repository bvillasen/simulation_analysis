import sys, os, time
import numpy as np
import h5py as h5
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from tools import *
from load_data import load_snapshot_data_distributed
from spectra_functions import compute_optical_depth


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


inDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/output_files_{1}/'.format(n_points, uvb, cosmo_name )
output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/skewers_grid_{1}/'.format(n_points, uvb, cosmo_name )
if rank == 0: create_directory( output_dir )

snaps = [ 83, 90,  96, 102,  119, 124, 130, 136, 143, 151, 159, 169, ]
snaps_boss = [  96,  102, 106, 110,  114, 119, 124, 130, 136, 143, 151, 159 ]
snapshots = list( set( snaps_boss ).union(set(snaps)))
snapshots.sort()
# print(snapshots)

n_snapshot = snapshots[0]
snapshot_dir = output_dir + 'snapshot_{0:03}/'.format( n_snapshot ) 
if rank == 0: create_directory( snapshot_dir )

Lbox = 50000.0     #kpc/h
n_points = 2048
proc_grid = [ 8, 8, 8]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ n_points, n_points, n_points ]


indices = range( n_points )
split_indices = split_indices( indices, rank, nprocs, adjacent=True )
index_start, index_end = split_indices[0], split_indices[-1]+1
# print( f'Rank:{rank}   start:{index_start}   end:{index_end}' )

axis = 'x'

if axis == 'x':
  subgrid_x = [ 0, n_points ]
  subgrid_y = [ 0, n_points ]
  subgrid_z = [ index_start, index_end ]

if axis == 'y':
  subgrid_x = [ 0, n_points ]
  subgrid_y = [ 0, n_points ]
  subgrid_z = [ index_start, index_end ]

if axis == 'z':
  subgrid_x = [ 0, n_points ]
  subgrid_y = [ index_start, index_end ]
  subgrid_z = [ 0, n_points ]


subgrid = [ subgrid_x, subgrid_y, subgrid_z ]
if axis == 'x': vel_field = 'momentum_x'
if axis == 'y': vel_field = 'momentum_y'
if axis == 'z': vel_field = 'momentum_z'


data_type = 'hydro'
precision = np.float64
fields = [ 'density', 'HI_density', 'temperature', vel_field  ]
data = load_snapshot_data_distributed( n_snapshot, inDir, data_type, fields, subgrid,  precision, proc_grid,  box_size, grid_size, show_progess=show_progess )

H0        = data['H0'] * 1000 #km/s/Mpc
Omega_L   = data['Omega_L']
Omega_M   = data['Omega_M']
current_z = data['Current_z']
density     = data[data_type]['density']
HI_density  = data[data_type]['HI_density']
velocity    = data[data_type][vel_field] / density
temperature = data[data_type]['temperature']
subgrid_shape = density.shape
F_subgrid = np.ones( subgrid_shape ) * -1 


if rank == 0: 
  print( '\nStarting Calculation') 
  print( f' H0:        {H0}' )
  print( f' Omega_L:   {Omega_L}' )
  print( f' Omega_M:   {Omega_M}' )
  print( f' Current_z: {current_z}' )
  print( f' Subgrid shape: {subgrid_shape}' )
  
box = { 'Lbox':[ Lbox, Lbox, Lbox ] }
cosmology = { 'H0':H0, 'Omega_L':Omega_L, 'Omega_M':Omega_M, 'current_z':current_z }


if axis == 'x': n_los, n_i, n_j = subgrid_shape
if axis == 'y': n_i, n_los, n_j = subgrid_shape
if axis == 'z': n_i, n_j, n_los = subgrid_shape

n_skewers = n_i * n_j
n_processed = 0
text = f'N processed: {n_processed} / {n_skewers}    {n_processed/n_skewers*100}%'
if rank == 0: print_line_flush( text )
for i in range( n_i ):
  for j in range( n_j ):
    
    n_processed += 1
    
    if n_processed % ( n_skewers//128 ) == 0:
      text = f'N processed: {n_processed} / {n_skewers}    {n_processed/n_skewers*100}%       '
      if rank == 0: print_line_flush( text )
    
    
    if axis == 'x':
      los_HI_density = HI_density[:, i, j]
      los_velocity   = velocity[:, i, j]
      los_temperature = temperature[:, i, j]
    
    if axis == 'y':
      los_HI_density = HI_density[i, :, j]
      los_velocity   = velocity[i, :, j]
      los_temperature = temperature[i, :, j]
    
    if axis == 'z':
      los_HI_density = HI_density[i, j, :]
      los_velocity   = velocity[i, j, :]
      los_temperature = temperature[i, j, :]
    
    skewer_data = {}  
    skewer_data['HI_density']  = los_HI_density
    skewer_data['velocity']    = los_velocity
    skewer_data['temperature'] = los_temperature
    
    if j > 0: los_F = np.ones_like( los_HI_density )
    
    else:
      tau_los_data = compute_optical_depth( cosmology, box, skewer_data, space='redshift', method='error_function' )
      los_vel_hubble = tau_los_data['vel_Hubble']
      los_tau = tau_los_data['tau']
      los_F = np.exp( -los_tau )
    
    if len( los_F ) != n_los: print ('ERROR: Length of array does not match size of box')
    
    if axis == 'x': F_subgrid[:, i, j] = los_F
    if axis == 'y': F_subgrid[i, :, j] = los_F
    if axis == 'z': F_subgrid[i, j, :] = los_F
    
    
neg_indices = F_subgrid < 0
if neg_indices.sum() > 0: print ('ERROR: Negative Values in F_subgrid')
      
    

file_name = snapshot_dir + f'skewers_subgrid_{axis}_{rank}.h5'
file = h5.File( file_name, 'w' )
file.attrs['current_z'] = current_z
file.attrs['subgrid_shape'] = subgrid_shape


file.create_dataset( 'vel_hubble', data=los_vel_hubble  )
file.create_dataset( 'F_subgrid',  data=F_subgrid )
file.close()

print ( f'Saved File: {file_name}' )
    
    







comm.Barrier()


if rank == 0: print( 'Finished Succesfully' )


