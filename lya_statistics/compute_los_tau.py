import os, sys
import numpy as np
import h5py as h5
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from flux_power_spectrum import get_skewer_flux_power_spectrum
from load_skewers import load_skewers_multiple_axis
from spectra_functions import compute_optical_depth
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
  
print_out = False
if rank == 0: print_out = True 

parameters = sys.argv
if print_out: print( parameters )
for option in parameters:
  if option.find("n_snap=") != -1: n_snap = int(option[option.find('=')+1:])

if print_out: print( f'Snapshot: {n_snap}' )

# 
uvb = 'pchw18'
# uvb = 'hm12'
# dataDir = '/home/bruno/Desktop/ssd_0/data/'
dataDir = '/raid/bruno/data/'
# dataDir = '/data/groups/comp-astro/bruno/'
simulation_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/'
input_dir = simulation_dir + 'skewers_{0}/'.format(uvb)
output_dir = simulation_dir + 'transmited_flux_{0}_review/los_F/'.format(uvb)
if rank == 0: create_directory( output_dir )


# Box parameters
Lbox = 50000.0 #kpc/h
nPoints = 2048
nx = nPoints
ny = nPoints
nz = nPoints
ncells = nx * ny * nz
box = {'Lbox':[ Lbox, Lbox, Lbox ] }

# Cosmology parameters
cosmology = {}
cosmology['H0'] = 67.66 
cosmology['Omega_M'] = 0.3111
cosmology['Omega_L'] = 0.6889

n_skewers_total = 60000
n_skewers_axis = n_skewers_total// 3 
n_skewers_list = [ n_skewers_axis, n_skewers_axis, n_skewers_axis ]
axis_list = [ 'x', 'y', 'z' ]
n_skewers_proc_list = [ ]


for i in range( len(axis_list) ):
  skewers_ids = range(n_skewers_list[i])
  skewers_ids_proc = split_indices( skewers_ids, rank, nprocs, adjacent=True )
  n_skewers_proc_list.append( len( skewers_ids_proc ))

if print_out: print(f"\nComputing LOS tau, s_nap:{n_snap}   n_skewers:{n_skewers_total}" )


# Load skewer data
skewer_dataset = load_skewers_multiple_axis( axis_list, n_skewers_proc_list, n_snap, input_dir, set_random_seed=False, print_out=print_out)
current_z = skewer_dataset['current_z']
cosmology['current_z'] = current_z
los_density = skewer_dataset['density']
los_HI_density = skewer_dataset['HI_density']
los_velocity = skewer_dataset['velocity']
los_temperature = skewer_dataset['temperature']


n_skewers = sum( n_skewers_proc_list )
skewers_ids = range(n_skewers)


processed_F = []
for i,skewer_id in enumerate(skewers_ids):

  if i%(n_skewers//10)==0: 
    text = ' Skewer {0}/{1}    {2:.0f} %'.format(i, n_skewers,  float(i)/n_skewers*100)
    # if rank == 0: print_line_flush( text )
    if rank == 0: print( text )

  skewer_data = {}  
  skewer_data['HI_density']  = los_HI_density[skewer_id]
  skewer_data['temperature'] = los_temperature[skewer_id]
  skewer_data['velocity']    = los_velocity[skewer_id]

  tau_los_data = compute_optical_depth( cosmology, box, skewer_data, space='redshift', method='error_function' )
  los_vel_hubble = tau_los_data['vel_Hubble']
  los_tau = tau_los_data['tau']
  los_F = np.exp( -los_tau )
  processed_F.append( los_F )


processed_F   = np.array( processed_F )

#Send the power spectrum to root process
if print_out: print( '\nGathering global data')
global_F   = comm.gather( processed_F, root=0 )



if rank == 0:

  global_F   = np.concatenate( global_F )
  n_processed = global_F.shape[0]
  print( f'n_processed: {n_processed},   ps_data shape: {global_F.shape}' )

  file_name = output_dir + f'los_transmitted_flux_{n_snap}.h5'
  file = h5.File( file_name, 'w')
  file.attrs['n_skewers'] = n_processed
  file.attrs['current_z'] = current_z
  file.create_dataset( 'los_F', data=global_F )
  file.create_dataset( 'vel_Hubble', data=los_vel_hubble )
  file.close()
  print( f'Saved File: {file_name}')



