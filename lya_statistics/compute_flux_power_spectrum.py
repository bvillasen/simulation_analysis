import os, sys
import numpy as np
import h5py as h5
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from flux_power_spectrum import get_skewer_flux_power_spectrum
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
input_dir = simulation_dir + 'transmited_flux_{0}_review/los_F/'.format(uvb)
output_dir = simulation_dir + 'transmited_flux_{0}_review/flux_power_spectrum/'.format(uvb)
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

snaps = [ 83, 90,  96, 102,  119, 124, 130, 136, 143, 151, 159, 169, ]
snaps_boss = [  96,  102, 106, 110,  114, 119, 124, 130, 136, 143, 151, 159 ]
snapshots = list( set( snaps_boss ).union(set(snaps)))
snapshots.sort()
print(snapshots)


for n_snap in snapshots:

  file_name = input_dir + f'los_transmitted_flux_{n_snap}.h5'
  print( f'Loading: {file_name}' )
  file = h5.File( file_name, 'r' )
  current_z = file.attrs['current_z']
  n_skewers = file.attrs['n_skewers']  
  los_F = file['los_F'][...]
  vel_Hubble = file['vel_Hubble'][...]
  file.close()

  if los_F.shape[0] != n_skewers: 
    print( 'ERROR: Size of array does not match n_skewers')
    exit(-1)

  F_mean = los_F.mean()
  power_spectrum_all = []
  for i in range( n_skewers ):
    
    if i%(n_skewers//10)==0: 
      text = ' Skewer {0}/{1}    {2:.0f} %'.format(i, n_skewers,  float(i)/n_skewers*100)
      if rank == 0: print_line_flush( text )
      
    F = los_F[i]
    # Flux fluctuations
    delta_F = ( F - F_mean ) / F_mean 
            
    d_log_k = 0.1
    bin_centers, skewer_power_spectrum = get_skewer_flux_power_spectrum(vel_Hubble, delta_F, d_log_k=d_log_k )
    power_spectrum_all.append( skewer_power_spectrum )

  power_spectrum_all = np.array( power_spectrum_all )

  print( f'Writing Data shape:{power_spectrum_all.shape} ')


  file_name = output_dir + f'flux_ps_{n_snap}.h5'
  file = h5.File( file_name, 'w')
  file.attrs['current_z']  = current_z
  file.create_dataset('vel_Hubble', data=vel_Hubble )
  file.create_dataset('k_vals', data=bin_centers )
  file.create_dataset('flux_power_spectrum', data=power_spectrum_all )

  file.close()
  print( f'Saved File: {file_name}' )

