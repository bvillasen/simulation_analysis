import os, sys
from pathlib import Path
import numpy as np
import h5py as h5
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from flux_power_spectrum import get_skewer_flux_power_spectrum
from tools import *


input_dir = '/raid/bruno/data/cosmo_sims/256_hydro_50Mpc/analysis_files/'

use_mpi = False
if use_mpi:
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  n_procs = comm.Get_size()
else:
  rank = 0
  n_procs = 1

output_dir = input_dir + 'lya_statistics/'
if rank == 0: create_directory( output_dir )
if use_mpi: comm.Barrier()


files = [f for f in listdir(input_dir) if (isfile(join(input_dir, f)) and ( f.find('_analysis') > 0) ) ]
indices = [ '{0:03}'.format( int(file.split('_')[0]) ) for file in files ]

indices.sort()
n_files = len( files )
if rank == 0: print( f' N_Analysis_Files: {n_files}' )

indices_to_generate = split_indices( indices, rank,  n_procs )
if len(indices_to_generate) == 0: exit()

n_file = 50
# for n_file in indices_to_generate:
lya_stats_file = output_dir + f'lya_stats_{n_file}.pkl'
file_path = Path(lya_stats_file)
if file_path.is_file():
  print( f' Skiping File: {n_file} ') 
  # continue
  
infile_name = f'{n_file}_analysis.h5'
infile = h5.File( input_dir + infile_name, 'r' )

header = infile.attrs
H0 =  header['H0'][0]
current_z =  header['current_z'][0]

skewers_F_x = infile['lya_statistics']['skewers_x']['los_transmitted_flux_HeII'][...]
skewers_F_y = infile['lya_statistics']['skewers_y']['los_transmitted_flux_HeII'][...]
skewers_F_z = infile['lya_statistics']['skewers_z']['los_transmitted_flux_HeII'][...]
skewers_F_HeII = np.concatenate( [ skewers_F_x, skewers_F_y, skewers_F_z ], axis=0 )
F_mean_sim_HeII = infile['lya_statistics'].attrs['Flux_mean_HeII'][0]
F_mean_HeII = skewers_F_HeII.mean()
diff = np.abs( F_mean_HeII - F_mean_sim_HeII ) / F_mean_sim_HeII
if diff > 1e-6: 
  print( 'ERROR: Skewers F_mmean doesnt match simulation F_mean' )
  exit(-1)


skewers_F_x = infile['lya_statistics']['skewers_x']['los_transmitted_flux_HI'][...]
skewers_F_y = infile['lya_statistics']['skewers_y']['los_transmitted_flux_HI'][...]
skewers_F_z = infile['lya_statistics']['skewers_z']['los_transmitted_flux_HI'][...]
skewers_F = np.concatenate( [ skewers_F_x, skewers_F_y, skewers_F_z ], axis=0 )
F_mean_sim = infile['lya_statistics'].attrs['Flux_mean_HI'][0]
F_mean = skewers_F.mean()
diff = np.abs( F_mean - F_mean_sim ) / F_mean_sim
if diff > 1e-6: 
  print( 'ERROR: Skewers F_mmean doesnt match simulation F_mean' )
  exit(-1)

vel_Hubble_x = infile['lya_statistics']['skewers_x']['vel_Hubble'][...]
vel_Hubble_y = infile['lya_statistics']['skewers_y']['vel_Hubble'][...]
vel_Hubble_z = infile['lya_statistics']['skewers_z']['vel_Hubble'][...]
diff_y = ( np.abs( vel_Hubble_x - vel_Hubble_y ) / vel_Hubble_x ).sum()
diff_z = ( np.abs( vel_Hubble_x - vel_Hubble_z ) / vel_Hubble_x ).sum()
if diff_y > 1e-6 or diff_z > 1e-6 : 
  print( 'ERROR: Skewers vel_Hubble' )
  exit(-1)
vel_Hubble = vel_Hubble_x

n_skewers = skewers_F.shape[0]
if rank == 0: print( f' N Skewers: {n_skewers}' )


ps_sim     = infile['lya_statistics']['power_spectrum']['p(k)'][...]
k_vals_sim = infile['lya_statistics']['power_spectrum']['k_vals'][...]
indices = ps_sim > 1e-20
ps_sim = ps_sim[indices]
k_vals_sim = k_vals_sim[indices]

d_log_k = 0.1
ps_all = []
for skewer_id in range( n_skewers ):
  los_F = skewers_F[skewer_id]
  delta_F = los_F / F_mean
  bin_centers, skewer_power_spectrum = get_skewer_flux_power_spectrum(vel_Hubble, delta_F, d_log_k=d_log_k )
  ps_all.append( skewer_power_spectrum )
  
ps_all = np.array( ps_all )
ps_mean = ps_all.mean( axis= 0 )

