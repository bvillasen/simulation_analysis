import os, sys, time
from pathlib import Path
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
cwd = os.getcwd()
cosmo_dir = cwd[: cwd.find('simulation_analysis')] + 'simulation_analysis/'
tools_dir = cosmo_dir + 'tools'
sys.path.append( tools_dir )
from tools import *
#Append analysis directories to path
extend_path()
from constants_cosmo import G_COSMO

use_mpi = True
if use_mpi:
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  n_procs = comm.Get_size()
else:
  rank = 0
  n_procs = 1

# data_dir = '/home/bruno/Desktop/data/'
# data_dir = '/home/bruno/Desktop/ssd_0/data/'
# data_dir = '/raid/bruno/data/'
# data_dir = '/data/groups/comp-astro/bruno/'
data_dir = '/gpfs/alpine/csc434/scratch/bvilasen/'
# input_dir = data_dir + 'cosmo_sims/2048_hydro_50Mpc/output_files_pchw18/'
# output_dir = data_dir + 'cosmo_sims/2048_hydro_50Mpc/neutral_fraction_pchw18/'
input_dir = data_dir + 'cosmo_sims/rescaled_P19/2048_50Mpc/snapshot_files/'
input_dir = data_dir + 'cosmo_sims/rescaled_P19/2048_50Mpc/neutral_fraction_P19m'
if rank == 0: create_directory( output_dir )



H0 = 67.66
Omega_b =  0.0497 
h = H0 / 100
X = 0.75984603480

rho_crit =  3*(H0*1e-3)**2/(8*np.pi*G_COSMO)/ h**2
rho_gas_mean = rho_crit * Omega_b 
dens_max = 3 * rho_gas_mean


files_per_snapshot = 512
indices_local = split_indices( range(files_per_snapshot), rank, n_procs )


snapshots = range( 0, 169 )


for n_snap in snapshots:

  if rank == 0: data_out = {}

  chem_types = [ 'HI', 'HeI', 'HeII'  ]
  for chem_indx, chem_type in enumerate(chem_types):
    chem_dens_name = f'{chem_type}_density'

    if chem_type == 'HI': chem_fraction = X 
    else: chem_fraction = 1 - X

    n_total_local = 0
    n_samples_local = 0
    fraction_sum = 0  

    for indx in indices_local:

      n_file = indx
      in_file_name = input_dir + f'{n_snap}.h5.{n_file}'
      in_file = h5.File( in_file_name, 'r' )
      # print( in_file.keys( ) )
      current_z = in_file.attrs['Current_z'][0]
      if rank == 0 and chem_indx == 0 and indx==0: print( f'\nn_sanp: {n_snap}    z:{current_z}' )

      density = in_file['density'][...]
      n_cells_local = np.prod( density.shape )
      chem_density = in_file[chem_dens_name][...]
      indices = density <= dens_max
      density = density[indices] * chem_fraction
      chem_density = chem_density[indices]
      dens_fraction = chem_density / density 
      fraction_sum += dens_fraction.sum()
      n_samples_local +=  indices.sum()
      n_total_local += n_cells_local



    #Send the phase diagram to root process
    fraction_all = comm.gather( fraction_sum, root=0 )
    n_local_all  = comm.gather( n_samples_local, root=0 )
    n_total_all  = comm.gather( n_total_local, root=0 )


    if rank == 0:
      fraction_all = np.array( fraction_all )
      n_local_all = np.array( n_local_all )
      n_total_all = np.array( n_total_all )
      fraction_sum_global = fraction_all.sum()
      n_local_global = n_local_all.sum()
      chem_fraction_global = fraction_sum_global / n_local_global
      n_total_global = n_total_all.sum()
      if chem_indx == 0: print( f'{chem_type} N_total:{n_total_global}    N_samples:{n_local_global}   Fraction:{n_local_global/n_total_global:.2f}  ' )  
      print( f'{chem_type} Fraction: {chem_fraction_global} ' )

      chem_data = { 'n_total':n_total_global, 'n_samples':n_local_global, 'dens_fraction':chem_fraction_global }
      data_out['current_z'] = current_z
      data_out[chem_type] = chem_data


  if rank == 0:
    output_file_name = output_dir + f'fraction_{n_snap}.pkl'
    # print( data_out )
    Write_Pickle_Directory( data_out, output_file_name )





