import os, sys
from os import listdir
from os.path import isfile, join
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from load_data import load_snapshot_data_distributed
from tools import *

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'

# data_dir = '/gpfs/alpine/csc434/scratch/bvilasen/'
# data_dir = '/data/groups/comp-astro/bruno/'
data_dir = '/raid/bruno/data/'
# input_dir_cpu = data_dir + 'cosmo_sims/256_dm_50Mpc/output_files_cpu/'
# input_dir_cpu = data_dir + 'cosmo_sims/256_hydro_test/data_gpu_part_grav/'
# input_dir_gpu = data_dir + 'cosmo_sims/256_hydro_test/data_summit_gpu_mpi/'

data_name = 'cpu'
input_dir_0 = data_dir + f'cosmo_sims/256_hydro_50Mpc/data_summit_{data_name}/'
input_dir_1 = data_dir + f'cosmo_sims/256_hydro_50Mpc/data_caar_{data_name}/'
# input_dir_cpu = data_dir + 'sphere/data_summit_gpu/'
# input_dir_gpu = data_dir + 'sphere/output_files/'
# output_dir = data_dir + 'cosmo_sims/256_hydro_50Mpc/figures/'
# create_directory( output_dir ) 

precision = np.float64
Lbox = 50000.0    #kpc/h
n_cells = 256
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ n_cells, n_cells, n_cells ] #Size of the simulation grid


n_snaps = 60

diff_all_gas = []
z_all = []

# n_snapshot = 1
for n_snapshot in range(n_snaps):

  #Load DM data
  fields = [ 'density'  ]
  data_gas = load_snapshot_data_distributed( 'hydro', fields, n_snapshot, input_dir_0, box_size, grid_size,  precision, show_progess=True, print_fields=True )
  dens_gas_0 = data_gas['density']          
  z_0 = data_gas['Current_z']

  fields = [ 'density'  ]
  data_gas = load_snapshot_data_distributed( 'hydro', fields, n_snapshot, input_dir_1, box_size, grid_size,  precision, show_progess=True, print_fields=True )
  z_1 = data_gas['Current_z']
  dens_gas_1 = data_gas['density']    

  if np.abs( z_1 - z_0 ) > 0.001: 
    print ( 'ERROR: redshifts dont math' )
    exit(-1)
    
  z_all.append(z_0)
  diff_gas = np.abs( dens_gas_1 - dens_gas_0 ) / dens_gas_0
  diff_all_gas.append( diff_gas.max() )
  
  # diff_max_indx = np.where( diff == diff.max() )
  # print( '')
  # print( diff_max_indx )
  # print( '')
  


# 
# 
# nrows = 1
# ncols = 1
# fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8,6))
# 
# 
# ax.plot( z_all, diff_all )
# 
# ax.set_yscale('log')
# 
# ax.set_xlabel( r'$z$' )
# ax.set_ylabel( r'max $\Delta \rho / \rho$' )
# 
# 
# fileName = output_dir + 'delta_rho_max_gravity_gpu.png'
# fig.savefig( fileName, bbox_inches='tight', dpi=300)
# print('Saved Image: ', fileName)
# 
# 

