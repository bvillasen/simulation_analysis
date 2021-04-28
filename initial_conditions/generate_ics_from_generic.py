import os, sys
from os import listdir
from os.path import isfile, join
import h5py as h5
import numpy as np
import subprocess
import yt
#Extend path to inclide local modules
root_dir = os.path.dirname(os.getcwd())
sub_directories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(sub_directories)
from tools import *
from constants_cosmo import G_COSMO
from ics_particles import generate_ics_particles
from ics_grid import expand_data_grid_to_cholla

# Box Size
Lbox = 50000.0    #kpc
nPoints = 256
nBoxes  = 8

current_z = 100.
H0 = 67.66
Omega_b = 0.0497


h = H0 / 100
current_a = 1. /(current_z + 1)
a3 = current_a **3

rho_crit = 3 * (H0*1e-3)**2 / ( 8 * np.pi * G_COSMO ) / h**2
rho_gas_mean = rho_crit * Omega_b



data_dir = '/raid/bruno/data/'
# data_dir = '/data/groups/comp-astro/bruno/'
# data_dir = '/home/bruno/Desktop/ssd_0/data/'
input_dir = data_dir + f'cosmo_sims/ics/generic/{nPoints}_hydro_50Mpc/'
output_dir = data_dir + f'cosmo_sims/{nPoints}_hydro_50Mpc/ics_{nBoxes}_generic/'
print(f'Input Dir: {input_dir}' )
print(f'Output Dir: {output_dir}' )
create_directory( output_dir )

hydro = False
particles = True

in_file_name = input_dir + 'ics_generic.h5'
print( f'Loading File: {in_file_name }' )
in_file = h5.File( in_file_name, 'r' )

header = in_file['header']
level_max, level_min = header.attrs['levelmax'], header.attrs['levelmin']
if level_max != level_min:
  print( 'Non Uniform Grid not suportted')
  exit(-1)
print ( f' Uniform Grid Level: {level_max}    ')
grid_len_x = header['grid_len_x'][...][0]
grid_len_y = header['grid_len_y'][...][0]
grid_len_z = header['grid_len_z'][...][0]
grid_off_x = header['grid_off_x'][...][0]
grid_off_y = header['grid_off_y'][...][0]
grid_off_z = header['grid_off_z'][...][0]
grid_size = [ grid_len_x, grid_len_y, grid_len_z ]
grid_offset = [ grid_off_x, grid_off_y, grid_off_z ]
print( f' Grid Size: {grid_size } ' )
print( f' Grid Offset: {grid_offset } ' )

level = f'{level_max:03}'
print( f' Loading Level: {level} ' )

data_key = f'level_{level}_BA_rho'
print ( f'  Loading Field: {data_key}')
data = in_file[data_key][...] 
data_size = data.shape
n_ghost = [ int((data_size[i] - grid_size[i])/2) for i in range(3) ] 
data = data[ n_ghost[0]:-n_ghost[0], n_ghost[1]:-n_ghost[1], n_ghost[2]:-n_ghost[2]  ]
# data = data * 1e9
# 
# 
# 
# 
# 
# 
# 
# 
