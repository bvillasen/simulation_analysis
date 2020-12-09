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

# Box Size
Lbox = 50000.0    #kpc
nPoints = 256
nBoxes  = 8

# data_dir = '/raid/bruno/data/'
# data_dir = '/data/groups/comp-astro/bruno/'
data_dir = '/home/bruno/Desktop/ssd_0/data/'
input_dir = data_dir + f'cosmo_sims/nyx/{nPoints}_hydro_50Mpc/'
output_dir = data_dir + f'cosmo_sims/nyx/{nPoints}_hydro_50Mpc/h5_files/'
print(f'Input Dir: {input_dir}' )
print(f'Output Dir: {output_dir}' )
create_directory( output_dir )

file_names = [f for f in listdir(input_dir) if (f.find('plt') >= 0 )  ]
file_names.sort()

for nSnap,file_name in enumerate(file_names):

  ds = yt.load( input_dir + file_name )
  data = ds.all_data()
  h = ds.hubble_constant
  current_z = np.float(ds.current_redshift)
  current_a = 1./(current_z + 1)
  
  data_grid = ds.covering_grid( level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions )
  gas_dens = data_grid[('boxlib', 'density')] 
  gas_temp = data_grid[('boxlib', 'Temp')]    
  # gas_vel_x = data_grid[('boxlib', 'xmom')] / gas_dens 
  # gas_vel_y = data_grid[('boxlib', 'ymom')] / gas_dens
  # gas_vel_z = data_grid[('boxlib', 'zmom')] / gas_dens
  gas_dens = gas_dens / h / h / 1e9
  # gas_U = get_internal_energy( gas_temp ) * gas_dens 
  # gas_E = 0.5 * gas_dens * ( gas_vel_x**2 + gas_vel_y**2 +gas_vel_z**2 )  + gas_U
  # gas_mom_x = gas_vel_x * gas_dens
  # gas_mom_y = gas_vel_y * gas_dens
  # gas_mom_z = gas_vel_z * gas_dens
  
  
  
  # p_mass = data[('all', 'particle_mass')] * h
  # p_pos_x = data[('all', 'particle_position_x')].in_units('kpc')/current_a * h
  # p_pos_y = data[('all', 'particle_position_y')].in_units('kpc')/current_a * h
  # p_pos_z = data[('all', 'particle_position_z')].in_units('kpc')/current_a * h
  # p_vel_x = data[('all', 'particle_xvel')]
  # p_vel_y = data[('all', 'particle_yvel')]
  # p_vel_z = data[('all', 'particle_zvel')]
  
  
  out_file_name = output_dir + 'snapshot_{0:03}.h5'.format( nSnap )
  out_file = h5.File( out_file_name , 'w' )
  
  out_file.attrs['current_a'] = current_a
  out_file.attrs['current_z'] = current_z
  
  
  gas_group = out_file.create_group( 'gas' )
  gas_group.create_dataset( 'density', data=gas_dens )
  gas_group.create_dataset( 'temperature', data=gas_temp ) 
  
  
  
  
  print( f'Saved File: {out_file_name}')