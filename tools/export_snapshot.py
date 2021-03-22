import os, sys
from os import listdir
from os.path import isfile, join
import h5py as h5
import numpy as np
from load_data import load_snapshot_data_distributed
from tools import *
#Example for Loading  Snapshot Data (Below)

data_dir = '/data/groups/comp-astro/bruno/'
input_dir = data_dir + 'cosmo_sims/2048_hydro_50Mpc/output_files_pchw18/'
output_dir = data_dir + 'cosmo_sims/2048_hydro_50Mpc/snapshot_files/'
create_directory( output_dir )

precision = np.float32
Lbox = 50000.0    #kpc/h
n_cells = 2048
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ n_cells, n_cells, n_cells ] #Size of the simulation grid
proc_grid = [ 8, 8, 8 ]

n_snapshot = 169


out_file_name = output_dir + f'snapshot_{n_snapshot:03}.h5'
# out_file = h5.File( out_file_name, 'w' )
# 
# 
# fields_to_load = [ 'density', 'HI_density', 'temperature']
# for field in fields_to_load:
#   #Load Gas data
#   fields = [ field ]
#   data_gas = load_snapshot_data_distributed( 'hydro', fields, n_snapshot, input_dir, box_size, grid_size,  precision, show_progess=True, proc_grid=proc_grid, print_fields=True )
#   current_z = data_gas['Current_z']  #redshift
#   field_gas = data_gas[field]  # h^2 Msun / kpc^3
# 
#   if field == 'density': field = 'gas_density'
#   out_file.create_dataset( field, data=field_gas )
#   print( f'Saved Field: {field}' )
# 
# 
# 
# fields_to_load = [ 'density']
# for field in fields_to_load:
#   #Load DM data
#   fields = [ field ]
#   data_dm = load_snapshot_data_distributed( 'particles', fields, n_snapshot, input_dir, box_size, grid_size,  precision, show_progess=True, proc_grid=proc_grid, print_fields=True )
#   field_dm = data_dm[field]  # h^2 Msun / kpc^3
# 
#   if field == 'density': field = 'dm_density'
#   out_file.create_dataset( field, data=field_dm )
#   print( f'Saved Field: {field}' )
# 
# 
# out_file.attrs['current_z'] = current_z
# out_file.close()
# print( f'Saved File: {out_file}')



file = h5.File( out_file_name, 'r' )
print( f'Loading: {out_file_name} ' )

key = 'gas_density'
print( f'Loading key: {key}' )
# for key in file.keys():
field = file[key][...]
min  = field.min()
max  = field.max()
mean = field.mean()
print( f'{key}  min:{min}  max:{max}  mean:{mean}')


log_field = np.log10(field)
out_file_name = output_dir + f'{n_snapshot}.{key}_2048x2048x2048_float.raw'
log_field.astype(‘float32’).tofile(out_file_name)
print( f'Saved File: {out_file_name}')
