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
out_file = h5.File( out_file_name, 'w' )


fields_to_load = [ 'density', 'HI_density', 'temperature']
for field in fields_to_load:
  #Load Gas data
  fields = [ field ]
  data_gas = load_snapshot_data_distributed( 'hydro', fields, n_snapshot, input_dir, box_size, grid_size,  precision, show_progess=True, proc_grid=proc_grid, print_fields=True )
  current_z = data_gas['Current_z']  #redshift
  field_gas = data_gas[field]  # h^2 Msun / kpc^3

  if field == 'density': field = 'gas_density'
  out_file.create_dataset( field, data=field_gas )
  print( f'Saved Field: {field}' )



fields_to_load = [ 'density']
for field in fields_to_load:
  #Load DM data
  fields = [ field ]
  data_dm = load_snapshot_data_distributed( 'particles', fields, n_snapshot, input_dir, box_size, grid_size,  precision, show_progess=True, proc_grid=proc_grid, print_fields=True )
  field_dm = data_dm[field]  # h^2 Msun / kpc^3

  if field == 'density': field = 'dm_density'
  out_file.create_dataset( field, data=field_dm )
  print( f'Saved Field: {field}' )


out_file.attrs['current_z'] = current_z
out_file.close()
print( f'Saved File: {out_file}')

#Load DM data
# fields = [ 'density', 'pos_x', 'pos_y', 'pos_z', 'particle_IDs' ]
# fields = [ 'density',  ]
# data_dm = load_snapshot_data_distributed( 'particles', fields, n_snapshot, input_dir, box_size, grid_size,  precision, show_progess=True, print_fields=True )
# particle_mass = data_dm['particle_mass'] #h^-1 Msun 
# density_dm = data_dm['density']          # h^2 Msun / kpc^3
# pos_x = data_dm['pos_x']                 #h^-1 kpc
# pos_y = data_dm['pos_y']                 #h^-1 kpc
# pos_z = data_dm['pos_z']                 #h^-1 kpc
# p_ids = data_dm['particle_IDs']      

