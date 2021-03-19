import os, sys
from os import listdir
from os.path import isfile, join
import h5py as h5
import numpy as np
from load_data import load_snapshot_data_distributed

#Example for Loading  Snapshot Data (Below)

data_dir = '/data/groups/comp-astro/bruno/'
input_dir = data_dir + 'cosmo_sims/halo_tests/256_hydro_50Mpc/output_files/'
# input_dir = data_dir + 'cosmo_sims/2048_hydro_50Mpc/output_files_pchw18/'


precision = np.float64
Lbox = 50000.0    #kpc/h
n_cells = 256
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ n_cells, n_cells, n_cells ] #Size of the simulation grid

n_snapshot = 169

#Load Gas data
fields = [ 'density' ]
data_gas = load_snapshot_data_distributed( 'hydro', fields, n_snapshot, input_dir, box_size, grid_size,  precision, show_progess=True )
current_z = data_gas['Current_z']  #redshift
density_gas = data_gas['density']  # h^2 Msun / kpc^3

Load DM data
fields = [ 'density', 'pos_x', 'pos_y', 'pos_z', 'particle_IDs' ]
fields = [ 'density',  ]
data_dm = load_snapshot_data_distributed( 'particles', fields, n_snapshot, input_dir, box_size, grid_size,  precision, show_progess=True, print_fields=True )
particle_mass = data_dm['particle_mass'] #h^-1 Msun 
density_dm = data_dm['density']          # h^2 Msun / kpc^3
pos_x = data_dm['pos_x']                 #h^-1 kpc
pos_y = data_dm['pos_y']                 #h^-1 kpc
pos_z = data_dm['pos_z']                 #h^-1 kpc
p_ids = data_dm['particle_IDs']      

