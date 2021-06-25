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
data_dir = '/data/groups/comp-astro/bruno/'
# data_dir = '/raid/bruno/data/'
input_dir_0 = data_dir + 'cosmo_sims/256_hydro_50Mpc/output_original/'
input_dir_1 = data_dir + 'cosmo_sims/256_hydro_50Mpc/output_new/'
output_dir = data_dir + 'cosmo_sims/256_hydro_50Mpc/figures/'
create_directory( output_dir ) 

precision = np.float64
Lbox = 50000.0    #kpc/h
n_cells = 256
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ n_cells, n_cells, n_cells ] #Size of the simulation grid


n_snaps = 60

fields = [ 'density'  ]

diff_all = []
z_all = []
for n_snapshot in range(n_snaps):

  #Load DM data
  data = load_snapshot_data_distributed( 'hydro', fields, n_snapshot, input_dir_0 + 'snapshot_files/', box_size, grid_size,  precision, show_progess=True, print_fields=True )
  dens_0 = data['density']          
  
  data = load_snapshot_data_distributed( 'hydro', fields, n_snapshot, input_dir_1 + 'snapshot_files', box_size, grid_size,  precision, show_progess=True, print_fields=True )
  dens_1 = data['density']    

  diff = np.abs( dens_0 - dens_1 ) / dens_0
  diff_all.append( diff.max() )
  
print( diff_all )