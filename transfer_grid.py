import os, sys
from os import listdir
from os.path import isfile, join, isdir
import numpy
import h5py as h5


data_dir = '/data/groups/comp-astro/bruno/cosmo_sims/sim_grid/'
src_dir = data_dir + '1024_P19m_np4_nsim256/' 
dst_dir = data_dir + '1024_P19m_np4_nsim320/'



grid_dir = src_dir
sim_dirs = [f for f in listdir(grid_dir) if (isdir(join(grid_dir, f)) and (f[0] == 'S' ) ) ]
sim_dirs.sort()

# for sim_dir in sim_dirs:
sim_dir = grid_dir + sim_dirs[0]
params_file = sim_dir + '/uvb_params.txt'
file = open( params_file, 'r' )
lines = file.read()

