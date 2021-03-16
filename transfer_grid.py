import os, sys
from os.path import isfile, join
import numpy
import h5py as h5


data_dir = '/data/groups/comp-astro/bruno/cosmo_sims/sim_grid/'
src_dir = data_dir + '1024_P19m_np4_nsim256/' 
dst_dir = data_dir + '1024_P19m_np4_nsim320/'



grid_dir = src_dir
sim_dirs = [f for f in listdir(grid_dir) if (isdir(join(grid_dir, f)) and (f[0] == 'S' ) ) ]