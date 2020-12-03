import os, sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from load_data import load_snapshot_data_distributed


data_name = 'SIMPLE_PPMP_eta035'

# data_dir = '/home/bruno/Desktop/data/'
# data_dir = '/raid/bruno/data/'
data_dir = '/home/bruno/Desktop/ssd_0/data/'
enzo_dir = data_dir + 'cosmo_sims/enzo/256_hydro_50Mpc/h5_files/'
cholla_dir = data_dir + f'cosmo_sims/256_hydro_50Mpc/{data_name}/'



n_snapshot = 0
precision = np.float64
data_type = 'hydro'
fields = [ 'density' ]

Lbox = 5000    #kpc/h
proc_grid = [ 2, 2, 2]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 256, 256, 256 ] #Size of the simulation grid
subgrid = [ [0, 256], [0, 256], [0, 256] ] #Size of the volume to load


data = load_snapshot_data_distributed( n_snapshot, cholla_dir, data_type, fields, subgrid,  precision, proc_grid,  box_size, grid_size, show_progess=True )

