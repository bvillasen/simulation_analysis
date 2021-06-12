import os, sys
from os import listdir
from os.path import isfile, join
from shutil import copyfile
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import ImageGrid
import pickle
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from load_data import load_snapshot_data_distributed
from tools import *

use_mpi = True
if use_mpi :
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  nprocs = comm.Get_size()
else:
  rank = 0
  nprocs = 1
  
show_progess = False
if rank == 0: show_progess = True

sim_id = rank

# data_dir = '/raid/bruno/data/'
data_dir = '/data/groups/comp-astro/bruno/'
input_dir  = data_dir + f'cosmo_sims/sim_grid/1024_P19m_np4_nsim400/selected_snapshot_files/sim_{sim_id}/'
output_dir = data_dir + f'cosmo_sims/sim_grid/1024_P19m_np4_nsim400/selected_snapshot_files/slices/sim_{sim_id}/'
create_directory( output_dir )

n_points = 1024
Lbox = 50000.0 #kpc/h
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ n_points, n_points, n_points ]
precision = np.float32

n_snap = 6

field =  'temperature'  
data_gas = load_snapshot_data_distributed( 'hydro', [field], n_snap, input_dir, box_size, grid_size,  precision, show_progess=show_progess )
current_z = data_gas['Current_z']
data = data_gas['temperature']

slice_depth = 64
n_slices = n_points // slice_depth
print( f'N slices: {n_slices}' )


# for slice_id in range( n_slices ): 
slice_id = 3
slice_start = slice_id * slice_depth

start = max( 0, slice_start )
end   = min( n_points, slice_start+slice_depth )
print( f' Slice:  start:{start}   end:{end}' )
data_slice = data[slice_start:end, :, :] 

out_file_name = output_dir + f'slice_{n_snap}_start{slice_start}_depth{slice_depth}.h5'
outfile = h5.File( out_file_name, 'w' )
outfile.attrs['current_z'] = current_z
outfile.create_dataset( field, data=data_slice )
outfile.close()
print( f'Saved File: {out_file_name}' )


src = input_dir  + 'data_sim.pkl'
dst = output_dir + 'data_sim.pkl'
copyfile(src, dst)













