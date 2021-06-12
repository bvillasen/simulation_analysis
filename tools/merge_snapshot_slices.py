import os, sys, time
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


data_dir = '/raid/bruno/data/' 
# data_dir = '/data/groups/comp-astro/bruno/'
root_dir  = data_dir + f'cosmo_sims/sim_grid/1024_P19m_np4_nsim400/selected_snapshot_files_params_H/slices/'
output_dir = root_dir
create_directory( output_dir )

scale_He_vals = []
deltaZ_He_vals = []

n_snap = 2
slice_depth = 64
slice_start = 192

params_type = 'H' 
# params_type = 'He'

field = 'temperature'
field = 'HI_density'

extend = True

n_sims = 20

data_sims = {}
for sim_id in range( n_sims ):
  input_dir = root_dir + f'sim_{sim_id}/'
  file_name = input_dir  + 'data_sim.pkl'
  sim_data = Load_Pickle_Directory( file_name )
  data_sims[sim_id] = sim_data
  if params_type == 'He':
    scale_He  = sim_data['parameters']['scale_He']
    deltaZ_He = sim_data['parameters']['deltaZ_He']
    param_vals_x.append( scale_He )
    param_vals_y.append( deltaZ_He )
  in_file_name = input_dir + f'slice_{n_snap}_start{slice_start}_depth{slice_depth}.h5'
  in_file = h5.File( in_file_name,  'r' )
  current_z = in_file.attrs['current_z']
  slice = in_file[field][...]
  if extend:
    nz, ny, nx = slice.shape
    slice_extended = np.zeros( (nz, 2*ny, 2*nx) )
    for i in range( 2):
      for j in range(2):
        slice_extended[:, i*ny:(i+1)*ny,  j*nx:(j+1)*nx] = slice
    slice = slice_extended
  data_sims[sim_id]['slice'] = slice
  in_file.close()

param_vals_x  = np.array( list( set( param_vals_x ) ) )
param_vals_y = np.array( list( set( param_vals_y) ) )
param_vals_x.sort()
param_vals_y.sort()


# nrows, ncols = len(param_vals_y), len(param_vals_x)
# indices_2D = np.ones( (nrows, ncols), dtype=np.int ) * -1  
# params_x = np.zeros( (nrows, ncols) )
# params_y = np.zeros( (nrows, ncols) )
# 
# 
# for sim_id in data_sims:
#   sim_data = data_sims[sim_id]
#   scale_He  = sim_data['parameters']['scale_He']
#   deltaZ_He = sim_data['parameters']['deltaZ_He']
#   param_x = scale_He
#   param_y = deltaZ_He
#   indx_x = np.where( param_vals_x == param_x )[0]
#   indx_y = np.where( param_vals_y == param_y )[0]
#   if len(indx_x) != 1: print(f'ERROR: indx_x: {indx_x}' )
#   if len(indx_y) != 1: print(f'ERROR: indx_x: {indx_y}' ) 
#   indices_2D[indx_y,indx_x] = sim_id   
#   params_x[indx_y,indx_x] = param_x
#   params_y[indx_y,indx_x] = param_y
# 
# 
# slice = data_sims[0]['slice']
# nz, ny, nx = slice.shape
# n_vals_x = len(param_vals_x)
# n_vals_y = len(param_vals_y)
# 
# n_vals_x, n_vals_y = 2, 2
# 
# n_split_x, n_split_y = n_vals_x-1, n_vals_y-1 
# delta_x, delta_y = nx//n_split_x, ny//n_split_y
# print( f'Slice Shape: {slice.shape} ' )
# print( f'N split   x:{n_split_x}   y:{n_split_y}' )
# 
# slice_new = image = np.zeros_like( slice )
# 
# 
# time_start = time.time()
# for i in range( ny ):
#   for j in range( nx ):
#     id = i*nx + j
#     pos_x, pos_y = j, i
#     indx_x, indx_y = pos_x // delta_x, pos_y//delta_y
# 
#     if indx_x == n_vals_x - 1: indx_x -= 1
#     if indx_y == n_vals_y - 1: indx_y -= 1
# 
# 
#     id_ld = indices_2D[ indx_y,   indx_x   ]
#     id_lu = indices_2D[ indx_y+3, indx_x   ]
#     id_rd = indices_2D[ indx_y,   indx_x+4 ]
#     id_ru = indices_2D[ indx_y+3, indx_x+4 ] 
# 
#     val_ld = data_sims[id_ld]['slice'][:, i, j]
#     val_lu = data_sims[id_lu]['slice'][:, i, j]
#     val_rd = data_sims[id_rd]['slice'][:, i, j]
#     val_ru = data_sims[id_ru]['slice'][:, i, j]
# 
#     x, y = pos_x % delta_x, pos_y % delta_y
#     alpha_x, alpha_y = x / delta_x, y / delta_y
# 
# 
#     diff_x_d = val_rd - val_ld 
#     diff_x_u = val_ru - val_lu 
#     interp_x_d = val_ld + alpha_x*diff_x_d 
#     interp_x_u = val_lu + alpha_x*diff_x_u
#     interp = interp_x_d + alpha_y * ( interp_x_u - interp_x_d )
#     slice_new[:, i, nx-j-1] = interp
# 
#     if id%16 == 0: print_progress( id+1, nx*ny, time_start )
# 
# 
# out_file_name = output_dir + f'slice_interp_{n_snap}_start{slice_start}_depth{slice_depth}'
# 
# if extend: out_file_name += '_extended'
# out_file_name += '.h5'
# 
# outfile = h5.File( out_file_name, 'w' )
# outfile.attrs['current_z'] = current_z
# outfile.create_dataset( field, data=slice_new )
# outfile.close()
# print( f'Saved File: {out_file_name}' )
# 








