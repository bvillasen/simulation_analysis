import os, sys
from os import listdir
from os.path import isfile, join, isdir
import numpy as np



def Link_Simulation_dirctories( src_params, dst_params ):
  param_names = dst_params[0]['parameters'].keys()
  n_params = len( param_names )
  n_src_sims = len( src_params )
  n_dst_sims = len( dst_params )
  dst_array = np.zeros([ n_dst_sims, n_params ])
  src_array = np.zeros([ n_src_sims, n_params ]) 

  for sim_id in range( n_src_sims ):
    sim_params = src_params[sim_id]['parameters']
    for param_id, param_name in enumerate(param_names):
      src_array[sim_id, param_id] = sim_params[param_name]
      
  for sim_id in range( n_dst_sims ):
    sim_params = dst_params[sim_id]['parameters']
    for param_id, param_name in enumerate(param_names):
      dst_array[sim_id, param_id] = sim_params[param_name]
      
  for sim_id in range( n_dst_sims ):
    dst_param_vals = dst_array[sim_id]
    diff = np.abs( src_array - dst_param_vals ).sum(axis=1)
    indx_src = np.where(diff == 0)[0]
    if len( indx_src ) == 0: src_sim_sir, src_id = None, None
    if len( indx_src ) == 1: src_sim_sir, src_id = src_params[indx_src[0]]['dir'], src_id
    dst_params[sim_id]['src_dir'] = src_sim_sir
    dst_params[sim_id]['src_id'] = src_sim_id
    



def Get_Grid_Parameter_Values( grid_dir ):
  sim_dirs = [f for f in listdir(grid_dir) if (isdir(join(grid_dir, f)) and (f[0] == 'S' ) ) ]
  sim_dirs.sort()
  grid_params = {}  
  for sim_id, sim_dir in enumerate(sim_dirs):
    grid_params[sim_id] = {}
    sim_dir = grid_dir + sim_dir
    grid_params[sim_id]['dir'] = sim_dir
    params_file = sim_dir + '/uvb_params.txt'
    file = open( params_file, 'r' )
    lines = file.readlines()
    params = {}
    for line in lines:
      param_name, param_val = line.split('=') 
      params[param_name] = float( param_val )
    grid_params[sim_id]['parameters'] = params
    file.close()
    
  return grid_params

