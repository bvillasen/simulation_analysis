import os, sys
from os import listdir
from os.path import isfile, join, isdir




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

