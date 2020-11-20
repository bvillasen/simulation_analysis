import sys, os, time
import numpy as np
import h5py as h5
import pymc
import matplotlib.pyplot as plt
analysis_dir = os.path.dirname(os.getcwd()) + '/'
sys.path.append(analysis_dir + 'phase_diagram')
sys.path.append(analysis_dir + 'lya_statistics')
sys.path.append(analysis_dir + 'tools')
from tools import *
from data_thermal_history import data_thermal_history_Gaikwad_2020a, data_thermal_history_Gaikwad_2020b
from data_optical_depth import *




def are_floats_equal( a, b, epsilon=1e-10 ):
  if np.abs( a - b ) < epsilon: return True
  else: return False
    
  
  
def Find_Parameter_Value_Near_IDs( param_id, param_value, parameters ):
  param_name = parameters[param_id]['name']    
  grid_param_values = np.array(parameters[param_id]['values'])
  param_min = grid_param_values.min() 
  param_max = grid_param_values.max()
  if param_value < param_min or param_value > param_max:
    print( f'ERROR: Paramneter Value outside {param_name} Range: [ {param_min} , {param_max} ] ')
    exit(-1)
  p_val_id_l, p_val_id_r = 0, 0
  diff_l, diff_r = -np.inf, np.inf
  for v_id, p_val in enumerate(grid_param_values):
    diff = p_val - param_value
    if diff > 0 and diff < diff_r: p_val_id_r, diff_r = v_id, diff
    if diff < 0 and diff > diff_l: p_val_id_l, diff_l = v_id, diff  
    if diff == 0: 
      if i < len(grid_param_values) -1:
        param_id_l = v_id
        param_id_r = v_id + 1
      else:
        param_id_l = v_id - 1
        param_id_r = v_id
      break
  return p_val_id_l, p_val_id_r
        


def Get_Parameter_Grid( param_values, parameters ):
  parameter_grid = {}
  for p_id, p_val in enumerate(param_values):
    parameter_grid[p_id] = {}
    v_id_l, v_id_r = Find_Parameter_Value_Near_IDs( p_id, p_val, parameters )
    parameter_grid[p_id]['v_id_l'] = v_id_l
    parameter_grid[p_id]['v_id_r'] = v_id_r
    parameter_grid[p_id]['v_l'] = parameters[p_id]['values'][v_id_l]
    parameter_grid[p_id]['v_r'] = parameters[p_id]['values'][v_id_r]
  return parameter_grid
  


def Get_Simulation_ID_From_Coordinates( sim_coords, SG ):
  grid = SG.Grid
  parameters = SG.parameters
  param_ids = parameters.keys()
  key = ''
  for param_id in param_ids:
    p_key = parameters[param_id]['key']
    key += f'_{p_key}{sim_coords[param_id]}'
  key = key[1:]
  sim_id = SG.coords[key]
  # print(sim_coords, key, sim_id)
  return sim_id

# linear_m = np.array([ -1.5, 2.5, -1.1, 3.0])
# linear_b = np.array([ 2.0, -0.5, 5.1, -4.0])
def Get_Value_From_Simulation( sim_coords, data_to_interpolate, field, sub_field, SG ):
  sim_id = Get_Simulation_ID_From_Coordinates( sim_coords, SG )
  sim = SG.Grid[sim_id]
  param_values = sim['parameter_values']
  # value = (linear_m * param_values + linear_b).sum()
  value = data_to_interpolate[sim_id][field][sub_field]
  return value
 


def Interpolate_MultiDim( p0, p1, p2, p3, data_to_interpolate, field, sub_field, SG, parameter_grid=None, param_id=None, sim_coords_before=None ):
  param_values = np.array([ p0, p1, p2, p3 ])
  n_param = len(param_values)
  if param_id == None: param_id = n_param - 1
  if sim_coords_before == None:  sim_coords_before = [ -1 for param_id in range(n_param)] 
  if parameter_grid == None: parameter_grid = Get_Parameter_Grid( param_values, SG.parameters )
  
  sim_coords_l = sim_coords_before.copy()
  sim_coords_r = sim_coords_before.copy()
  
  v_id_l = parameter_grid[param_id]['v_id_l']
  v_id_r = parameter_grid[param_id]['v_id_r']
  p_val_l = parameter_grid[param_id]['v_l']
  p_val_r = parameter_grid[param_id]['v_r']
  sim_coords_l[param_id] = v_id_l
  sim_coords_r[param_id] = v_id_r
  p_val = param_values[param_id]
  if p_val < p_val_l or p_val > p_val_r:
    print( ' ERROR: Parameter outside left and right values')
    exit()
  delta = ( p_val - p_val_l ) / ( p_val_r - p_val_l )  
  if param_id == 0:
    value_l = Get_Value_From_Simulation( sim_coords_l, data_to_interpolate, field, sub_field, SG )
    value_r = Get_Value_From_Simulation( sim_coords_r, data_to_interpolate, field, sub_field, SG )
    value_interp = delta*value_r + (1-delta)*value_l 
    return value_interp
  
  value_l = Interpolate_MultiDim( p0, p1, p2, p3, data_to_interpolate, field, sub_field, SG, parameter_grid=parameter_grid, param_id=param_id-1, sim_coords_before=sim_coords_l )
  value_r = Interpolate_MultiDim( p0, p1, p2, p3, data_to_interpolate, field, sub_field, SG, parameter_grid=parameter_grid, param_id=param_id-1, sim_coords_before=sim_coords_r )
  value_interp = delta*value_r + (1-delta)*value_l
  return value_interp



def Interpolate_Comparable_1D( param_id, param_value,  comparable_grid, SG ):
  parameters = SG.parameters
  param_name = parameters[param_id]['name']
  param_vals = parameters[param_id]['values']
  param_max = max(param_vals)
  param_min = min(param_vals) 
  if param_value < param_min or param_value > param_max:
    print( f'ERROR: Parameter Value outside {param_name} Range: [ {param_min} , {param_max} ] ')
    exit(-1)
  #find closest simulations 
  sim_ids = SG.sim_ids
  sim_id_l, sim_id_r = 0, 0
  diff_l, diff_r = -np.inf, np.inf
  for sim_id in sim_ids:
    sim_params = SG.Grid[sim_id]
    sim_param_val = sim_params['parameters'][param_name]
    diff = sim_param_val - param_value
    if diff > 0 and diff < diff_r: sim_id_r, diff_r = sim_id, diff
    if diff < 0 and diff > diff_l: sim_id_l, diff_l = sim_id, diff  
  param_l = SG.Grid[sim_id_l]['parameters'][param_name]
  param_r = SG.Grid[sim_id_r]['parameters'][param_name]    
  delta = param_value - param_l 
  mean_l = comparable_grid[sim_id_l]['mean']
  mean_r = comparable_grid[sim_id_r]['mean'] 
  mean = (mean_r - mean_l ) / ( param_r - param_l ) * delta  + mean_l 
  return mean


  
def Interpolate_Observable_1D( param_id, observable, param_value,  SG ):
  parameters = SG.parameters
  param_name = parameters[param_id]['name']
  param_vals = parameters[param_id]['values']
  param_max = max(param_vals)
  param_min = min(param_vals) 
  if param_value < param_min or param_value > param_max:
    print( f'ERROR: Paramneter Value outside {param_name} Range: [ {param_min} , {param_max} ] ')
    exit(-1)
  #find closest simulations 
  sim_ids = SG.sim_ids
  sim_id_l, sim_id_r = 0, 0
  diff_l, diff_r = -np.inf, np.inf
  for sim_id in sim_ids:
    sim_params = SG.Grid[sim_id]
    sim_param_val = sim_params['parameters'][param_name]
    diff = sim_param_val - param_value
    if diff > 0 and diff < diff_r: sim_id_r, diff_r = sim_id, diff
    if diff < 0 and diff > diff_l: sim_id_l, diff_l = sim_id, diff  
  param_l = SG.Grid[sim_id_l]['parameters'][param_name]
  param_r = SG.Grid[sim_id_r]['parameters'][param_name]    
  delta = param_value - param_l 
  mean_l = SG.Grid[sim_id_l]['analysis'][observable]
  mean_r = SG.Grid[sim_id_r]['analysis'][observable] 
  mean = (mean_r - mean_l ) / ( param_r - param_l ) * delta  + mean_l 
  return mean,  SG.Grid[sim_id_r]['analysis']['z'] 

def Get_Comparable_Tau():
  comparable_z, comparable_tau, comparable_sigma = [], [], []
  # Add data Becker 2013
  data_set = data_optical_depth_Becker_2013
  z   = data_set['z']
  tau = data_set['tau']
  sigma = data_set['tau_sigma']
  indices = z < 4.3
  comparable_z.append(z[indices])
  comparable_tau.append(tau[indices])
  comparable_sigma.append(sigma[indices])
  # Add data Keating 2020
  data_set = data_optical_depth_Keating_2020
  z   = data_set['z']
  tau = data_set['tau']
  sigma = data_set['tau_sigma']
  indices = z > 4.3
  comparable_z.append(z[indices])
  comparable_tau.append(tau[indices])
  comparable_sigma.append(sigma[indices])
  comparable = {}
  comparable['z'] = np.concatenate( comparable_z )
  comparable['mean'] = np.concatenate( comparable_tau )
  comparable['sigma'] = np.concatenate( comparable_sigma )
  return comparable

def Get_Comparable_T0_Gaikwad():
  data_sets = [ data_thermal_history_Gaikwad_2020b, data_thermal_history_Gaikwad_2020a ]
  z = np.concatenate( [ds['z'] for ds in data_sets ])
  data_mean  = np.concatenate( [ds['T0'] for ds in data_sets ])
  data_sigma = np.concatenate( [( ds['T0_sigma_plus'] + ds['T0_sigma_minus'] )*0.5  for ds in data_sets ] )
  comparable = {}
  comparable['z'] = z
  comparable['mean'] = data_mean
  comparable['sigma'] = data_sigma
  return comparable

def Get_Comparable_Composite_T0_tau():
  comparable_T0 = Get_Comparable_T0_Gaikwad()
  comparable_tau = Get_Comparable_Tau()
  comparable = {}
  comparable['T0'] = comparable_T0
  comparable['tau'] = comparable_tau
  comparable['T0+tau'] = {}
  for key in ['z', 'mean', 'sigma']:
    comparable['T0+tau'][key] = np.concatenate( [ comparable_T0[key], comparable_tau[key] ])
  return comparable


def Get_Comparable_Composite_T0_tau_from_Grid( comparable_data, SG ):
  comparable_T0_data = comparable_data['T0']
  comparable_tau_data = comparable_data['tau']
  comparable_T0_grid = Get_Comparable_T0_from_Grid( comparable_T0_data, SG )
  comparable_tau_grid = Get_Comparable_Tau_from_Grid( comparable_tau_data, SG )
  comparable_grid = {}
  sim_ids = SG.sim_ids
  for sim_id in sim_ids:
    comparable_grid[sim_id] = {}
    comparable_grid[sim_id]['T0'] = comparable_T0_grid[sim_id]
    comparable_grid[sim_id]['tau'] = comparable_tau_grid[sim_id]
    comparable_grid[sim_id]['T0+tau'] = {}
    for key in ['z', 'mean']:
      comparable_grid[sim_id]['T0+tau'][key] = np.concatenate( [ comparable_T0_grid[sim_id][key], comparable_tau_grid[sim_id][key] ])
  return comparable_grid

def Get_Comparable_Tau_from_Grid( comparable_data, SG ):
  z_data = comparable_data['z']
  sim_ids = SG.sim_ids
  comparable_grid = {}
  for sim_id in sim_ids:
    comparable_grid[sim_id] = {}
    sim_analysis = SG.Grid[sim_id]['analysis']
    z_sim = sim_analysis['z']
    indices = []
    for z in z_data:
      diff = np.abs( z_sim - z )
      indx = np.where(diff == diff.min())[0]
      indices.append( indx )
    indices = np.array( indices )
    comparable_grid[sim_id]['z'] = sim_analysis['z'][indices].flatten()
    comparable_grid[sim_id]['mean'] = sim_analysis['tau'][indices].flatten()    
  return comparable_grid

def Get_Comparable_T0_from_Grid( comparable_data, SG ):
  z_data = comparable_data['z']
  sim_ids = SG.sim_ids
  comparable_grid = {}
  for sim_id in sim_ids:
    comparable_grid[sim_id] = {}
    sim_analysis = SG.Grid[sim_id]['analysis']
    z_sim = sim_analysis['z']
    indices = []
    for z in z_data:
      diff = np.abs( z_sim - z )
      indx = np.where(diff == diff.min())[0]
      indices.append( indx )
    indices = np.array( indices )
    comparable_grid[sim_id]['z'] = sim_analysis['z'][indices].flatten()
    comparable_grid[sim_id]['mean'] = sim_analysis['T0'][indices].flatten()    
  return comparable_grid