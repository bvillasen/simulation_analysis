import sys, os, time
import numpy as np
import h5py as h5
import pymc
import matplotlib.pyplot as plt
analysis_dir = os.path.dirname(os.getcwd()) + '/'
sys.path.append(analysis_dir + 'phase_diagram')
sys.path.append(analysis_dir + 'tools')
from tools import *
from data_thermal_history import data_thermal_history_Gaikwad_2020a, data_thermal_history_Gaikwad_2020b



def mcmc_model_1D( param_to_fit, comparable_data, comparable_grid, SG):
  parameters = SG.parameters
  param_name = parameters[param_to_fit]['name']
  param_vals = parameters[param_to_fit]['values']
  print(f' Fitting: {param_name}  {param_vals}')
  param_min = min(param_vals)
  param_max = max(param_vals)
  param_mid = ( param_max + param_min ) / 2.
  param_mcmc  = pymc.Uniform(param_name, param_min, param_max, value=param_mid )
  @pymc.deterministic( plot=False )
  def mcmc_model_1D( param_to_fit=param_to_fit, comparable_grid=comparable_grid, SG=SG, param_value=param_mcmc   ):
    mean_interp = Interpolate_Comparable_1D( param_to_fit, param_value,  comparable_grid, SG ) 
    return mean_interp
  densObsrv = pymc.Normal('T0', mu=mcmc_model_1D, tau=1./(comparable_data['sigma']**2), value=comparable_data['mean'], observed=True)
  return locals()
 
  

def Interpolate_Comparable_1D( param_id, param_value,  comparable_grid, SG ):
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