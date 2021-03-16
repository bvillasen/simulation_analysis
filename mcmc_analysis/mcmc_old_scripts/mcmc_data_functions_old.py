import sys, os, time
import numpy as np
import h5py as h5
import pymc
import pickle
import matplotlib.pyplot as plt
analysis_dir = os.path.dirname(os.getcwd()) + '/'
sys.path.append(analysis_dir + 'phase_diagram')
sys.path.append(analysis_dir + 'lya_statistics')
sys.path.append(analysis_dir + 'tools')
from tools import *
from data_thermal_history import data_thermal_history_Gaikwad_2020a, data_thermal_history_Gaikwad_2020b
from data_optical_depth import *
from load_tabulated_data import load_power_spectrum_table, load_tabulated_data_boera, load_tabulated_data_viel, load_data_boss
from stats_functions import compute_distribution, get_highest_probability_interval

def Write_MCMC_Results( stats, MDL, params_mcmc,  stats_file, samples_file,  output_dir  ):
  cwd = os.getcwd()
  os.chdir( output_dir )

  f = open( stats_file, 'wb' )
  pickle.dump( stats, f)
  print ( f'Saved File: {stats_file}' )
  
  samples = {} 
  for p_id in params_mcmc.keys():
    param = params_mcmc[p_id]
    samples[p_id] = {}
    samples[p_id]['name'] = param['name']
    samples[p_id]['trace'] = param['sampler'].trace() 
  
  
  f = open( samples_file, 'wb' )
  pickle.dump( samples, f)
  print ( f'Saved File: {samples_file}' )
  os.chdir( cwd )
  return samples

def Get_Comparable_Composite_T0_tau(  factor_sigma_tau_becker=1, factor_sigma_tau_keating=1 ):
  comparable_T0 = Get_Comparable_T0_Gaikwad()
  comparable_tau = Get_Comparable_Tau(  factor_sigma_tau_becker=factor_sigma_tau_becker, factor_sigma_tau_keating=factor_sigma_tau_keating )
  comparable = {}
  comparable['T0'] = comparable_T0
  comparable['tau'] = comparable_tau
  comparable['T0+tau'] = {}
  for key in ['z', 'mean', 'sigma']:
    comparable['T0+tau'][key] = np.concatenate( [ comparable_T0[key], comparable_tau[key] ])
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
  

def Get_Comparable_Power_Spectrum_T0_from_Grid( comparable_data, SG ):
  comparable_ps = Get_Comparable_Power_Spectrum_from_Grid( comparable_data['P(k)']['separate'], SG )
  comparable_T0 = Get_Comparable_T0_from_Grid( comparable_data['T0'], SG )
  sim_ids = SG.sim_ids
  comparable = {}
  for sim_id in sim_ids:
    comparable[sim_id] = {}
    comparable[sim_id]['T0'] = comparable_T0[sim_id]
    comparable[sim_id]['P(k)'] = comparable_ps[sim_id]['P(k)']
    comparable[sim_id]['P(k)+T0'] = {}
    for key in ['mean']:
      comparable[sim_id]['P(k)+T0'][key] = np.concatenate([ comparable[sim_id]['P(k)'][key], comparable[sim_id]['T0'][key] ]) 
  return comparable

def Get_Comparable_Power_Spectrum_T0( ps_data_dir, z_min, z_max, data_sets, ps_range  ):
  comparable_ps = Get_Comparable_Power_Spectrum(  ps_data_dir, z_min, z_max, data_sets, ps_range )
  comparable_T0 = Get_Comparable_T0_Gaikwad()
  comparable = {}
  comparable['T0'] = comparable_T0
  comparable['P(k)'] = comparable_ps['P(k)']
  comparable['P(k)']['separate'] = comparable_ps['separate']
  comparable['P(k)+T0'] = {}
  for key in ['mean', 'sigma' ]:
      comparable['P(k)+T0'][key] = np.concatenate( ( comparable['P(k)'][key], comparable['T0'][key] ) )
  return comparable  
 
def Get_Comparable_Power_Spectrum( ps_data_dir, z_min, z_max, data_sets, ps_range ):
  print( f'Loading P(k) Data:' )
  dir_boss = ps_data_dir + 'data_power_spectrum_boss/'
  data_filename = dir_boss + 'data_table.py'
  data_boss = load_data_boss( data_filename )

  data_filename = ps_data_dir + 'data_power_spectrum_walther_2019/data_table.txt'
  data_walther = load_power_spectrum_table( data_filename )

  dir_data_boera = ps_data_dir + 'data_power_spectrum_boera_2019/'
  data_boera = load_tabulated_data_boera( dir_data_boera )

  data_dir_viel = ps_data_dir + 'data_power_spectrum_viel_2013/'
  data_viel = load_tabulated_data_viel( data_dir_viel)

  data_dir = { 'Boss':data_boss, 'Walther':data_walther, 'Boera':data_boera, 'Viel':data_viel }


  data_kvals, data_ps, data_ps_sigma, data_indices, data_z  = [], [], [], [], []

  sim_z, sim_kmin, sim_kmax = ps_range['z'], ps_range['k_min'], ps_range['k_max']

  ps_data = {}
  data_id = 0
  for data_index, data_name in enumerate(data_sets):
    print( f' Loading P(k) Data: {data_name}' )
    data_set = data_dir[data_name]
    keys = data_set.keys()
    n_indices = len(keys) - 1
    for index in range(n_indices):
      data = data_set[index]
      z = data['z']
      if z >= z_min and z <= z_max:
        diff = np.abs( sim_z - z )
        id_min = np.where( diff == diff.min() )[0][0]
        z_sim = sim_z[id_min]
        # print( f'data_z: {z:.1f} sim_z: {z_sim:.1f}')
        kmin = sim_kmin[id_min]
        kmax = sim_kmax[id_min]
        k_vals = data['k_vals']
        k_indices = np.where( (k_vals >= kmin) & (k_vals <= kmax) )
        k_vals = k_vals[k_indices]
        delta_ps = data['delta_power'][k_indices]
        delta_ps_sigma = data['delta_power_error'][k_indices]
        ps_data[data_id] = {'z':z, 'k_vals':k_vals, 'delta_ps':delta_ps, 'delta_ps_sigma':delta_ps_sigma }
        data_z.append( z )
        data_kvals.append( k_vals )
        data_ps.append( delta_ps )
        data_ps_sigma.append( delta_ps_sigma )
        data_id += 1
  k_vals_all         = np.concatenate( data_kvals )
  delta_ps_all       = np.concatenate( data_ps )
  delta_ps_sigma_all = np.concatenate( data_ps_sigma )
  ps_data_out = {'P(k)':{}, 'separate':ps_data }
  ps_data_out['P(k)']['k_vals'] = k_vals_all
  ps_data_out['P(k)']['mean']   = delta_ps_all
  ps_data_out['P(k)']['sigma']  = delta_ps_sigma_all

  n_data_points = len( k_vals_all )
  print( f' N data points: {n_data_points}' )
  return ps_data_out

def Get_Comparable_Power_Spectrum_from_Grid( comparable_data, SG ):
  
  print( 'Generating Simulation P(k) comparable data:')
  indices = comparable_data.keys()
  
  
  comparable_grid = {}
  sim_ids = SG.sim_ids
  for sim_id in sim_ids:
    comparable_grid[sim_id] = {}
    comparable_grid[sim_id]['P(k)'] = {}
    
    
    sim_data = SG.Grid[sim_id]['analysis']['power_spectrum']
    sim_z_all = sim_data['z']
    sim_k_vals_all = sim_data['k_vals']
    sim_ps_all = sim_data['ps_mean']
    
    sim_delta_all = []
    
    for index in indices:
      data = comparable_data[index]
      data_z = data['z']
      data_kvals = data['k_vals']
      data_delta_vals = data['delta_ps']
      diff = np.abs( sim_z_all - data_z )
      id_sim = np.where( diff == diff.min() )[0][0]
      sim_z = sim_z_all[id_sim]
      sim_kvals = sim_k_vals_all[id_sim]
      sim_ps = sim_ps_all[id_sim]
      sim_delta = sim_ps * sim_kvals / np.pi 
      sim_delta_interp = np.interp( data_kvals, sim_kvals, sim_delta )
      diff = ( sim_delta_interp - data_delta_vals ) / data_delta_vals
      sim_delta_all.append( sim_delta_interp )
      
    comparable_grid[sim_id]['P(k)']['mean'] = np.concatenate( sim_delta_all )
  
  n_points = len(  comparable_grid[0]['P(k)']['mean'] )
  print ( f' N sim points: {n_points}' )
  return comparable_grid  
    
    


def Get_Chi2( observables, params, comparable_grid, comparable_data, SG ):
  chi2_vals = {}
  for field in observables:
    obs_mean = Interpolate_3D(  params[0]['mean'], params[1]['mean'], params[2]['mean'], comparable_grid, field, 'mean', SG, clip_params=True ) 
    data_z     = comparable_data[field]['z']
    data_mean  = comparable_data[field]['mean']  
    data_sigma = comparable_data[field]['sigma']
    chi2 =  np.sum( ( ( obs_mean - data_mean ) / data_sigma )**2  )
    chi2_vals[field] = chi2
  return chi2_vals


def Sample_Power_Spectrum_from_Trace( param_samples, data_grid, SG, hpi_sum=0.7, n_samples=None ):

  print(f'\nSampling Power Spectrum')
  ps_samples = {}

  param_ids = param_samples.keys()
  n_param = len(param_ids )
  if not n_samples: n_samples = len( param_samples[0]['trace'] )
  print(f' N Samples: {n_samples}')

  n_z_ids = len( data_grid.keys() )

  for id_z in range( n_z_ids ):
    ps_data = data_grid[id_z]
    ps_samples[id_z] = {}
    ps_samples[id_z]['z'] = ps_data['z']

    samples = []
    for i in range( n_samples ):
      p_vals = []
      for p_id in param_ids:
        p_vals.append( param_samples[p_id]['trace'][i] )

      if n_param == 3: ps_interp = Interpolate_3D(  p_vals[0], p_vals[1], p_vals[2], ps_data, 'P(k)', 'mean', SG, clip_params=True ) 
      if n_param == 4: ps_interp = Interpolate_4D(  p_vals[0], p_vals[1], p_vals[2], p_vals[3], ps_data, 'P(k)', 'mean', SG, clip_params=True ) 
      samples.append( ps_interp )
    samples = np.array( samples ).T
    ps_mean = np.array([ ps_vals.mean() for ps_vals in samples ])
    ps_sigma = [ ]
    ps_lower, ps_higher = [], []
    for i in range( len( samples ) ):
      ps_sigma.append( np.sqrt(  ( (samples[i] - ps_mean[i])**2).mean()  ) )
      values = samples[i]
      n_bins = 100
      distribution, bin_centers = compute_distribution( values, n_bins, log=True )
      fill_sum = hpi_sum
      v_l, v_r, v_max, sum = get_highest_probability_interval( bin_centers, distribution, fill_sum, log=True, n_interpolate=1000)
      ps_lower.append( v_l )
      ps_higher.append( v_r )
    ps_sigma  = np.array( ps_sigma )
    ps_lower  = np.array( ps_lower )
    ps_higher = np.array( ps_higher )
    ps_samples[id_z]['mean'] = ps_mean
    ps_samples[id_z]['sigma'] = ps_sigma
    ps_samples[id_z]['k_vals'] = ps_data[0]['P(k)']['k_vals']
    ps_samples[id_z]['lower'] = ps_lower
    ps_samples[id_z]['higher'] = ps_higher
  return ps_samples


def Sample_T0_from_Trace( param_samples, data_grid, SG, hpi_sum=0.7, n_samples=None ):

  print(f'\nSampling Temperateure')
  temp_samples = {}

  param_ids = param_samples.keys()
  n_param = len(param_ids )
  if not n_samples: n_samples = len( param_samples[0]['trace'] )
  print(f' N Samples: {n_samples}')
  # 
  samples = []
  for i in range( n_samples ):
    p_vals = []
    for p_id in param_ids:
      p_vals.append( param_samples[p_id]['trace'][i] )

    if n_param == 3: temp_interp = Interpolate_3D(  p_vals[0], p_vals[1], p_vals[2], data_grid, 'T0', 'mean', SG, clip_params=True ) 
    if n_param == 4: temp_interp = Interpolate_4D(  p_vals[0], p_vals[1], p_vals[2], p_vals[3], data_grid, 'T0', 'mean', SG, clip_params=True ) 
    samples.append( temp_interp )
  samples = np.array( samples ).T
  temp_mean = np.array([ temp_vals.mean() for temp_vals in samples ])
  temp_sigma = [ ]
  temp_lower, temp_higher = [], []
  for i in range( len( samples ) ):
    temp_sigma.append( np.sqrt(  ( (samples[i] - temp_mean[i])**2).mean()  ) )
    values = samples[i]
    n_bins = 100
    distribution, bin_centers = compute_distribution( values, n_bins, log=True )
    fill_sum = hpi_sum
    v_l, v_r, v_max, sum = get_highest_probability_interval( bin_centers, distribution, fill_sum, log=True, n_interpolate=1000)
    temp_lower.append( v_l )
    temp_higher.append( v_r )
  temp_sigma  = np.array( temp_sigma )
  temp_lower  = np.array( temp_lower )
  temp_higher = np.array( temp_higher )
  temp_samples['mean']   = temp_mean
  temp_samples['sigma']  = temp_sigma
  temp_samples['z'] = data_grid[0]['z']
  temp_samples['lower']  = temp_lower
  temp_samples['higher'] = temp_higher
  return temp_samples






def Sample_Power_Spectrum_Multiple_Params( n_samples, params_all, data_grid, SG, sampling='gaussian', hpi_sum=0.7,  ):
  print(f'\nSampling Power Spectrum')
  ps_samples = {}
  
  n_param = len(params_all[0]['params'])
  
  n_z_ids = len( data_grid.keys() )

  for id_z in range( n_z_ids ):
    ps_data = data_grid[id_z]
    ps_samples[id_z] = {}
    ps_samples[id_z]['z'] = ps_data['z']

    samples = []
    for i in range( n_samples ):
      p_rand = []
      params_z = params_all[id_z]
      for p_id in params_z['params'].keys():
        if sampling == 'gaussian': p_rand.append( np.random.normal( params_z['params'][p_id]['mean'], params_z['params'][p_id]['sigma']  ) )
      if n_param == 3: ps_interp = Interpolate_3D(  p_rand[0], p_rand[1], p_rand[2], ps_data, 'P(k)', 'mean', SG, clip_params=True ) 
      if n_param == 4: ps_interp = Interpolate_3D(  p_rand[0], p_rand[1], p_rand[2], p_rand[3], ps_data, 'P(k)', 'mean', SG, clip_params=True ) 
      samples.append( ps_interp )
    samples = np.array( samples ).T
    ps_mean = np.array([ ps_vals.mean() for ps_vals in samples ])
    ps_sigma = [ ]
    ps_lower, ps_higher = [], []
    for i in range( len( samples ) ):
      ps_sigma.append( np.sqrt(  ( (samples[i] - ps_mean[i])**2).mean()  ) )
      values = samples[i]
      n_bins = 100
      distribution, bin_centers = compute_distribution( values, n_bins, log=True )
      fill_sum = hpi_sum
      v_l, v_r, v_max, sum = get_highest_probability_interval( bin_centers, distribution, fill_sum, log=True, n_interpolate=1000)
      ps_lower.append( v_l )
      ps_higher.append( v_r )
    ps_sigma  = np.array( ps_sigma )
    ps_lower  = np.array( ps_lower )
    ps_higher = np.array( ps_higher )
    ps_samples[id_z]['mean'] = ps_mean
    ps_samples[id_z]['sigma'] = ps_sigma
    ps_samples[id_z]['k_vals'] = ps_data[0]['P(k)']['k_vals']
    ps_samples[id_z]['lower'] = ps_lower
    ps_samples[id_z]['higher'] = ps_higher
  return ps_samples


def Sample_Power_Spectrum( n_samples, params, data_grid, SG, sampling='gaussian', hpi_sum=0.7 ):
  print(f'\nSampling Power Spectrum')
  ps_samples = {}
  n_param = len( params.keys() )

  n_z_ids = len( data_grid.keys() )

  for id_z in range( n_z_ids ):
    ps_data = data_grid[id_z]
    ps_samples[id_z] = {}
    ps_samples[id_z]['z'] = ps_data['z']

    samples = []
    for i in range( n_samples ):
      p_rand = []
      for p_id in params.keys():
        if sampling == 'gaussian': p_rand.append( np.random.normal( params[p_id]['mean'], params[p_id]['sigma']  ) )
        if sampling == 'uniform':
          p_min = params[p_id]['min']
          p_max = params[p_id]['max']
          p_rand.append( np.random.rand() * ( p_max - p_min) + p_min )
      if n_param == 3: ps_interp = Interpolate_3D(  p_rand[0], p_rand[1], p_rand[2], ps_data, 'P(k)', 'mean', SG, clip_params=True ) 
      if n_param == 4: ps_interp = Interpolate_4D(  p_rand[0], p_rand[1], p_rand[2], p_rand[3], ps_data, 'P(k)', 'mean', SG, clip_params=True ) 
      samples.append( ps_interp )
    samples = np.array( samples ).T
    ps_mean = np.array([ ps_vals.mean() for ps_vals in samples ])
    ps_sigma = [ ]
    ps_lower, ps_higher = [], []
    for i in range( len( samples ) ):
      ps_sigma.append( np.sqrt(  ( (samples[i] - ps_mean[i])**2).mean()  ) )
      values = samples[i]
      n_bins = 100
      distribution, bin_centers = compute_distribution( values, n_bins, log=True )
      fill_sum = hpi_sum
      v_l, v_r, v_max, sum = get_highest_probability_interval( bin_centers, distribution, fill_sum, log=True, n_interpolate=1000)
      ps_lower.append( v_l )
      ps_higher.append( v_r )
    ps_sigma  = np.array( ps_sigma )
    ps_lower  = np.array( ps_lower )
    ps_higher = np.array( ps_higher )
    ps_samples[id_z]['mean'] = ps_mean
    ps_samples[id_z]['sigma'] = ps_sigma
    ps_samples[id_z]['k_vals'] = ps_data[0]['P(k)']['k_vals']
    ps_samples[id_z]['lower'] = ps_lower
    ps_samples[id_z]['higher'] = ps_higher
  return ps_samples
  

def Sample_Observables( n_samples, observables, params, data_grid, SG  ):
  print(f'\nSampling Observables: {observables}')
  observables_samples = {}
  n_param = len( params.keys() )
  
  for observable in observables:
    observables_samples[observable] = {}
    observables_samples[observable]['samples'] = []

  for i in range(n_samples):
    p_rand = []
    for p_id in params.keys():
      p_rand.append( np.random.normal( params[p_id]['mean'], params[p_id]['sigma']  ) )
    for observable in observables:
      # obs_interp = Interpolate_MultiDim(  p_rand[0], p_rand[1], p_rand[2], p_rand[3], data_grid, observable, 'mean', SG, clip_params=True ) 
      if n_param == 3: obs_interp = Interpolate_3D(  p_rand[0], p_rand[1], p_rand[2], data_grid, observable, 'mean', SG, clip_params=True ) 
      if n_param == 4: obs_interp = Interpolate_3D(  p_rand[0], p_rand[1], p_rand[2], p_rand[3], data_grid, observable, 'mean', SG, clip_params=True ) 
      observables_samples[observable]['samples'].append( obs_interp )
  for observable in observables:
    obs_all = np.array(observables_samples[observable]['samples']).T
    obs_mean = [ obs_vals.mean() for obs_vals in obs_all ]
    obs_sigma = []
    for i in range( len(obs_all) ):
      obs_sigma.append( np.sqrt( (( obs_all[i] - obs_mean[i] )**2 ).mean() ) )
    observables_samples[observable]['mean'] = np.array( obs_mean )
    observables_samples[observable]['sigma'] = np.array( obs_sigma )   
    observables_samples[observable]['z'] = data_grid[0]['z']   
  
  return observables_samples

def Get_Data_Grid( fields, SG ):
  sim_ids = SG.sim_ids
  data_grid = {}
  for sim_id in sim_ids:
    data_grid[sim_id] = {}
    data_grid[sim_id]['z'] = SG.Grid[sim_id]['analysis']['z']
    for field in fields:
      data_grid[sim_id][field] = {}
      data_grid[sim_id][field]['mean'] = SG.Grid[sim_id]['analysis'][field]
  return data_grid

def Get_Data_Grid_Power_spectrum( z_vals, SG ):
  data_grid = {}
  for id_z, z_val in enumerate( z_vals ):
    data_grid[id_z] = {}
    data_grid[id_z]['z'] = z_val
    sim_ids = SG.sim_ids
    for sim_id in sim_ids:
      sim_ps_data = SG.Grid[sim_id]['analysis']['power_spectrum']
      sim_z_vals = sim_ps_data['z']
      diff_z = np.abs( sim_z_vals - z_val )
      diff_min = diff_z.min()
      if diff_min > 0.05: print( f'Warning: Large Z difference: {diff_min}')
      index = np.where( diff_z == diff_min )[0][0]
      k_vals = sim_ps_data['k_vals'][index]
      ps_vals = sim_ps_data['ps_mean'][index]
      delta_ps = ps_vals * k_vals / np.pi
      data_grid[id_z][sim_id] = { 'P(k)':{} } 
      data_grid[id_z][sim_id]['P(k)']['mean'] = delta_ps
      data_grid[id_z][sim_id]['P(k)']['k_vals'] = k_vals
  return data_grid

def are_floats_equal( a, b, epsilon=1e-10 ):
  if np.abs( a - b ) < epsilon: return True
  else: return False
    
  
  
def Find_Parameter_Value_Near_IDs( param_id, param_value, parameters, clip_params=False ):
  param_name = parameters[param_id]['name']    
  grid_param_values = np.array(parameters[param_id]['values'])
  param_min = grid_param_values.min() 
  param_max = grid_param_values.max()
  n_param_values = len( grid_param_values )
  # print( f' Param_id:{param_id}   value:{param_value}' )
  if n_param_values == 1:
    p_val_id_l,  p_val_id_r = 0, 0
    return p_val_id_l, p_val_id_r  
  if clip_params:
    if param_value < param_min: param_value = param_min
    if param_value > param_max: param_value = param_max
  else:  
    if param_value < param_min or param_value > param_max:
      print( f'ERROR: Paramneter Value outside {param_name} Range: [ {param_min} , {param_max} ] value:{param_value}')
      exit(-1)
  if are_floats_equal( param_value, param_min ):
    p_val_id_l,  p_val_id_r = 0, 1
    return p_val_id_l, p_val_id_r
  if are_floats_equal( param_value, param_max ):
    p_val_id_l,  p_val_id_r = n_param_values-2, n_param_values-1
    return p_val_id_l, p_val_id_r
  p_val_id_l, p_val_id_r = 0, 0
  diff_l, diff_r = -np.inf, np.inf
  for v_id, p_val in enumerate(grid_param_values):
    diff = p_val - param_value
    if diff > 0 and diff < diff_r: p_val_id_r, diff_r = v_id, diff
    if diff < 0 and diff > diff_l: p_val_id_l, diff_l = v_id, diff  
  if p_val_id_l == p_val_id_r: print('ERROR: Same values for left and right')
  return p_val_id_l, p_val_id_r
        


def Get_Parameter_Grid( param_values, parameters, clip_params=False ):
  parameter_grid = {}
  for p_id, p_val in enumerate(param_values):
    parameter_grid[p_id] = {}
    # print( f' Param_id:{p_id}   value:{p_val}' )
    v_id_l, v_id_r = Find_Parameter_Value_Near_IDs( p_id, p_val, parameters, clip_params=clip_params )
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
  return sim_id


def Get_Value_From_Simulation( sim_coords, data_to_interpolate, field, sub_field, SG ):
  sim_id = Get_Simulation_ID_From_Coordinates( sim_coords, SG )
  sim = SG.Grid[sim_id]
  param_values = sim['parameter_values']
  value = data_to_interpolate[sim_id][field][sub_field]
  return value
 


def Interpolate_MultiDim( param_values, data_to_interpolate, field, sub_field, SG, clip_params=False, parameter_grid=None, param_id=None, sim_coords_before=None ):
  n_param = len(param_values)
  if param_id == None: param_id = n_param - 1
  if sim_coords_before == None:  sim_coords_before = [ -1 for param_id in range(n_param)] 
  print( f' n_param: {n_param}  Param_id:{param_id}   value:{param_values}' )
  if parameter_grid == None: parameter_grid = Get_Parameter_Grid( param_values, SG.parameters, clip_params=clip_params )
  
  sim_coords_l = sim_coords_before.copy()
  sim_coords_r = sim_coords_before.copy()
  
  v_id_l = parameter_grid[param_id]['v_id_l']
  v_id_r = parameter_grid[param_id]['v_id_r']
  p_val_l = parameter_grid[param_id]['v_l']
  p_val_r = parameter_grid[param_id]['v_r']
  sim_coords_l[param_id] = v_id_l
  sim_coords_r[param_id] = v_id_r
  p_val = param_values[param_id]
  
  if clip_params:
    if p_val < p_val_l: p_val = p_val_l
    if p_val > p_val_r: p_val = p_val_r
  else:      
    if p_val < p_val_l or p_val > p_val_r:
      print( ' ERROR: Parameter outside left and right values')
      exit()
  if p_val_l == p_val_r: delta = 0.5
  else: delta = ( p_val - p_val_l ) / ( p_val_r - p_val_l )  
  if param_id == 0:
    value_l = Get_Value_From_Simulation( sim_coords_l, data_to_interpolate, field, sub_field, SG )
    value_r = Get_Value_From_Simulation( sim_coords_r, data_to_interpolate, field, sub_field, SG )
    value_interp = delta*value_r + (1-delta)*value_l 
    return value_interp
  
  value_l = Interpolate_MultiDim( param_values, data_to_interpolate, field, sub_field, SG, parameter_grid=parameter_grid, param_id=param_id-1, sim_coords_before=sim_coords_l, clip_params=clip_params )
  value_r = Interpolate_MultiDim( param_values, data_to_interpolate, field, sub_field, SG, parameter_grid=parameter_grid, param_id=param_id-1, sim_coords_before=sim_coords_r, clip_params=clip_params )
  value_interp = delta*value_r + (1-delta)*value_l
  return value_interp


def Interpolate_4D( p0, p1, p2, p3, data_to_interpolate, field, sub_field, SG, clip_params=False, parameter_grid=None, param_id=None, sim_coords_before=None ):
  param_values = np.array([ p0, p1, p2, p3 ])
  n_param = len(param_values)
  if param_id == None: param_id = n_param - 1
  if sim_coords_before == None:  sim_coords_before = [ -1 for param_id in range(n_param)] 
  if parameter_grid == None: parameter_grid = Get_Parameter_Grid( param_values, SG.parameters, clip_params=clip_params )
  
  sim_coords_l = sim_coords_before.copy()
  sim_coords_r = sim_coords_before.copy()
  
  v_id_l = parameter_grid[param_id]['v_id_l']
  v_id_r = parameter_grid[param_id]['v_id_r']
  p_val_l = parameter_grid[param_id]['v_l']
  p_val_r = parameter_grid[param_id]['v_r']
  sim_coords_l[param_id] = v_id_l
  sim_coords_r[param_id] = v_id_r
  p_val = param_values[param_id]
  
  if clip_params:
    if p_val < p_val_l: p_val = p_val_l
    if p_val > p_val_r: p_val = p_val_r
  else:      
    if p_val < p_val_l or p_val > p_val_r:
      print( ' ERROR: Parameter outside left and right values')
      exit()
  if p_val_l == p_val_r: delta = 0.5
  else: delta = ( p_val - p_val_l ) / ( p_val_r - p_val_l )  
  if param_id == 0:
    value_l = Get_Value_From_Simulation( sim_coords_l, data_to_interpolate, field, sub_field, SG )
    value_r = Get_Value_From_Simulation( sim_coords_r, data_to_interpolate, field, sub_field, SG )
    value_interp = delta*value_r + (1-delta)*value_l 
    return value_interp
  
  value_l = Interpolate_4D( p0, p1, p2, p3, data_to_interpolate, field, sub_field, SG, parameter_grid=parameter_grid, param_id=param_id-1, sim_coords_before=sim_coords_l, clip_params=clip_params )
  value_r = Interpolate_4D( p0, p1, p2, p3, data_to_interpolate, field, sub_field, SG, parameter_grid=parameter_grid, param_id=param_id-1, sim_coords_before=sim_coords_r, clip_params=clip_params )
  value_interp = delta*value_r + (1-delta)*value_l
  return value_interp


def Interpolate_3D( p0, p1, p2, data_to_interpolate, field, sub_field, SG, clip_params=False, parameter_grid=None, param_id=None, sim_coords_before=None ):
  param_values = np.array([ p0, p1, p2 ])
  n_param = len(param_values)
  if param_id == None: param_id = n_param - 1
  if sim_coords_before == None:  sim_coords_before = [ -1 for param_id in range(n_param)] 
  if parameter_grid == None: parameter_grid = Get_Parameter_Grid( param_values, SG.parameters, clip_params=clip_params )
  
  sim_coords_l = sim_coords_before.copy()
  sim_coords_r = sim_coords_before.copy()
  
  v_id_l = parameter_grid[param_id]['v_id_l']
  v_id_r = parameter_grid[param_id]['v_id_r']
  p_val_l = parameter_grid[param_id]['v_l']
  p_val_r = parameter_grid[param_id]['v_r']
  sim_coords_l[param_id] = v_id_l
  sim_coords_r[param_id] = v_id_r
  p_val = param_values[param_id]
  
  if clip_params:
    if p_val < p_val_l: p_val = p_val_l
    if p_val > p_val_r: p_val = p_val_r
  else:      
    if p_val < p_val_l or p_val > p_val_r:
      print( ' ERROR: Parameter outside left and right values')
      exit()
  if p_val_l == p_val_r: delta = 0.5
  else: delta = ( p_val - p_val_l ) / ( p_val_r - p_val_l )  
  if param_id == 0:
    value_l = Get_Value_From_Simulation( sim_coords_l, data_to_interpolate, field, sub_field, SG )
    value_r = Get_Value_From_Simulation( sim_coords_r, data_to_interpolate, field, sub_field, SG )
    value_interp = delta*value_r + (1-delta)*value_l 
    return value_interp
  
  value_l = Interpolate_3D( p0, p1, p2, data_to_interpolate, field, sub_field, SG, parameter_grid=parameter_grid, param_id=param_id-1, sim_coords_before=sim_coords_l, clip_params=clip_params )
  value_r = Interpolate_3D( p0, p1, p2, data_to_interpolate, field, sub_field, SG, parameter_grid=parameter_grid, param_id=param_id-1, sim_coords_before=sim_coords_r, clip_params=clip_params )
  value_interp = delta*value_r + (1-delta)*value_l
  return value_interp



def Interpolate_Comparable_1D( param_id, param_value,  comparable_grid, field, SG ):
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
  mean_l = comparable_grid[sim_id_l][field]['mean']
  mean_r = comparable_grid[sim_id_r][field]['mean'] 
  mean = (mean_r - mean_l ) / ( param_r - param_l ) * delta  + mean_l 
  return mean


  
def Interpolate_Observable_1D( param_id, observable, param_value,  SG, clip_values=True ):
  parameters = SG.parameters
  param_name = parameters[param_id]['name']
  param_vals = parameters[param_id]['values']
  param_max = max(param_vals)
  param_min = min(param_vals) 
  # if param_value < param_min or param_value > param_max:
  #   print( f'ERROR: Paramneter Value outside {param_name} {param_value} Range: [ {param_min} , {param_max} ] ')
  if clip_values:
    if param_value < param_min: param_value = param_min
    if param_value > param_max: param_value = param_max 
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

def Get_Comparable_Tau( factor_sigma_tau_becker=1, factor_sigma_tau_keating=1  ):
  comparable_z, comparable_tau, comparable_sigma = [], [], []

  # Add data Becker 2013
  data_set = data_optical_depth_Becker_2013
  z   = data_set['z']
  tau = data_set['tau']
  sigma = data_set['tau_sigma'] * factor_sigma_tau_becker
  indices = z < 4.3
  comparable_z.append(z[indices][::2])
  comparable_tau.append(tau[indices][::2])
  comparable_sigma.append(sigma[indices][::2])

  # Add data Jiani
  # data_set = data_optical_depth_Jiani
  # z   = data_set['z']
  # tau = data_set['tau']
  # sigma = data_set['tau_sigma'] * factor_sigma_tau
  # indices = z < 4.3
  # comparable_z.append(z[indices])
  # comparable_tau.append(tau[indices])
  # comparable_sigma.append(sigma[indices])
  
  # Add data Bosman
  # data_set = data_optical_depth_Bosman_2018
  # z   = data_set['z']
  # tau = data_set['tau']
  # sigma = data_set['tau_sigma'] * factor_sigma_tau
  # indices = z > 4.3
  # comparable_z.append(z[indices])
  # comparable_tau.append(tau[indices])
  # comparable_sigma.append(sigma[indices])

  # # Add data Keating 2020
  data_set = data_optical_depth_Keating_2020
  z   = data_set['z']
  tau = data_set['tau']
  sigma = data_set['tau_sigma'] * factor_sigma_tau_keating
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
  


