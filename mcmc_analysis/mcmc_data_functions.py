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
from data_optical_depth_HeII import data_tau_HeII_Worserc_2019
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

def Get_Data_Grid_Composite( fields_list,  SG, z_vals=None, load_normalized_ps=False, sim_ids=None, load_uvb_rates=False ):
  # fields_list = fields.split('+')
  data_grid_all = {}
  for field in fields_list:
    if field == 'T0':   data_grid_all[field] = Get_Data_Grid( [field], SG, sim_ids=sim_ids ) 
    if field == 'tau':  data_grid_all[field] = Get_Data_Grid( [field], SG, sim_ids=sim_ids ) 
    if field == 'P(k)': data_grid_all[field] = Get_Data_Grid_Power_spectrum( z_vals, SG, normalized_ps=load_normalized_ps, sim_ids=sim_ids )
    if field == 'tau_HeII':  data_grid_all[field] = Get_Data_Grid( [field], SG, sim_ids=sim_ids ) 

  data_grid = {}
  if not sim_ids: sim_ids = SG.sim_ids
  
  for sim_id in sim_ids:  
    data_grid[sim_id] = {}
    for field in fields_list:
     if field == 'P(k)' or field == '': continue
     z = data_grid_all[field][sim_id]['z']
     mean = data_grid_all[field][sim_id][field]['mean']
     data_grid[sim_id][field] = {'mean':mean, 'z':z }

  if load_uvb_rates:
    photoheating_keys    = [ 'piHI', 'piHeI', 'piHeII' ]
    photoionization_keys = [ 'k24', 'k25', 'k26' ]

    key_names = { 'piHI': 'photoheating_HI',   'piHeI': 'photoheating_HeI',  'piHeII': 'photoheating_HeII', 
    'k24': 'photoionization_HI', 'k26': 'photoionization_HeI', 'k25': 'photoionization_HeII' }
     
    for sim_id in sim_ids:
      for key in photoheating_keys:
        z = SG.Grid[sim_id]['UVB_rates']['z']
        rates = SG.Grid[sim_id]['UVB_rates']['Photoheating'][key]
        key_name = key_names[key]
        data_grid[sim_id][key_name] = {}
        data_grid[sim_id][key_name]['mean'] = rates
        data_grid[sim_id][key_name]['z'] = z
        
      for key in photoionization_keys:
        z = SG.Grid[sim_id]['UVB_rates']['z']
        rates = SG.Grid[sim_id]['UVB_rates']['Chemistry'][key]
        key_name = key_names[key]
        data_grid[sim_id][key_name] = {}
        data_grid[sim_id][key_name]['mean'] = rates
        data_grid[sim_id][key_name]['z'] = z
    
  if 'P(k)' in fields_list:  return data_grid, data_grid_all['P(k)']
  else:                      return data_grid



def Get_Data_Grid_Power_spectrum( z_vals, SG, normalized_ps=False, sim_ids=None ):
  data_grid = {}
  print_norm = True
  for id_z, z_val in enumerate( z_vals ):
    data_grid[id_z] = {}
    data_grid[id_z]['z'] = z_val
    if not sim_ids: sim_ids = SG.sim_ids
    for sim_id in sim_ids:
      if normalized_ps:
        sim_ps_data = SG.Grid[sim_id]['analysis']['power_spectrum_normalized']
        label = sim_ps_data['normalization_key']
        if print_norm: print( f'Loading Normalized: {label}' )
        print_norm = False 
      else: sim_ps_data = SG.Grid[sim_id]['analysis']['power_spectrum']
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


def Get_Data_Grid( fields, SG, sim_ids=None ):
  if not sim_ids: sim_ids = SG.sim_ids
  data_grid = {}
  for sim_id in sim_ids:
    data_grid[sim_id] = {}
    data_grid[sim_id]['z'] = SG.Grid[sim_id]['analysis']['z']
    for field in fields:
      data_grid[sim_id][field] = {}
      data_grid[sim_id][field]['mean'] = SG.Grid[sim_id]['analysis'][field]
  return data_grid

def Get_Comparable_Composite_from_Grid( fields, comparable_data, SG, log_ps=False, load_normalized_ps=False ):
  fields_list = fields.split('+')

  sim_ids = SG.sim_ids
  comparable_grid_all = {}
  for field in fields_list:
    if field == 'T0':   comparable_grid_all[field] = Get_Comparable_T0_from_Grid(  comparable_data[field], SG )
    if field == 'tau':  comparable_grid_all[field] = Get_Comparable_tau_from_Grid( comparable_data[field], SG )
    if field == 'P(k)': comparable_grid_all[field] = Get_Comparable_Power_Spectrum_from_Grid( comparable_data[field]['separate'], SG, log_ps=log_ps, normalized_ps=load_normalized_ps )
    if field == 'tau_HeII':  comparable_grid_all[field] = Get_Comparable_tau_HeII_from_Grid( comparable_data[field], SG )
    
  comparable_grid = {}
  for sim_id in sim_ids:
    comparable_grid[sim_id] = {}
    mean_all = []
    for field in fields_list:
      if field == '': continue
      comparable_grid[sim_id][field] = comparable_grid_all[field][sim_id]
      mean_all.append( comparable_grid_all[field][sim_id]['mean'] )
    comparable_grid[sim_id][fields] = {'mean': np.concatenate( mean_all ) }
  return comparable_grid


def Get_Comparable_Power_Spectrum_from_Grid( comparable_data, SG, log_ps=False, normalized_ps=False ):
  
  print( 'Generating Simulation P(k) comparable data:')
  indices = comparable_data.keys()
  comparable_grid = {}
  sim_ids = SG.sim_ids
  print_norm = True
  for sim_id in sim_ids:
    comparable_grid[sim_id] = {}  
    if normalized_ps:
      sim_data = SG.Grid[sim_id]['analysis']['power_spectrum_normalized']
      if print_norm:
        label = sim_data['normalization_key']
        print( f'Loading Normalized PS: {label}' )
        print_norm = False
    else:sim_data = SG.Grid[sim_id]['analysis']['power_spectrum']
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
      if log_ps: sim_delta = np.log( sim_delta ) 
      sim_delta_interp = np.interp( data_kvals, sim_kvals, sim_delta )
      diff = ( sim_delta_interp - data_delta_vals ) / data_delta_vals
      sim_delta_all.append( sim_delta_interp )
      
    comparable_grid[sim_id]['mean'] = np.concatenate( sim_delta_all )
  n_points = len(  comparable_grid[0]['mean'] )
  print ( f' N sim points: {n_points}' )
  return comparable_grid  

def Get_Comparable_T0_from_Grid( comparable_data, SG ):
  return Get_Comparable_Field_from_Grid( 'T0', comparable_data, SG )

def Get_Comparable_tau_from_Grid( comparable_data, SG ):
  return Get_Comparable_Field_from_Grid( 'tau', comparable_data, SG )
  
def Get_Comparable_tau_HeII_from_Grid( comparable_data, SG ):
  return Get_Comparable_Field_from_Grid( 'tau_HeII', comparable_data, SG )

def Get_Comparable_Field_from_Grid( field, comparable_data, SG, interpolate=True ):
  print( f' Loading Comparabe from Grid: {field}')
  z_data = comparable_data['z']
  sim_ids = SG.sim_ids
  comparable_grid = {}
  for sim_id in sim_ids:
    comparable_grid[sim_id] = {}
    sim_analysis = SG.Grid[sim_id]['analysis']
    z_sim = sim_analysis['z']
    if interpolate:
      mean_sim = sim_analysis[field]
      if z_sim[0] > z_sim[-1]:  
        z_sim_sorted = z_sim[::-1]
        mean_sim = mean_sim[::-1]
      else: z_sim_sorted = z_sim
      mean_interp = np.interp( z_data, z_sim_sorted, mean_sim )
      comparable_grid[sim_id]['z'] = z_data
      comparable_grid[sim_id]['mean'] = mean_interp    
    else:  
      indices = []
      for z in z_data:
        diff = np.abs( z_sim - z )
        diff_min = diff.min()
        if diff_min > 0.03:
          print( f'Warning Z diff is large: diff:{diff_min}')
        indx = np.where(diff == diff_min)[0]
        indices.append( indx )
      indices = np.array( indices )
      comparable_grid[sim_id]['z'] = sim_analysis['z'][indices].flatten()
      comparable_grid[sim_id]['mean'] = sim_analysis[field][indices].flatten()    
  return comparable_grid
  

def Get_Comparable_Composite( fields, z_min, z_max, ps_extras=None, tau_extras=None, log_ps=False, rescale_tau_HeII_sigma=1.0 ):
  
  rescaled_walther = False
  rescale_walter_file = None
  if ps_extras is not None:
    ps_data_dir = ps_extras['data_dir']
    data_ps_sets = ps_extras['data_sets'] 
    ps_range = ps_extras['range']
    rescaled_walther = ps_extras['rescaled_walther'] 
    rescale_walter_file = ps_extras['rescale_walter_file']

  if tau_extras is not None:
    factor_sigma_tau_becker  = tau_extras['factor_sigma_becker']
    factor_sigma_tau_keating = tau_extras['factor_sigma_keating']
    
  fields_list = fields.split('+')
  mean_all, sigma_all = [], []
  comparable_all = {}
  for field in fields_list:
    append_comparable = False
    if field == 'P(k)':
      comparable_ps = Get_Comparable_Power_Spectrum( ps_data_dir, z_min, z_max, data_ps_sets, ps_range, log_ps=log_ps, rescaled_walther=rescaled_walther, rescale_walter_file=rescale_walter_file )
      comparable_ps_all = comparable_ps['P(k)']
      comparable_ps_separate = comparable_ps['separate']
      comparable_all['P(k)'] = { 'all':comparable_ps_all, 'separate':comparable_ps_separate }
      print('Added comparable P(k) separate')
      mean_all.append( comparable_ps_all['mean'] )
      sigma_all.append( comparable_ps_all['sigma'] )
    if field == 'T0':  
      comparable_field = Get_Comparable_T0_Gaikwad()
      append_comparable = True
    if field == 'tau': 
      comparable_field = Get_Comparable_Tau( z_min, z_max, factor_sigma_tau_becker=factor_sigma_tau_becker, factor_sigma_tau_keating=factor_sigma_tau_keating )
      append_comparable = True
    if field == 'tau_HeII': 
      comparable_field = Get_Comparable_Tau_HeII( rescale_tau_HeII_sigma=rescale_tau_HeII_sigma )
      append_comparable = True
    if append_comparable:
      comparable_all[field] = comparable_field
      mean_all.append( comparable_field['mean'] )
      sigma_all.append( comparable_field['sigma'] )
  comparable_all[fields] = { 'mean':np.concatenate(mean_all), 'sigma':np.concatenate(sigma_all) }  
  return comparable_all
  


def Get_Comparable_Power_Spectrum( ps_data_dir, z_min, z_max, data_sets, ps_range, log_ps=False, rescaled_walther=False, rescale_walter_file=None ):
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
  log_data_ps, log_data_ps_sigma = [], []
  sim_z, sim_kmin, sim_kmax = ps_range['z'], ps_range['k_min'], ps_range['k_max']

  ps_data = {}
  data_id = 0
  for data_index, data_name in enumerate(data_sets):
    print( f' Loading P(k) Data: {data_name}' )
    data_set = data_dir[data_name]
    keys = data_set.keys()
    n_indices = len(keys) - 1
    if data_name == 'Walther':
      if rescaled_walther:
        print(f" Loading Walther rescale values: {rescale_walter_file}")
        rescale_walter_alphas = Load_Pickle_Directory( rescale_walter_file)
    for index in range(n_indices):
      data = data_set[index]
      z = data['z']
      if z >= z_min and z <= z_max:
        diff = np.abs( sim_z - z )
        id_min = np.where( diff == diff.min() )[0][0]
        z_sim = sim_z[id_min]
        kmin = sim_kmin[id_min]
        kmax = sim_kmax[id_min]
        k_vals = data['k_vals']
        k_indices = np.where( (k_vals >= kmin) & (k_vals <= kmax) )
        k_vals = k_vals[k_indices]
        delta_ps = data['delta_power'][k_indices]
        delta_ps_sigma = data['delta_power_error'][k_indices]
        log_delta_ps = np.log( delta_ps )
        log_delta_ps_sigma = 1/delta_ps * delta_ps_sigma
        if data_name == 'Walther' and rescaled_walther:
          rescale_z = rescale_walter_alphas[index]['z']
          rescale_alpha = rescale_walter_alphas[index]['alpha']
          print( f'  Rescaling z={rescale_z:.1f}    alpha={rescale_alpha:.3f} ')
          delta_ps *= rescale_alpha
        ps_data[data_id] = {'z':z, 'k_vals':k_vals, 'delta_ps':delta_ps, 'delta_ps_sigma':delta_ps_sigma }
        data_z.append( z )
        data_kvals.append( k_vals )
        data_ps.append( delta_ps )
        data_ps_sigma.append( delta_ps_sigma )
        log_data_ps.append( log_delta_ps )
        log_data_ps_sigma.append( log_delta_ps_sigma )
        data_id += 1
  k_vals_all         = np.concatenate( data_kvals )
  delta_ps_all       = np.concatenate( data_ps )
  delta_ps_sigma_all = np.concatenate( data_ps_sigma )
  log_delta_ps_all       = np.concatenate( log_data_ps )
  log_delta_ps_sigma_all = np.concatenate( log_data_ps_sigma )
  ps_data_out = {'P(k)':{}, 'separate':ps_data }
  ps_data_out['P(k)']['k_vals'] = k_vals_all
  if log_ps:
    ps_data_out['P(k)']['mean']   = log_delta_ps_all
    ps_data_out['P(k)']['sigma']  = log_delta_ps_sigma_all
  else:
    ps_data_out['P(k)']['mean']   = delta_ps_all
    ps_data_out['P(k)']['sigma']  = delta_ps_sigma_all

  n_data_points = len( k_vals_all )
  print( f' N data points: {n_data_points}' )
  return ps_data_out



def Get_Comparable_T0_Gaikwad():
  print( 'Loading T0 Data: ')
  data_sets = [ data_thermal_history_Gaikwad_2020b, data_thermal_history_Gaikwad_2020a ]
  z = np.concatenate( [ds['z'] for ds in data_sets ])
  data_mean  = np.concatenate( [ds['T0'] for ds in data_sets ])
  data_sigma = np.concatenate( [( ds['T0_sigma_plus'] + ds['T0_sigma_minus'] )*0.5  for ds in data_sets ] )
  comparable = {}
  comparable['z'] = z
  comparable['mean'] = data_mean
  comparable['sigma'] = data_sigma
  print( f' N data points: {len(data_mean)} ' )
  return comparable

def Get_Comparable_Tau_HeII( rescale_tau_HeII_sigma = 1.0 ):
  comparable_z, comparable_tau, comparable_sigma = [], [], []
  data_set = data_tau_HeII_Worserc_2019
  z   = data_set['z']
  tau = data_set['tau']
  sigma = data_set['tau_sigma'] 
  if rescale_tau_HeII_sigma != 1.0:
    print( f' Rescaling tau HeII sigma by {rescale_tau_HeII_sigma} ')
    sigma *= rescale_tau_HeII_sigma 
  comparable = {}
  comparable['z']     = z
  comparable['mean']  = tau
  comparable['sigma'] = sigma
  return comparable

  
def Get_Comparable_Tau( z_min, z_max, factor_sigma_tau_becker=1, factor_sigma_tau_keating=1  ):
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
  
  z     = np.concatenate( comparable_z )
  mean  = np.concatenate( comparable_tau )
  sigma = np.concatenate( comparable_sigma )
  
  indices = ( z >= z_min ) * ( z <= z_max )
  comparable = {}
  comparable['z']     = z[indices]
  comparable['mean']  = mean[indices]
  comparable['sigma'] = sigma[indices]
  return comparable
  

def Interpolate_Power_Spectrum( p_vals, data_grid, SG ):
  ps_output = {}
  n_param = len( p_vals )
  n_z_ids = len( data_grid.keys() )
  for id_z in range( n_z_ids ):
    ps_data = data_grid[id_z]
    ps_output[id_z] = {}
    ps_output[id_z]['z'] = ps_data['z']
    ps_output[id_z]['k_vals'] = ps_data[id_z]['P(k)']['k_vals']

    if n_param == 3: ps_interp = Interpolate_3D(  p_vals[0], p_vals[1], p_vals[2], ps_data, 'P(k)', 'mean', SG, clip_params=True ) 
    if n_param == 4: ps_interp = Interpolate_4D(  p_vals[0], p_vals[1], p_vals[2], p_vals[3], ps_data, 'P(k)', 'mean', SG, clip_params=True ) 
    ps_output[id_z]['mean'] = ps_interp
  return ps_output


def Interpolate_4D( p0, p1, p2, p3, data_to_interpolate, field, sub_field, SG, clip_params=False, parameter_grid=None, param_id=None, sim_coords_before=None, interp_log=False ):
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
    value_l = Get_Value_From_Simulation( sim_coords_l, data_to_interpolate, field, sub_field, SG, interp_log=interp_log )
    value_r = Get_Value_From_Simulation( sim_coords_r, data_to_interpolate, field, sub_field, SG, interp_log=interp_log )
    value_interp = delta*value_r + (1-delta)*value_l 
    return value_interp
  
  value_l = Interpolate_4D( p0, p1, p2, p3, data_to_interpolate, field, sub_field, SG, parameter_grid=parameter_grid, param_id=param_id-1, sim_coords_before=sim_coords_l, clip_params=clip_params, interp_log=interp_log )
  value_r = Interpolate_4D( p0, p1, p2, p3, data_to_interpolate, field, sub_field, SG, parameter_grid=parameter_grid, param_id=param_id-1, sim_coords_before=sim_coords_r, clip_params=clip_params, interp_log=interp_log )
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





def Interpolate_1D( param_id, param_value,  comparable_grid, field, SG ):
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


def Get_Value_From_Simulation( sim_coords, data_to_interpolate, field, sub_field, SG, interp_log=False ):
  sim_id = Get_Simulation_ID_From_Coordinates( sim_coords, SG )
  sim = SG.Grid[sim_id]
  param_values = sim['parameter_values']
  value = data_to_interpolate[sim_id][field][sub_field]
  if interp_log: value = np.log10(value)
  return value
