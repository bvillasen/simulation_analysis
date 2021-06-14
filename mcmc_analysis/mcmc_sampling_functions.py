import sys, os, time
import numpy as np
import h5py as h5
import pymc
analysis_dir = os.path.dirname(os.getcwd()) + '/'
sys.path.append(analysis_dir + 'tools')
from tools import *
from mcmc_data_functions import *


def Get_Highest_Likelihood_Params( param_samples, n_bins=100 ):
  param_ids = param_samples.keys()
  n_param = len(param_ids )
  param_samples_array = np.array([ param_samples[i]['trace'] for i in range(n_param) ] ).T

  hist_4D, bin_edges = np.histogramdd( param_samples_array, bins=n_bins )
  bin_centers = [ (edges[1:] + edges[:-1])/2 for edges in bin_edges ]
  hist_max = hist_4D.max()
  max_id = np.where( hist_4D == hist_max  )
  p_vals = np.array([ bin_centers[i][max_id[i]] for i in range(n_param) ])
  print( f"Highest_Likelihood: {hist_max} {p_vals}")
  while( len(p_vals.flatten()) > n_param ):
    n_bins = np.int( n_bins * 0.9 )
    hist_4D, bin_edges = np.histogramdd( param_samples_array, bins=n_bins )
    bin_centers = [ (edges[1:] + edges[:-1])/2 for edges in bin_edges ]
    hist_max = hist_4D.max()
    max_id = np.where( hist_4D == hist_max  )
    p_vals = np.array([ bin_centers[i][max_id[i]] for i in range(n_param) ])
    print( f"Highest_Likelihood: {hist_max} {p_vals}")
  return p_vals 


def Sample_Power_Spectrum_from_Trace( param_samples, data_grid, SG, hpi_sum=0.7, n_samples=None, params_HL=None ):  
  print(f'\nSampling Power Spectrum')
  ps_samples = {}
  param_ids = param_samples.keys()
  n_param = len(param_ids )
  param_samples_array = np.array([ param_samples[i]['trace'] for i in range(n_param) ] ).T
  if not n_samples: n_samples = len( param_samples[0]['trace'] )
  print(f' N Samples: {n_samples}')
  n_z_ids = len( data_grid.keys() )
  for id_z in range( n_z_ids ):
    print_line_flush( f' Sampling z_id: {id_z+1}/{n_z_ids}')
    ps_data = data_grid[id_z]
    ps_samples[id_z] = {}
    ps_samples[id_z]['z'] = ps_data['z']

    samples = []
    for i in range( n_samples ):
      p_vals = param_samples_array[i]
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
      v_l, v_r, v_max, sum = get_highest_probability_interval( bin_centers, distribution, fill_sum, log=True, n_interpolate=500)
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
    if params_HL is not None:
      if n_param == 3: ps_HL = Interpolate_3D(  params_HL[0], params_HL[1], params_HL[2], ps_data, 'P(k)', 'mean', SG, clip_params=True ) 
      if n_param == 4: ps_HL = Interpolate_4D(  params_HL[0], params_HL[1], params_HL[2], params_HL[3], ps_data, 'P(k)', 'mean', SG, clip_params=True ) 
      ps_samples[id_z]['Highest_Likelihood'] = ps_HL
  print( '\n' )
  return ps_samples
  



def Sample_Fields_from_Trace( fields_list, param_samples, data_grid, SG, hpi_sum=0.7, n_samples=None, params_HL=None, sample_log=False ):

  if sample_log: print( 'WARNING: Sampling Log Space')
  fields_list = [ field for field in fields_list if field != 'P(k)' ]
  print(f'\nSampling Fields: {fields_list}')
  param_ids = param_samples.keys()
  n_param = len(param_ids )
  param_samples_array = np.array([ param_samples[i]['trace'] for i in range(n_param) ] ).T

  if not n_samples: n_samples = len( param_samples[0]['trace'] )
  print(f' N Samples: {n_samples}')


  samples_out = {}
  for field in fields_list:
    print(f' Sampling Field: {field}')
    samples = []
    for i in range( n_samples ):
      p_vals = param_samples_array[i]
      if n_param == 3: interp = Interpolate_3D(  p_vals[0], p_vals[1], p_vals[2], data_grid, field, 'mean', SG, clip_params=True, interp_log=sample_log ) 
      if n_param == 4: interp = Interpolate_4D(  p_vals[0], p_vals[1], p_vals[2], p_vals[3], data_grid, field, 'mean', SG, clip_params=True, interp_log=sample_log ) 
      samples.append( interp )
    samples = np.array( samples ).T
    mean = np.array([ vals.mean() for vals in samples ])
    sigma = [ ]
    lower, higher = [], []
    for i in range( len( samples ) ):
      sigma.append( np.sqrt(  ( (samples[i] - mean[i])**2).mean()  ) )
      values = samples[i]
      n_bins = 1000
      distribution, bin_centers = compute_distribution( values, n_bins, log=False )
      fill_sum = hpi_sum
      log_hpi = True
      if field in  ['gamma', 'T0' ] : log_hpi = False
      if sample_log: log_hpi = False
      v_l, v_r, v_max, sum = get_highest_probability_interval( bin_centers, distribution, fill_sum, log=log_hpi, n_interpolate=1000)
      lower.append( v_l )
      higher.append( v_r )
    sigma  = np.array( sigma )
    lower  = np.array( lower )
    higher = np.array( higher )
    samples_stats = {}
    samples_stats['mean']   = mean
    # samples_stats['sigma']  = sigma
    samples_stats['z'] = data_grid[0][field]['z']
    samples_stats['lower']  = lower
    samples_stats['higher'] = higher
    if params_HL is not None:
      if n_param == 3: interp_HL = Interpolate_3D(  params_HL[0], params_HL[1], params_HL[2], data_grid, field, 'mean', SG, clip_params=True, interp_log=sample_log ) 
      if n_param == 4: interp_HL = Interpolate_4D(  params_HL[0], params_HL[1], params_HL[2], params_HL[3], data_grid, field, 'mean', SG, clip_params=True, interp_log=sample_log ) 
      samples_stats['Highest_Likelihood'] = interp_HL
    
    if sample_log:
      for key in samples_stats:
        if key == 'z': continue
        samples_stats[key] = 10**samples_stats[key]
    samples_out[field] = samples_stats
  return samples_out

