import os, sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from parameters_UVB_rates import param_UVB_Rates
from simulation_grid import Simulation_Grid
from simulation_parameters import *
from mcmc_functions import *
from mcmc_data_functions import *
from data_thermal_history import *
from mcmc_plotting_functions import *

# data_sets = [ 'Boss', 'Walther', 'Boera', 'Viel' ]
data_sets = [ 'Boss' ]
# data_sets = [ 'Walther' ]
# data_sets = [ 'Boera' ]
# data_sets = [ 'Boss', 'Walther' ]
# data_sets = [ 'Walther', 'Boera' ]
# data_sets = [ 'Walther', 'Viel' ]


name = ''
for data_set in data_sets:
  name += data_set + '_'
name = name[:-1] 

field = 'P(k)'
ps_data_dir = 'lya_statistics/data/'
output_dir = root_dir + f'fit_results_{field}_{name}/'
create_directory( output_dir )

load_mcmc_stats = False
# load_mcmc_stats = True


SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Load_Grid_Analysis_Data()
ps_range = SG.Get_Power_Spectrum_Range( kmax=0.01 )
sim_ids = SG.sim_ids

z_min = 2.5
z_max = 5.0 
comparable_data = Get_Comparable_Power_Spectrum(  ps_data_dir, z_min, z_max, data_sets, ps_range )
comparable_grid = Get_Comparable_Power_Spectrum_from_Grid( comparable_data['separate'], SG )


stats_file = output_dir + 'fit_mcmc.pkl'

fields = [ 'P(k)' ]
z_vals = [ 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.6, 5.0,  ]
data_grid = Get_Data_Grid_Power_spectrum( z_vals, SG )

params = SG.parameters

if load_mcmc_stats:
  print( f'Loading File: {stats_file}')
  stats = pickle.load( open( stats_file, 'rb' ) )
  param_stats = {}
  for p_id in params.keys():
    p_name = params[p_id]['name']
    p_stats = stats[p_name]
    params[p_id]['mean'] = p_stats['mean']
    params[p_id]['sigma'] = p_stats['standard deviation']

else:
  nIter = 200000
  # nIter = 100000
  nBurn = nIter / 5
  nThin = 1
  # model, params_mcmc = mcmc_model_3D( comparable_data, comparable_grid, field, 'mean', SG )
  model, params_mcmc = mcmc_model_4D( comparable_data, comparable_grid, field, 'mean', SG )
  MDL = pymc.MCMC( model )  
  MDL.sample( iter=nIter, burn=nBurn, thin=nThin )
  stats = MDL.stats()
  param_stats = {}
  for p_id in params.keys():
    p_name = params[p_id]['name']
    p_stats = stats[p_name]
    params[p_id]['mean'] = p_stats['mean']
    params[p_id]['sigma'] = p_stats['standard deviation']
  Plot_MCMC_Stats( stats, MDL, params_mcmc,  stats_file, output_dir, plot_corner=False )

  samples = {} 
  for p_id in params_mcmc.keys():
    param = params_mcmc[p_id]
    samples[p_id] = {}
    samples[p_id]['name'] = param['name']
    samples[p_id]['trace'] = param['sampler'].trace() 
  
  # labels = { 'scale_He':r'$\beta_{\mathrm{He}}$', 'scale_H':r'$\beta_{\mathrm{H}}$', 'deltaZ_He':r'$\Delta z_{\mathrm{He}}$'    }
  labels = { 'scale_He':r'$\beta_{\mathrm{He}}$', 'scale_H':r'$\beta_{\mathrm{H}}$', 'deltaZ_He':r'$\Delta z_{\mathrm{He}}$', 'deltaZ_H':r'$\Delta z_{\mathrm{H}}$'    }
  Plot_Corner( samples, labels, output_dir  )




for p_id in params.keys():
  p_name = params[p_id]['name']
  p_stats = stats[p_name]
  params[p_id]['mean'] = p_stats['mean']
  params[p_id]['sigma'] = p_stats['standard deviation']

# 
# # Obtain distribution of the power spectrum
# n_samples = 1000
# ps_samples = Sample_Power_Spectrum( n_samples, params, data_grid, SG  )
# # chi2_vals = Get_Chi2( observables, params, comparable_grid, comparable_data, SG )
# # Plot_Observables( observables_samples, comparable_data, params, SG, 'sampling', output_dir, chi2=chi2_vals)
# # Plot_Observables( observables_samples, comparable_data, params, SG, 'grid', output_dir, chi2=None)
# # 

# Obtain distribution of the power spectrum
n_samples = 1000
ps_samples = Sample_Power_Spectrum( n_samples, params, data_grid, SG, hpi_sum=0.98  )
Plot_Power_Spectrum_Sampling( ps_samples, ps_data_dir, output_dir, scales='large', system=system )