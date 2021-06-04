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
from mcmc_sampling_functions import *


ps_data_dir = 'lya_statistics/data/'
mcmc_dir = root_dir + 'fit_mcmc/'
output_dir = mcmc_dir + f'observable_figures/'
create_directory( output_dir )
rescale_walter_file = ps_data_dir + 'rescale_walther_to_boss.pkl' 

data_boss = 'fit_results_P(k)+tau_HeII_Boss'
data_boss_irsic = 'fit_results_P(k)+tau_HeII_Boss_Irsic'
data_boss_boera = 'fit_results_P(k)+tau_HeII_Boss_Boera'
data_boss_irsic_boera = 'fit_results_P(k)+tau_HeII_Boss_Irsic_Boera'

data_sets = [ data_boss, data_boss_irsic, data_boss_boera, data_boss_irsic_boera ]
data_labels = [ 'BOSS', 'BOSS + Irsic', 'BOSS + Boera', 'BOSS + Irsic + Boera'  ]

samples_all = {}
samples_all['param'] = {}
samples_all['P(k)'] = {}
field_list = ['T0', 'tau', 'tau_HeII']
for field in field_list:
  samples_all[field] = {}


for data_id, data_name in enumerate(data_sets):
  
  print(f'Loading Dataset: {data_name}' )
  input_dir = mcmc_dir + f'{data_name}/observable_samples/' 
  stats_file = input_dir + 'fit_mcmc.pkl'
  samples_file = input_dir + 'samples_mcmc.pkl'

  params = Load_Pickle_Directory( input_dir + 'params.pkl' )

  print( f'Loading File: {stats_file}')
  stats = pickle.load( open( stats_file, 'rb' ) )
  for p_id in params.keys():
    p_name = params[p_id]['name']
    p_stats = stats[p_name]
    params[p_id]['mean'] = p_stats['mean']
    params[p_id]['sigma'] = p_stats['standard deviation']
  print( f'Loading File: {samples_file}')
  param_samples = pickle.load( open( samples_file, 'rb' ) )
  samples_all['param'][data_id] = param_samples

  # Get the Highest_Likelihood parameter values 
  params_HL = Get_Highest_Likelihood_Params( param_samples, n_bins=100 )
  # 
  # Obtain distribution of the power spectrum
  file_name = input_dir + 'samples_power_spectrum.pkl'
  samples_ps = Load_Pickle_Directory( file_name )
  
  # Obtain distribution of the other fields
  file_name = input_dir + 'samples_fields.pkl'
  samples_fields = Load_Pickle_Directory( file_name )
  
  samples_all['P(k)'][data_id] = samples_ps
  for field in field_list:
    samples_all[field][data_id] = samples_fields[field] 


# corner_labels = { 'scale_He':r'$\beta_{\mathrm{He}}$', 'scale_H':r'$\beta_{\mathrm{H}}$', 'deltaZ_He':r'$\Delta z_{\mathrm{He}}$', 'deltaZ_H':r'$\Delta z_{\mathrm{H}}$'    }
# Plot_Corner( samples_all['param'], data_labels, corner_labels, output_dir, n_bins_1D=40, n_bins_2D=40, lower_mask_factor=500, multiple=True  )


# Plot_Power_Spectrum_Sampling( samples_all['P(k)'], ps_data_dir, output_dir, scales='large', linewidth=2, system=system, label=data_labels,  multiple=True )
# Plot_Power_Spectrum_Sampling( samples_all['P(k)'], ps_data_dir, output_dir, scales='middle', linewidth=2, system=system, label=data_labels, rescaled_walther=True, rescale_walter_file=rescale_walter_file, multiple=True )
# Plot_Power_Spectrum_Sampling( samples_all['P(k)'], ps_data_dir, output_dir, scales='all', linewidth=2, system=system, label=data_labels, rescaled_walther=True, rescale_walter_file=rescale_walter_file, multiple=True )
# 
Plot_T0_Sampling( samples_all['T0'], output_dir, system=system, label=data_labels, plot_splines=True, multiple=True)

# 
# Plot_tau_HeII_Sampling( samples_all['tau'], samples_all['tau_HeII'], output_dir, system=system, label=data_labels, multiple=True )