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

# data_sets = [ 'Boss', 'Walther', 'Boera', 'Viel' ]
data_ps_sets = [ 'Boss' ]
# data_sets = [ 'Walther' ]
# data_sets = [ 'Boera' ]
# data_sets = [ 'Boss', 'Walther' ]
# data_sets = [ 'Walther', 'Boera' ]
# data_sets = [ 'Walther', 'Viel' ]


name = ''
for data_set in data_ps_sets:
  name += data_set + '_'
name = name[:-1]  


field = 'P(k)+'

ps_data_dir = 'lya_statistics/data/'
mcmc_dir = root_dir + 'fit_mcmc/'
create_directory( mcmc_dir )
output_dir = mcmc_dir + f'fit_results_{field}_{name}/'
create_directory( output_dir )

# load_mcmc_results = False
load_mcmc_results = True


SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Load_Grid_Analysis_Data()
ps_range = SG.Get_Power_Spectrum_Range( kmax=0.01 )
sim_ids = SG.sim_ids

z_min = 2.0
z_max = 5.0 
ps_extras = { 'range':ps_range, 'data_dir':ps_data_dir, 'data_sets':data_ps_sets }
comparable_data = Get_Comparable_Composite( field,  z_min, z_max, ps_extras=ps_extras )
comparable_grid = Get_Comparable_Composite_from_Grid( field, comparable_data, SG )
# Plot_Comparable_Data( field, comparable_data, comparable_grid, output_dir  )

z_vals = [ 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.6, 5.0,  ]
data_grid, data_grid_power_spectrum = Get_Data_Grid_Composite( field, SG, z_vals=z_vals )


stats_file = output_dir + 'fit_mcmc.pkl'
samples_file = output_dir + 'samples_mcmc.pkl'


params = SG.parameters

if load_mcmc_results:
  print( f'Loading File: {stats_file}')
  stats = pickle.load( open( stats_file, 'rb' ) )
  param_stats = {}
  for p_id in params.keys():
    p_name = params[p_id]['name']
    p_stats = stats[p_name]
    params[p_id]['mean'] = p_stats['mean']
    params[p_id]['sigma'] = p_stats['standard deviation']
  print( f'Loading File: {samples_file}')
  param_samples = pickle.load( open( samples_file, 'rb' ) )

else:
  nIter = 200000 
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
  param_samples = Write_MCMC_Results( stats, MDL, params_mcmc,  stats_file, samples_file,  output_dir  )


# Make Corner plot from posteriors 
labels = { 'scale_He':r'$\beta_{\mathrm{He}}$', 'scale_H':r'$\beta_{\mathrm{H}}$', 'deltaZ_He':r'$\Delta z_{\mathrm{He}}$', 'deltaZ_H':r'$\Delta z_{\mathrm{H}}$'    }
Plot_Corner( param_samples, labels, output_dir  )


# Get the Highest_Likelihood parameter values 
params_HL = Get_Highest_Likelihood_Params( param_samples, n_bins=100 )

hpi_sum = 0.95
n_samples = 1000

# Obtain distribution of the power spectrum
samples_ps = Sample_Power_Spectrum_from_Trace( param_samples, data_grid_power_spectrum, SG, hpi_sum=hpi_sum, n_samples=n_samples, params_HL=params_HL )
Plot_Power_Spectrum_Sampling( samples_ps, ps_data_dir, output_dir, scales='large', system=system )
# 
# # Obtain distribution of the other fields 
# field_list = ['T0']
# samples_fields = Sample_Fields_from_Trace( field_list, param_samples, data_grid, SG, hpi_sum=hpi_sum, n_samples=n_samples, params_HL=params_HL )
# Plot_T0_Sampling( samples_fields['T0'], comparable_data, output_dir, system=system )


