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


data_name = 'fit_results_P(k)+tau_HeII_Boss'
# data_name = 'fit_results_P(k)+tau_HeII_Walther_kmax0.02_rescaleTauHeII1.0'
# data_name = 'fit_results_P(k)+tau_HeII_Walther_kmax0.10_rescaleTauHeII0.8'
# data_name = 'fit_results_P(k)+tau_HeII_Boss_Walther_kmax0.10_rescaleTauHeII0.3'

ps_data_dir = 'lya_statistics/data/'
mcmc_dir = root_dir + 'fit_mcmc/'
input_dir = mcmc_dir + f'{data_name}/' 
output_dir = input_dir + 'observable_samples/'
create_directory( output_dir )

kmax = 0.1

# sim_ids = range(10)
sim_ids = None
SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Load_Grid_UVB_Rates()
SG.Load_Grid_Analysis_Data( sim_ids=sim_ids )
sim_ids = SG.sim_ids
ps_range = SG.Get_Power_Spectrum_Range( kmax=kmax )



z_vals = [ 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.6   ]
data_grid, data_grid_power_spectrum = Get_Data_Grid_Composite(  ['P(k)', 'T0', 'tau', 'tau_HeII'], SG, z_vals=z_vals, sim_ids=sim_ids, load_uvb_rates=True )


stats_file = input_dir + 'fit_mcmc.pkl'
samples_file = input_dir + 'samples_mcmc.pkl'

params = SG.parameters

print( f'Loading File: {stats_file}')
stats = pickle.load( open( stats_file, 'rb' ) )
for p_id in params.keys():
  p_name = params[p_id]['name']
  p_stats = stats[p_name]
  params[p_id]['mean'] = p_stats['mean']
  params[p_id]['sigma'] = p_stats['standard deviation']
print( f'Loading File: {samples_file}')
param_samples = pickle.load( open( samples_file, 'rb' ) )



Write_Pickle_Directory( params, output_dir + 'params.pkl' )
Write_Pickle_Directory( stats, output_dir + 'fit_mcmc.pkl' )
Write_Pickle_Directory( param_samples, output_dir + 'samples_mcmc.pkl' )


n_bins = 100
for p_id in param_samples:
  trace = param_samples[p_id]['trace']
  distribution, bin_centers = compute_distribution( trace, n_bins )
  max_index = np.where( distribution == distribution.max() )[0][0]
  bin_max = bin_centers[max_index]
  param_samples[p_id]['max'] = bin_max
  print( f"{param_samples[p_id]['name']}   max {bin_max} " )

# rescale_width = [ 1.0, 1.0, 1.0, 1.0 ]
# for p_id in param_samples:
#   trace = param_samples[p_id]['trace']
#   max = param_samples[p_id]['max']
#   rescale = rescale_width[p_id]
#   trace_rescaled = ( trace - max )*rescale + max
#   param_samples[p_id]['trace'] = trace_rescaled
# 
# Get the Highest_Likelihood parameter values 
params_HL = Get_Highest_Likelihood_Params( param_samples, n_bins=100 )


# Make Corner plot from posteriors 
labels = { 'scale_He':r'$\beta_{\mathrm{He}}$', 'scale_H':r'$\beta_{\mathrm{H}}$', 'deltaZ_He':r'$\Delta z_{\mathrm{He}}$', 'deltaZ_H':r'$\Delta z_{\mathrm{H}}$'    }
data_labels = ''
Plot_Corner( param_samples, data_labels, labels, output_dir, n_bins_1D=40, n_bins_2D=40, lower_mask_factor=70  )


hpi_sum = 0.95
n_samples = 400000

load_samples = False


# # Obtain distribution of the power spectrum
# file_name = output_dir + 'samples_power_spectrum.pkl'
# if load_samples:
#   samples_ps = Load_Pickle_Directory( file_name )
# else:  
#   samples_ps = Sample_Power_Spectrum_from_Trace( param_samples, data_grid_power_spectrum, SG, hpi_sum=hpi_sum, n_samples=n_samples, params_HL=params_HL )
#   Write_Pickle_Directory( samples_ps, file_name )



# # Obtain distribution of the other fields
# file_name = output_dir + 'samples_fields.pkl' 
# field_list = ['T0', 'tau', 'tau_HeII']
# if load_samples:
#   samples_fields = Load_Pickle_Directory( file_name )
# else:  
#   samples_fields = Sample_Fields_from_Trace( field_list, param_samples, data_grid, SG, hpi_sum=hpi_sum, n_samples=n_samples, params_HL=params_HL )
#   Write_Pickle_Directory( samples_fields, file_name )


# Obtain distribution of the UVBRates
file_name = output_dir + 'samples_uvb_rates.pkl' 
field_list = ['photoheating_HI', 'photoheating_HeI', 'photoheating_HeII', 'photoionization_HI', 'photoionization_HeI', 'photoionization_HeII' ]
if load_samples:
  samples_uvb_rates = Load_Pickle_Directory( file_name )
else:  
  samples_uvb_rates = Sample_Fields_from_Trace( field_list, param_samples, data_grid, SG, hpi_sum=hpi_sum, n_samples=n_samples, params_HL=params_HL, sample_log=False )
  Write_Pickle_Directory( samples_uvb_rates, file_name )
