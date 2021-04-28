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
from mcmc_sampling_functions import *
from generate_grackle_uvb_file import Load_Grackle_File, Modify_UVB_Rates, Extend_Rates_Redshift, Copy_Grakle_UVB_Rates


data_name = 'fit_results_P(k)+tau_HeII_Boss'
# data_name = 'fit_results_P(k)+tau_HeII_Walther_kmax0.02_rescaleTauHeII1.0'
# data_name = 'fit_results_P(k)+tau_HeII_Walther_kmax0.10_rescaleTauHeII0.8'
# data_name = 'fit_results_P(k)+tau_HeII_Boss_Walther_kmax0.10_rescaleTauHeII0.3'

ps_data_dir = 'lya_statistics/data/'
mcmc_dir = root_dir + 'fit_mcmc/'
input_dir = mcmc_dir + f'{data_name}/' 
output_dir = input_dir + 'observable_samples/'
create_directory( output_dir )


SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
params = SG.parameters

stats_file = input_dir + 'fit_mcmc.pkl'
samples_file = input_dir + 'samples_mcmc.pkl'



print( f'Loading File: {stats_file}')
stats = pickle.load( open( stats_file, 'rb' ) )
for p_id in params.keys():
  p_name = params[p_id]['name']
  p_stats = stats[p_name]
  params[p_id]['mean'] = p_stats['mean']
  params[p_id]['sigma'] = p_stats['standard deviation']
print( f'Loading File: {samples_file}')
param_samples = pickle.load( open( samples_file, 'rb' ) )


# Get the Highest_Likelihood parameter values 
params_HL = Get_Highest_Likelihood_Params( param_samples, n_bins=100 )


hpi_sum = 0.95
n_samples = 100000

load_samples = False


grackle_file_name = 'rates_uvb/CloudyData_UVB_Puchwein2019_cloudy.h5'
grackle_data = Load_Grackle_File( grackle_file_name )
max_delta_z = 0.01
rates_data = Extend_Rates_Redshift( max_delta_z, grackle_data )


field_list = ['photoheating_HI', 'photoheating_HeI', 'photoheating_HeII', 'photoionization_HI', 'photoionization_HeI', 'photoionization_HeII' ]
rates_samples = {}
for field in field_list:
  rates_samples[field] = []

for sample_id in range(n_samples ):
  parameter_values = {}
  for p_id in param_samples:
    p_name = param_samples[p_id]['name']
    p_val = param_samples[p_id]['trace'][sample_id]
    parameter_values[p_name] = p_val
    # p_val = params_HL[p_id][0]
    parameter_values[p_name] = p_val
  # print( parameter_values )
  input_rates = Copy_Grakle_UVB_Rates( rates_data )
  rates_modified = Modify_UVB_Rates( parameter_values, input_rates )
  for field in field_list:
    rates_samples[field].append(  rates_modified[field] )

if params_HL is not None:
  parameter_values_HL = {}
  for p_id in param_samples:
    p_name = param_samples[p_id]['name']
    p_val = params_HL[p_id][0]
    parameter_values_HL[p_name] = p_val
  input_rates = Copy_Grakle_UVB_Rates( rates_data )  
  rates_modified_HL = Modify_UVB_Rates( parameter_values_HL, input_rates )


samples_uvb_rates = {}
for field in field_list:
  print( f'Obtaining Distribution: {field} ' )
  samples_stats = {}
  samples = rates_samples[field]
  samples = np.array( samples ).T
  mean = np.array([ vals.mean() for vals in samples ])
  sigma = [ ]
  lower, higher = [], []
  for i in range( len( samples ) ):
    sigma.append( np.sqrt(  ( (samples[i] - mean[i])**2).mean()  ) )
    values = samples[i]
    n_bins = 100
    distribution, bin_centers = compute_distribution( values, n_bins, log=False )
    fill_sum = hpi_sum
    log_hpi = True
    v_l, v_r, v_max, sum = get_highest_probability_interval( bin_centers, distribution, fill_sum, log=log_hpi, n_interpolate=1000)
    lower.append( v_l )
    higher.append( v_r )
  sigma  = np.array( sigma )
  lower  = np.array( lower )
  higher = np.array( higher )
  samples_stats = {}
  samples_stats['z'] = rates_modified['z']
  samples_stats['mean']   = mean
  samples_stats['sigma']  = sigma
  samples_stats['lower']  = lower
  samples_stats['higher'] = higher
  samples_stats['Highest_Likelihood'] = rates_modified_HL[field]  
  samples_uvb_rates[field] = samples_stats

    
file_name = output_dir + 'samples_uvb_rates_new.pkl' 
Write_Pickle_Directory( samples_uvb_rates, file_name )



