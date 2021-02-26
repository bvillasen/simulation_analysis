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



SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Load_Grid_Analysis_Data()
ps_range = SG.Get_Power_Spectrum_Range( kmax=0.01 )
sim_ids = SG.sim_ids
parameters = SG.parameters

z_vals_all = [ 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.6,   ]
n_zvals = len( z_vals_all )
data_grid = Get_Data_Grid_Power_spectrum( z_vals_all, SG )

params_fit = {}

z_indx = 0
for z_indx in range( n_zvals ):
  input_dir = root_dir + f'fit_results_{field}_{name}_single_{z_indx}/'

  z_val = z_vals_all[ z_indx ]
  z_vals = [ z_val ]

  stats_file = input_dir + 'fit_mcmc.pkl'

  fields = [ 'P(k)' ]

  params = {}

  print( f'Loading File: {stats_file}')
  stats = pickle.load( open( stats_file, 'rb' ) )
  param_stats = {}
  for p_id in parameters.keys():
    p_name = parameters[p_id]['name']
    p_stats = stats[p_name]
    params[p_id] = {}
    params[p_id]['mean'] = p_stats['mean']
    params[p_id]['sigma'] = p_stats['standard deviation']
    params[p_id]['quantiles'] = p_stats['quantiles']
    mean = params[p_id]['mean']

  params_fit[z_indx] = {}
  params_fit[z_indx]['z'] = z_val
  params_fit[z_indx]['params'] = params
  
  
  
n_samples = 1000
ps_samples = Sample_Power_Spectrum_Multiple_Params( n_samples, params_fit, data_grid, SG, hpi_sum=0.98  )
Plot_Power_Spectrum_Sampling( ps_samples, ps_data_dir, output_dir, scales='large', system=system )

  
  
  
  
  
# 
# z_vals, p_vals = [], {}
# for p_id in parameters.keys():
#   p_vals[p_id] = { 'low':[], 'high':[], 'mean':[] }
# for z_indx in range( n_zvals ):
#   z_vals.append( params_fit[z_indx]['z'] )  
#   for p_id in params.keys():
#     p_vals[p_id]['low'].append( params_fit[z_indx]['params'][p_id]['quantiles'][2.5]  )
#     p_vals[p_id]['high'].append( params_fit[z_indx]['params'][p_id]['quantiles'][97.5]  )
#     p_vals[p_id]['mean'].append( params_fit[z_indx]['params'][p_id]['mean']  )
# 
# 
# nrows, ncols = 3, 1
# fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20*ncols,5*nrows))
# 
# import matplotlib
# matplotlib.rcParams['mathtext.fontset'] = 'cm'
# matplotlib.rcParams['mathtext.rm'] = 'serif'
# 
# 
# fsize = 20
# alpha = 0.4
# 
# p_id = 0
# ax = ax_l[p_id]
# low  = p_vals[p_id]['low']
# high = p_vals[p_id]['high']
# mean = p_vals[p_id]['mean']
# ax.plot( z_vals, mean )
# ax.fill_between( z_vals, high, low, alpha= alpha)
# ax.set_xlabel( r'$z$', fontsize=fsize )
# ax.set_ylabel( r'$\beta_{\mathrm{He}}$', fontsize=fsize)
# 
# p_id = 1
# ax = ax_l[p_id]
# low  = p_vals[p_id]['low']
# high = p_vals[p_id]['high']
# mean = p_vals[p_id]['mean']
# ax.plot( z_vals, mean )
# ax.fill_between( z_vals, high, low, alpha= alpha)
# ax.set_xlabel( r'$z$', fontsize=fsize )
# ax.set_ylabel( r'$\beta_{\mathrm{H}}$', fontsize=fsize)
# 
# p_id = 2
# ax = ax_l[p_id]
# low  = p_vals[p_id]['low']
# high = p_vals[p_id]['high']
# mean = p_vals[p_id]['mean']
# ax.plot( z_vals, mean )
# ax.fill_between( z_vals, high, low, alpha= alpha)
# ax.set_xlabel( r'$z$', fontsize=fsize )
# ax.set_ylabel( r'$\Delta z_{\mathrm{He}}$', fontsize=fsize)
# 
# 
# figure_name = output_dir + 'params_fit_single.png'
# fig.savefig( figure_name, bbox_inches='tight', dpi=500 )
# print( f'Saved Figure: {figure_name}' )
# 
# 
# # Plot
# 
