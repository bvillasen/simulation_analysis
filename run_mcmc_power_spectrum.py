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

# field = 'T0'
# field = 'tau'
# field = 'T0+tau'
field = 'P(k)'
output_dir = root_dir + f'fit_results_{field}/'
create_directory( output_dir )

load_mcmc_stats = False


SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Load_Grid_Analysis_Data()
sim_ids = SG.sim_ids

# k_min = SG.Grid[0]['analysis']['power_spectrum']['k_min']
# for sim_id in sim_ids:
#   comp = SG.Grid[sim_id]['analysis']['power_spectrum']['k_min'] - k_min
#   print( comp)

# comparable_data = Get_Comparable_Composite_T0_tau( factor_sigma_tau_becker=6.0, factor_sigma_tau_keating=4.0, )
# comparable_grid = Get_Comparable_Composite_T0_tau_from_Grid( comparable_data, SG )
# 
# 
# stats_file = output_dir + 'fit_mcmc.pkl'
# 
# fields = [ 'T0', 'tau' ]
# data_grid = Get_Data_Grid( fields, SG )
# 
# params = SG.parameters
# 
# if load_mcmc_stats:
#   print( f'Loading File: {stats_file}')
#   stats = pickle.load( open( stats_file, 'rb' ) )
#   param_stats = {}
#   for p_id in params.keys():
#     p_name = params[p_id]['name']
#     p_stats = stats[p_name]
#     params[p_id]['mean'] = p_stats['mean']
#     params[p_id]['sigma'] = p_stats['standard deviation']
# 
# else:
#   nIter = 200000
#   # nIter = 100000
#   nBurn = nIter / 5
#   nThin = 1
#   model, params_mcmc = mcmc_model_3D( comparable_data, comparable_grid, field, 'mean', SG )
#   MDL = pymc.MCMC( model )  
#   MDL.sample( iter=nIter, burn=nBurn, thin=nThin )
#   stats = MDL.stats()
#   param_stats = {}
#   for p_id in params.keys():
#     p_name = params[p_id]['name']
#     p_stats = stats[p_name]
#     params[p_id]['mean'] = p_stats['mean']
#     params[p_id]['sigma'] = p_stats['standard deviation']
#   Plot_MCMC_Stats( stats, MDL, params_mcmc,  stats_file, output_dir, plot_corner=False )
# 
#   samples = {} 
#   for p_id in params_mcmc.keys():
#     param = params_mcmc[p_id]
#     samples[p_id] = {}
#     samples[p_id]['name'] = param['name']
#     samples[p_id]['trace'] = param['sampler'].trace() 
# 
# 
#   # labels = { 'scale_He':r'$\beta_{\mathrm{He}}$', 'scale_H':r'$\beta_{\mathrm{H}}$', 'deltaZ_He':r'$\Delta z_{\mathrm{He}}$', 'deltaZ_H':r'$\Delta z_{\mathrm{H}}$'    }
#   labels = { 'scale_He':r'$\beta_{\mathrm{He}}$', 'scale_H':r'$\beta_{\mathrm{H}}$', 'deltaZ_He':r'$\Delta z_{\mathrm{He}}$'    }
# 
# 
#   Plot_Corner( samples, labels, output_dir  )
# 
# 
# 
# 
# param_stats = {}
# for p_id in params.keys():
#   p_name = params[p_id]['name']
#   p_stats = stats[p_name]
#   params[p_id]['mean'] = p_stats['mean']
#   params[p_id]['sigma'] = p_stats['standard deviation']
# 
# 
# # Obtain distribution of observables
# n_samples = 10000
# observables = [ 'T0', 'tau' ]
# observables_samples = Sample_Observables( n_samples, observables, params, data_grid, SG  )
# chi2_vals = Get_Chi2( observables, params, comparable_grid, comparable_data, SG )
# Plot_Observables( observables_samples, comparable_data, params, SG, 'sampling', output_dir, chi2=chi2_vals)
# Plot_Observables( observables_samples, comparable_data, params, SG, 'grid', output_dir, chi2=None)
# 
