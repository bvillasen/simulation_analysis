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
# data_ps_sets = [ 'Walther' ]
# data_ps_sets = [ 'Boera' ]
# data_ps_sets = [ 'Boss', 'Walther' ]
# data_ps_sets = [ 'Walther', 'Boera' ]
# data_ps_sets = [ 'Walther', 'Viel' ]
# data_ps_sets = [ 'Boss', 'Walther', 'Boera' ]
# data_ps_sets = [ 'Boss', 'Walther', 'Viel' ]


name = ''
for data_set in data_ps_sets:
  name += data_set + '_'
name = name[:-1] 

# field = 'P(k)+T0'
# field = 'P(k)+'
field = 'P(k)+tau_HeII'

# fit_log_power_spectrum =  True
fit_log_power_spectrum =  False
if fit_log_power_spectrum: name += '_log'

ps_data_dir = 'lya_statistics/data/'
mcmc_dir = root_dir + 'fit_mcmc/'
create_directory( mcmc_dir )
output_dir = mcmc_dir + f'fit_results_{field}_{name}/'
create_directory( output_dir )

# load_mcmc_results = False
load_mcmc_results = True


SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Load_Grid_Analysis_Data()
ps_range = SG.Get_Power_Spectrum_Range( kmax=1e-1 )
sim_ids = SG.sim_ids

z_min = 2.0
z_max = 5.0 
ps_extras = { 'range':ps_range, 'data_dir':ps_data_dir, 'data_sets':data_ps_sets }
comparable_data = Get_Comparable_Composite( field,  z_min, z_max, ps_extras=ps_extras, log_ps=fit_log_power_spectrum )
comparable_grid = Get_Comparable_Composite_from_Grid( field, comparable_data, SG, log_ps=fit_log_power_spectrum )
Plot_Comparable_Data( field, comparable_data, comparable_grid, output_dir, log_ps=fit_log_power_spectrum  )

z_vals = [ 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.6, 5.0,  ]
data_grid, data_grid_power_spectrum = Get_Data_Grid_Composite( ['P(k)', 'T0', 'tau', 'tau_HeII'], SG, z_vals=z_vals )


stats_file = output_dir + 'fit_mcmc.pkl'
samples_file = output_dir + 'samples_mcmc.pkl'


params = SG.parameters

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




# Get the Highest_Likelihood parameter values 
params_HL = Get_Highest_Likelihood_Params( param_samples, n_bins=100 )
scale_He_HL, scale_H_HL, deltaZ_He_HL, deltaZ_H_HL = params_HL.flatten()


scale_He_values  = params[0]['values']
scale_H_values   = params[1]['values']
deltaZ_He_values = params[2]['values']
deltaZ_H_values  = params[3]['values']


scale_He = scale_He_values[0]




samples_ps = {}
for indx, scale_He in enumerate( scale_He_values  ):
  p_vals = [ scale_He, scale_H_HL, deltaZ_He_HL, deltaZ_H_HL ]
  ps_interp = Interpolate_Power_Spectrum( p_vals, data_grid_power_spectrum, SG )  
  ps_interp['label'] = r'$\beta_{\mathrm{He}}:$' + f'{scale_He:.2f}' 
  samples_ps[indx] = ps_interp
Plot_Power_Spectrum_Sampling( samples_ps, ps_data_dir, output_dir, scales='large', system=system, label=None, plot_type='multiple_lines', name='scale_He' )
Plot_Power_Spectrum_Sampling( samples_ps, ps_data_dir, output_dir, scales='small', system=system, label=None, plot_type='multiple_lines', name='scale_He' )
Plot_Power_Spectrum_Sampling( samples_ps, ps_data_dir, output_dir, scales='middle', system=system, label=None, plot_type='multiple_lines', name='scale_He' )
  


samples_ps = {}
for indx, scale_H in enumerate( scale_H_values  ):
  p_vals = [ scale_He_HL, scale_H, deltaZ_He_HL, deltaZ_H_HL ]
  ps_interp = Interpolate_Power_Spectrum( p_vals, data_grid_power_spectrum, SG )  
  ps_interp['label'] = r'$\beta_{\mathrm{H}}:$' + f'{scale_H:.2f}' 
  samples_ps[indx] = ps_interp
Plot_Power_Spectrum_Sampling( samples_ps, ps_data_dir, output_dir, scales='large', system=system, label=None, plot_type='multiple_lines', name='scale_H' )
Plot_Power_Spectrum_Sampling( samples_ps, ps_data_dir, output_dir, scales='small', system=system, label=None, plot_type='multiple_lines', name='scale_H' )
Plot_Power_Spectrum_Sampling( samples_ps, ps_data_dir, output_dir, scales='middle', system=system, label=None, plot_type='multiple_lines', name='scale_H' )


samples_ps = {}
for indx, deltaZ_He in enumerate( deltaZ_He_values  ):
  p_vals = [ scale_He_HL, scale_H_HL, deltaZ_He, deltaZ_H_HL ]
  ps_interp = Interpolate_Power_Spectrum( p_vals, data_grid_power_spectrum, SG )  
  ps_interp['label'] = r'$\Delta z_{\mathrm{He}}:$' + f'{deltaZ_He:.2f}' 
  samples_ps[indx] = ps_interp
Plot_Power_Spectrum_Sampling( samples_ps, ps_data_dir, output_dir, scales='large', system=system, label=None, plot_type='multiple_lines', name='deltaZ_He' )
Plot_Power_Spectrum_Sampling( samples_ps, ps_data_dir, output_dir, scales='small', system=system, label=None, plot_type='multiple_lines', name='deltaZ_He' )
Plot_Power_Spectrum_Sampling( samples_ps, ps_data_dir, output_dir, scales='middle', system=system, label=None, plot_type='multiple_lines', name='deltaZ_He' )


samples_ps = {}
for indx, deltaZ_H in enumerate( deltaZ_H_values  ):
  p_vals = [ scale_He_HL, scale_H_HL, deltaZ_He_HL, deltaZ_H ]
  ps_interp = Interpolate_Power_Spectrum( p_vals, data_grid_power_spectrum, SG )  
  ps_interp['label'] = r'$\Delta z_{\mathrm{H}}:$' + f'{deltaZ_H:.2f}' 
  samples_ps[indx] = ps_interp
Plot_Power_Spectrum_Sampling( samples_ps, ps_data_dir, output_dir, scales='large', system=system, label=None, plot_type='multiple_lines', name='deltaZ_H' )
Plot_Power_Spectrum_Sampling( samples_ps, ps_data_dir, output_dir, scales='small', system=system, label=None, plot_type='multiple_lines', name='deltaZ_H' )
Plot_Power_Spectrum_Sampling( samples_ps, ps_data_dir, output_dir, scales='middle', system=system, label=None, plot_type='multiple_lines', name='deltaZ_H' )

