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

field = 'T0+tau'
output_dir = root_dir + f'fit_results_{field}/'
create_directory( output_dir )



SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Load_Grid_Analysis_Data()
sim_ids = SG.sim_ids

comparable_data = Get_Comparable_Composite_T0_tau( factor_sigma_tau_becker=6.0, factor_sigma_tau_keating=4.0,  )
comparable_grid = Get_Comparable_Composite_T0_tau_from_Grid( comparable_data, SG )

fields = [ 'T0', 'tau' ]
data_grid = Get_Data_Grid( fields, SG )

params = SG.parameters




nrows = 1
ncols = 2
color = 'C0'
data_color = 'C9'
font_size = 15
label_size = 14
alpha = 0.6
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))
ax = ax_l[0]


obs_name = 'T0'

params = {}
params[0] = {'mean':0.41}
params[1] = {'mean':0.82}
params[2] = {'mean':0.13}
params[3] = {}

# deltaZ_H_vals  = [ -0.2, -0.1, -0.05, 0.05, 0.1, 0.15, 0.2 ]
deltaZ_H_vals  = [ -0.15, 0.10, 0.15, 0.2 ]
z = data_grid[0]['z']
param_name = r'$\Delta z_{\mathrm{H}}=$'
for deltaZ_H in deltaZ_H_vals:
  params[3] = {'mean':deltaZ_H}
  obs_interp = Interpolate_MultiDim(  params[0]['mean'], params[1]['mean'], params[2]['mean'], params[3]['mean'], data_grid, obs_name, 'mean', SG, clip_params=True ) 
  chi2_vals = Get_Chi2( ['T0'], params, comparable_grid, comparable_data, SG )
  print(chi2_vals)
  label = param_name + '{0:.2f}   '.format(deltaZ_H) + r'$\chi ^2 =$' + '{0:.2f}   '.format(chi2_vals['T0'])
  ax.plot( z, obs_interp , label=label, zorder=1 )


data_set = comparable_data[obs_name]
data_z = data_set['z']
data_mean = data_set['mean'] 
data_error = data_set['sigma'] 
ax.errorbar( data_z, data_mean, yerr=data_error, fmt='none',  alpha=0.8, ecolor= data_color, zorder=2)
ax.scatter( data_z, data_mean, label='Data for MCMC fit', alpha=0.8, color= data_color, zorder=2) 

ax.tick_params(axis='both', which='major', direction='in', labelsize=label_size )
ax.tick_params(axis='both', which='minor', direction='in' )
ax.set_ylabel( r'$T_0   \,\,\, [\,\mathrm{K}\,]$', fontsize=font_size  )
ax.set_xlabel( r'$z$', fontsize=font_size )
leg = ax.legend(loc=1, frameon=False, fontsize=font_size)
ax.set_xlim( 2, 12 )
ax.set_ylim( 3000, 18000)

ax = ax_l[1]
obs_name = 'tau'

z = data_grid[0]['z']
param_name = r'$\Delta z_{\mathrm{H}}=$'
for deltaZ_H in deltaZ_H_vals:
  params[3] = {'mean':deltaZ_H}
  obs_interp = Interpolate_MultiDim(  params[0]['mean'], params[1]['mean'], params[2]['mean'], params[3]['mean'], data_grid, obs_name, 'mean', SG, clip_params=True ) 
  chi2_vals = Get_Chi2( ['tau'], params, comparable_grid, comparable_data, SG )
  print(chi2_vals)
  label = param_name + '{0:.2f}   '.format(deltaZ_H) + r'$\chi ^2 =$' + '{0:.2f}   '.format(chi2_vals['tau'])
  ax.plot( z, obs_interp , label=label, zorder=1 )


data_set = comparable_data[obs_name]
data_z = data_set['z']
data_mean = data_set['mean'] 
data_error = data_set['sigma'] 
ax.errorbar( data_z, data_mean, yerr=data_error, fmt='none',  alpha=0.8, ecolor= data_color, zorder=2)
ax.scatter( data_z, data_mean, label='Data for MCMC fit', alpha=0.8, color= data_color, zorder=2) 

ax.tick_params(axis='both', which='major', direction='in', labelsize=label_size )
ax.tick_params(axis='both', which='minor', direction='in' )
ax.set_ylabel( r'$\tau_{eff}$', fontsize=font_size  )
ax.set_xlabel( r'$z$', fontsize=font_size )
leg = ax.legend(loc=2, frameon=False, fontsize=font_size)
ax.set_xlim( 2, 6.3 )
ax.set_ylim( 0.1, 8)
ax.set_yscale('log')

figure_name = output_dir + f'fig_composite_interpolation.png'
fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
print( f'Saved Figure: {figure_name}' )

