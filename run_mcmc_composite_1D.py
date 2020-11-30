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



SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Load_Grid_Analysis_Data()
sim_ids = SG.sim_ids

comparable_data = Get_Comparable_Composite_T0_tau()
comparable_grid = Get_Comparable_Composite_T0_tau_from_Grid( comparable_data, SG )

field = 'T0+tau'
param_to_fit = 1
param_name = SG.parameters[param_to_fit]['name']

nIter = 100000 
nBurn = nIter / 5
nThin = 1
model = mcmc_model_1D( param_to_fit, comparable_data, comparable_grid, field, SG )
MDL = pymc.MCMC( model )  
MDL.sample( iter=nIter, burn=nBurn, thin=nThin )
stats = MDL.stats()

param_stats = stats[param_name]
param_mean  = param_stats['mean']
param_sigma = param_stats['standard deviation']
print( f'\n {param_name}: {param_mean}  +- {param_sigma}' )



n_samples = 10000
observables_to_sample = ['T0', 'tau']

observables = { observable:{} for observable in observables_to_sample }
for observable in observables_to_sample:
  observables[observable]['samples'] = []
  for i in range ( n_samples ):
    param_val = np.random.normal( param_mean, param_sigma )
    observable_interp, z_interp =  Interpolate_Observable_1D( param_to_fit, observable, param_val,  SG )
    observables[observable]['samples'].append( observable_interp )
    observables[observable]['z'] = z_interp
  obs_all = np.array(observables[observable]['samples']).T
  obs_mean = [ obs_vals.mean() for obs_vals in obs_all ]
  obs_sigma = []
  for i in range( len(obs_all) ):
    obs_sigma.append( np.sqrt( (( obs_all[i] - obs_mean[i] )**2 ).mean() ) )
  observables[observable]['mean'] = np.array( obs_mean )
  observables[observable]['sigma'] = np.array( obs_sigma ) 




output_dir = root_dir + 'fit_results_composite/'
create_directory( output_dir )



nrows = 1
ncols = 2

color = 'C0'
data_color = 'C9'
font_size = 15
label_size = 14
alpha = 0.8


for plot_type in ['grid', 'sampling']:

  fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))

  ax = ax_l[0]
  obs_name = 'T0'
  obs_z     = observables[obs_name]['z']
  obs_mean  = observables[obs_name]['mean']
  obs_sigma = observables[obs_name]['sigma'] 

  if plot_type == 'sampling':
    label = '{0} = {1:.2}'.format(param_name, param_mean) + r' $\pm$ '+ '{0:.2}'.format( param_sigma )
    ax.fill_between( obs_z, obs_mean + obs_sigma, obs_mean - obs_sigma, alpha=alpha, color=color, zorder=1)
    ax.plot( obs_z, obs_mean, c=color, label=label, zorder=1 )

  if plot_type == 'grid':
    for sim_id in sim_ids:
      data_sim = SG.Grid[sim_id]['analysis']
      z = data_sim['z']
      obs_vals = data_sim[obs_name]
      param_val = SG.Grid[sim_id]['parameters'][param_name]
      if param_name == 'scale_He': label_param = r'$\beta_{HeII}$' 
      if param_name == 'scale_H': label_param = r'$\beta_{HI}$' 
      if param_name == 'deltaZ_He': label_param = r'$\Delta z_{HeII}$' 
      if param_name == 'deltaZ_H': label_param = r'$\Delta z_{HI}$' 
      if param_name == 'scale_H_photoion': label_param = r'$\beta_{HI}$ ionization'
      if param_name == 'scale_H_photoheat': label_param = r'$\beta_{HI}$ heating' 
      
      label =  label_param + ' $= {0}$'.format(param_val)
      ax.plot( z, obs_vals , label=label, zorder=1 )

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
  obs_z     = observables[obs_name]['z']
  obs_mean  = observables[obs_name]['mean']
  obs_sigma = observables[obs_name]['sigma'] 

  if plot_type == 'sampling':
    label = '{0} = {1:.2}'.format(param_name, param_mean) + r' $\pm$ '+ '{0:.2}'.format( param_sigma )
    ax.fill_between( obs_z, obs_mean + obs_sigma, obs_mean - obs_sigma, alpha=alpha, color=color, zorder=1)
    ax.plot( obs_z, obs_mean, c=color, label=label, zorder=1 )

  if plot_type == 'grid':
    for sim_id in sim_ids:
      data_sim = SG.Grid[sim_id]['analysis']
      z = data_sim['z']
      obs_vals = data_sim[obs_name]
      param_val = SG.Grid[sim_id]['parameters'][param_name]
      if param_name == 'scale_H': label_param = r'$\beta_{HI}$' 
      if param_name == 'scale_He': label_param = r'$\beta_{HeII}$' 
      if param_name == 'scale_H_photoion': label_param = r'$\beta_{HI}$ ionization'
      if param_name == 'scale_H_photoheat': label_param = r'$\beta_{HI}$ heating' 
      label =  label_param + ' $= {0}$'.format(param_val)
      ax.plot( z, obs_vals , label=label, zorder=1 )

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
  ax.set_xlim( 2, 6 )
  ax.set_ylim( 0.1, 8)
  ax.set_yscale('log')

  figure_name = output_dir + f'fig_composite_{param_name}_{plot_type}.png'
  fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
  print( f'Saved Figure: {figure_name}' )





cwd = os.getcwd()
os.chdir( output_dir )

out_file_name = output_dir + 'fit_mcmc.pkl'
f = open( out_file_name, "wb")
pickle.dump( stats, f)
print ( f'Saved File: {out_file_name}' )

pymc.Matplot.plot(MDL)  

os.chdir( cwd )  
