import os, sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
import pymc 

def Plot_MCMC_Stats( stats, MDL, params_mcmc,  stats_file, output_dir ):
  cwd = os.getcwd()
  os.chdir( output_dir )

  f = open( stats_file, 'wb' )
  pickle.dump( stats, f)
  print ( f'Saved File: {stats_file}' )
  pymc.Matplot.plot(MDL)  

  labels, samples = [], []
  for p_id in params_mcmc.keys():
    param = params_mcmc[p_id]
    labels.append( param['name'] )
    samples.append( param['sampler'].trace())
  samples = np.array( samples ).T

  import corner
  corner_fig = corner.corner(samples[:,:], labels=labels )
  corner_fig.savefig( 'corner_fig.png' )  
  os.chdir( cwd )  


def Plot_Observables( observables_samples, comparable_data, params, SG, plot_type, output_dir):

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
  obs_z     = observables_samples[obs_name]['z']
  obs_mean  = observables_samples[obs_name]['mean']
  obs_sigma = observables_samples[obs_name]['sigma'] 

  if plot_type == 'sampling':
    label = ''
    for p_id in params.keys():
      p_name = params[p_id]['name']
      p_mean = params[p_id]['mean']
      p_sigma = params[p_id]['sigma']
      label += '{0} = {1:.2}'.format(p_name, p_mean) + r' $\pm$ '+ '{0:.2}\n'.format( p_sigma )
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
  obs_z     = observables_samples[obs_name]['z']
  obs_mean  = observables_samples[obs_name]['mean']
  obs_sigma = observables_samples[obs_name]['sigma'] 

  if plot_type == 'sampling':
    # label = '{0} = {1:.2}'.format(param_name, param_mean) + r' $\pm$ '+ '{0:.2}'.format( param_sigma )
    # label = ''
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

  figure_name = output_dir + f'fig_composite_{plot_type}.png'
  fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
  print( f'Saved Figure: {figure_name}' )
