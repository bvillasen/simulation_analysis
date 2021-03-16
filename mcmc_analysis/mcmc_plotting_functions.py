import os, sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
import pymc 
import palettable
from  load_tabulated_data import *



def Plot_T0_tau_Sampling( samples_fields, comparable_data, output_dir, system='Shamrock' ):
   
  import pylab
  import matplotlib
  import matplotlib.font_manager
  matplotlib.rcParams['mathtext.fontset'] = 'cm'
  matplotlib.rcParams['mathtext.rm'] = 'serif'

  if system == 'Lux':      prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/brvillas/fonts', "Helvetica.ttf"), size=12)
  if system == 'Shamrock': prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=12)

  nrows = 1
  ncols = 2

  font_size = 18
  label_size = 16
  alpha = 0.5

  c_pchw18 = pylab.cm.viridis(.7)
  c_hm12 = pylab.cm.cool(.3)

  c_boss = pylab.cm.viridis(.3)
  c_walther = pylab.cm.viridis(.3)
  c_viel = 'C1'
  c_boera = pylab.cm.Purples(.7)

  text_color  = 'black'
  color_line = c_pchw18
  color_data = c_boss

  fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))
  
  ax = ax_l[0]
  obs_name = 'T0'
  samples = samples_fields[obs_name]
  z = samples['z']
  mean = samples['mean']
  high = samples['higher']
  low = samples['lower']
  if 'Highest_Likelihood' in samples:
    print( 'Plotting Highest_Likelihood T0')
    mean = samples['Highest_Likelihood']
  ax.plot( z, mean, zorder=1 )
  ax.fill_between( z, high, low, alpha=alpha, zorder=1 )  

  data_set = comparable_data[obs_name]
  data_z = data_set['z']
  data_mean = data_set['mean'] 
  data_error = data_set['sigma'] 
  ax.errorbar( data_z, data_mean, yerr=data_error, fmt='none',  alpha=0.8, ecolor= color_data, zorder=2)
  ax.scatter( data_z, data_mean, label='Gaikwad et al.', alpha=0.8, color= color_data, zorder=2) 

  ax.tick_params(axis='both', which='major', direction='in', labelsize=label_size )
  ax.tick_params(axis='both', which='minor', direction='in' )
  ax.set_ylabel( r'$T_0   \,\,\, [\,\mathrm{K}\,]$', fontsize=font_size  )
  ax.set_xlabel( r'$z$', fontsize=font_size )
  leg = ax.legend(loc=1, frameon=False, fontsize=font_size, prop=prop)
  ax.set_xlim( 1.8, 12 )
  ax.set_ylim( 3000, 18000)

  ax = ax_l[1]
  obs_name = 'tau'
  samples = samples_fields[obs_name]
  z = samples['z']
  mean = samples['mean']
  high = samples['higher']
  low = samples['lower']
  if 'Highest_Likelihood' in samples:
    print( 'Plotting Highest_Likelihood T0')
    mean = samples['Highest_Likelihood']
  ax.plot( z, mean, zorder=1 )
  ax.fill_between( z, high, low, alpha=alpha, zorder=1 )  

  data_set = comparable_data[obs_name]
  data_z = data_set['z']
  data_mean = data_set['mean'] 
  data_error = data_set['sigma'] 
  ax.errorbar( data_z, data_mean, yerr=data_error, fmt='none',  alpha=0.8, ecolor= color_data, zorder=2)
  ax.scatter( data_z, data_mean, label='Gaikwad et al.', alpha=0.8, color= color_data, zorder=2) 

  ax.tick_params(axis='both', which='major', direction='in', labelsize=label_size )
  ax.tick_params(axis='both', which='minor', direction='in' )
  ax.set_ylabel( r'$\tau_{eff}$', fontsize=font_size  )
  ax.set_xlabel( r'$z$', fontsize=font_size )
  leg = ax.legend(loc=2, frameon=False, fontsize=font_size, prop=prop)
  ax.set_xlim( 2, 6 )
  ax.set_ylim( 0.1, 8)
  ax.set_yscale('log')
  
  
  figure_name = output_dir + f'fig_T0_tau_sampling.png'
  fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
  print( f'Saved Figure: {figure_name}' )


def Plot_T0_Sampling( samples, comparable_data, output_dir, system='Shamrock' ):
   
  import pylab
  import matplotlib
  import matplotlib.font_manager
  matplotlib.rcParams['mathtext.fontset'] = 'cm'
  matplotlib.rcParams['mathtext.rm'] = 'serif'

  if system == 'Lux':      prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/brvillas/fonts', "Helvetica.ttf"), size=12)
  if system == 'Shamrock': prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=12)

  nrows = 1
  ncols = 1

  font_size = 18
  label_size = 16
  alpha = 0.5

  c_pchw18 = pylab.cm.viridis(.7)
  c_hm12 = pylab.cm.cool(.3)

  c_boss = pylab.cm.viridis(.3)
  c_walther = pylab.cm.viridis(.3)
  c_viel = 'C1'
  c_boera = pylab.cm.Purples(.7)

  text_color  = 'black'
  color_line = c_pchw18
  color_data = c_boss

  fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))

  obs_name = 'T0'
  z = samples['z']
  mean = samples['mean']
  high = samples['higher']
  low = samples['lower']
  if 'Highest_Likelihood' in samples:
    print( 'Plotting Highest_Likelihood T0')
    mean = samples['Highest_Likelihood']
  ax.plot( z, mean, zorder=1 )
  ax.fill_between( z, high, low, alpha=alpha, zorder=1 )  

  data_set = comparable_data[obs_name]
  data_z = data_set['z']
  data_mean = data_set['mean'] 
  data_error = data_set['sigma'] 
  ax.errorbar( data_z, data_mean, yerr=data_error, fmt='none',  alpha=0.8, ecolor= color_data, zorder=2)
  ax.scatter( data_z, data_mean, label='Gaikwad et al.', alpha=0.8, color= color_data, zorder=2) 

  ax.tick_params(axis='both', which='major', direction='in', labelsize=label_size )
  ax.tick_params(axis='both', which='minor', direction='in' )
  ax.set_ylabel( r'$T_0   \,\,\, [\,\mathrm{K}\,]$', fontsize=font_size  )
  ax.set_xlabel( r'$z$', fontsize=font_size )
  leg = ax.legend(loc=1, frameon=False, fontsize=font_size, prop=prop)
  ax.set_xlim( 1.8, 12 )
  ax.set_ylim( 3000, 18000)

  figure_name = output_dir + f'fig_T0_sampling.png'
  fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
  print( f'Saved Figure: {figure_name}' )





def Plot_Comparable_Data( field, comparable_data, comparable_grid, output_dir  ):

  nrows, ncols = 1, 1
  fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20*ncols,5*nrows))

  data_mean = comparable_data[field]['mean']
  data_sigma = comparable_data[field]['sigma']
  n_points = len( data_mean )
  x = np.arange( 0, n_points, 1)

  ax.errorbar( x, data_mean, yerr=data_sigma, fmt='o', c='C0', label='Data', ms=1)

  sim_ids = comparable_grid.keys()
  for sim_id in sim_ids:
    sim_mean = comparable_grid[sim_id][field]['mean']
    ax.scatter(x, sim_mean, s=1 )


  ax.set_yscale('log')
  ax.legend( frameon=False )

  figure_name = output_dir + 'data_for_fit.png'
  fig.savefig( figure_name, bbox_inches='tight', dpi=500 )
  print( f'Saved Figure: {figure_name}' )


def Plot_Corner( samples, labels, output_dir  ):
  param_ids = samples.keys()
  n_param = len( param_ids )
  color = 'C0'
  data_color = 'C9'
  font_size = 16
  label_size = 20
  alpha = 0.6
  fig_size = 5
  space = 0.05
  tick_label_size = 13
  tick_length = 7
  tick_width = 2
  border_width = 2.0
  n_bins_1D = 20
  hist_1D_line_width = 2
  hist_1D_line_color = 'C0'

  n_bins_2D = 30
  hist_2D_colormap = palettable.cmocean.sequential.Ice_20_r.mpl_colormap

  fig, ax_l = plt.subplots(nrows=n_param, ncols=n_param, figsize=(fig_size*n_param,fig_size*n_param),  sharex=False,)
  fig.subplots_adjust( wspace=space, hspace=space )

  for j in range( n_param ):
    for i in range( n_param ):

      ax = ax_l[j][i]
      plot_y_lables, plot_x_lables = False, False
      plot_y_ticks  = True
      if i == 0: plot_y_lables = True
      if j == n_param-1: plot_x_lables = True 
      if i == j: 
        plot_y_lables = False
        plot_y_ticks = False

      if i > j:
        ax.axis("off")
        continue

      if plot_x_lables: ax.tick_params(axis='x', which='major', direction='in', labelsize=tick_label_size, length=tick_length, width=tick_width )
      else:             ax.tick_params(axis='x', which='major', direction='in', labelsize=0, length=tick_length, width=tick_width )
      if plot_y_lables: ax.tick_params(axis='y', which='major', direction='in', labelsize=tick_label_size, length=tick_length, width=tick_width )
      else:             ax.tick_params(axis='y', which='major', direction='in', labelsize=0, length=tick_length, width=tick_width )
      if not plot_y_ticks: ax.tick_params(axis='y', which='major', length=0 )

      if plot_y_lables:
        if j == 0: y_label = ''
        else: y_label = labels[samples[j]['name']]  
        ax.set_ylabel( y_label, fontsize=label_size )

      if plot_x_lables:
        x_label = labels[samples[i]['name']]  
        ax.set_xlabel( x_label, fontsize=label_size )

      if i == j: plot_type = '1D'
      if i < j:  plot_type = '2D'


      if plot_type == '1D':
        name  = samples[j]['name']
        trace = samples[j]['trace']
        hist, bin_edges = np.histogram( trace, bins=n_bins_1D ) 
        bin_centers = ( bin_edges[:-1] + bin_edges[1:] ) / 2.
        bin_width = bin_centers[0] - bin_centers[1]  
        ax.step( bin_centers, hist, where='mid', color=hist_1D_line_color, linewidth=hist_1D_line_width  )

      if plot_type == '2D':
        trace_y = samples[j]['trace']
        trace_x = samples[i]['trace']
        hist, x_edges, y_edges = np.histogram2d( trace_x, trace_y, bins=[n_bins_2D, n_bins_2D] )
        hist = hist.astype( np.float ) 
        hist = hist.T / hist.sum()
        hist_mean = hist.mean()
        hist_sigma = np.sqrt( ((hist-hist_mean)**2).mean() )
        extent = [ trace_x.min(), trace_x.max(), trace_y.min(), trace_y.max() ]
        ax.imshow( hist[::-1], cmap=hist_2D_colormap, extent=extent, aspect='auto' )
        ax.contour( hist, [hist_sigma, 2* hist_sigma], extent=extent )

      [sp.set_linewidth(border_width) for sp in ax.spines.values()]

  figure_name = output_dir + 'corner.png'
  fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
  print( f'Saved Figure: {figure_name}' )




def Plot_MCMC_Stats( stats, MDL, params_mcmc,  stats_file, output_dir, plot_corner=True ):
  cwd = os.getcwd()
  os.chdir( output_dir )
  
  pymc.Matplot.plot(MDL)  

  labels, samples = [], []
  for p_id in params_mcmc.keys():
    param = params_mcmc[p_id]
    labels.append( param['name'] )
    samples.append( param['sampler'].trace())
  samples = np.array( samples ).T

  if plot_corner:
    import corner
    corner_fig = corner.corner(samples[:,:], labels=labels )
    corner_fig.savefig( 'corner_fig.png' )  
  os.chdir( cwd )  


def Plot_Observables( observables_samples, comparable_data, params, SG, plot_type, output_dir, chi2=None):


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
  
  sim_ids = SG.sim_ids
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
      # param_val = SG.Grid[sim_id]['parameters'][param_name]
      # if param_name == 'scale_He': label_param = r'$\beta_{HeII}$' 
      # if param_name == 'scale_H': label_param = r'$\beta_{HI}$' 
      # if param_name == 'deltaZ_He': label_param = r'$\Delta z_{HeII}$' 
      # if param_name == 'deltaZ_H': label_param = r'$\Delta z_{HI}$' 
      # label =  label_param + ' $= {0}$'.format(param_val)
      label = ''
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

  if chi2 != None and plot_type=='sampling':
    chi2_val = chi2['T0']
    text = r'$\chi ^2 = {0:.2f}$ '.format( chi2_val) 
    ax.text(0.7, 0.65, text, horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=font_size, )

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
      # param_val = SG.Grid[sim_id]['parameters'][param_name]
      # if param_name == 'scale_H': label_param = r'$\beta_{HI}$' 
      # if param_name == 'scale_He': label_param = r'$\beta_{HeII}$' 
      # label =  label_param + ' $= {0}$'.format(param_val)
      lavel = ''
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

  if chi2 != None and plot_type=='sampling':
    chi2_val = chi2['tau']
    text = r'$\chi ^2 = {0:.2f}$ '.format( chi2_val) 
    ax.text(0.2, 0.65, text, horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=font_size, )
  figure_name = output_dir + f'fig_composite_{plot_type}.png'
  fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
  print( f'Saved Figure: {figure_name}' )




def Plot_Power_Spectrum_Sampling( ps_samples, ps_data_dir, output_dir, scales='small', system=None, name='' ):

  import matplotlib
  import pylab

  matplotlib.rcParams['mathtext.fontset'] = 'cm'
  matplotlib.rcParams['mathtext.rm'] = 'serif'


  if system == 'Lux': matplotlib.use('Agg')

  print( f'Plotting P(k) Scales: {scales}' )

  fig_height = 5
  fig_width = 8
  fig_dpi = 300


  label_size = 18
  figure_text_size = 18
  legend_font_size = 16
  tick_label_size_major = 15
  tick_label_size_minor = 13
  tick_size_major = 5
  tick_size_minor = 3
  tick_width_major = 1.5
  tick_width_minor = 1
  border_width = 1

  if system == 'Lux':      prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/brvillas/fonts', "Helvetica.ttf"), size=12)
  if system == 'Shamrock': prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=12)

  dir_boss = ps_data_dir + 'data_power_spectrum_boss/'
  data_filename = dir_boss + 'data_table.py'
  data_boss = load_data_boss( data_filename )
  data_z_boss = data_boss['z_vals']

  data_filename = ps_data_dir + 'data_power_spectrum_walther_2019/data_table.txt'
  data_walther = load_power_spectrum_table( data_filename )
  data_z_w = data_walther['z_vals']

  dir_data_boera = ps_data_dir + 'data_power_spectrum_boera_2019/'
  data_boera = load_tabulated_data_boera( dir_data_boera )
  data_z_b = data_boera['z_vals']

  data_dir_viel = ps_data_dir + 'data_power_spectrum_viel_2013/'
  data_viel = load_tabulated_data_viel( data_dir_viel)
  data_z_v = data_viel['z_vals']

  z_vals_small_scale  = [ 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 4.2, 4.6, 5.0, 5.4 ]
  z_vals_large_scale  = [ 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.6 ]
  z_vals_middle_scale = [ 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 4.2, 4.6 ]
  z_vals_all_scale    = z_vals_large_scale
  z_high = [ 5.0, 5.4 ]



  n_samples = len( ps_samples )
  z_samples = np.array([ ps_samples[i]['z'] for i in range(n_samples) ])

  if scales == 'large':    z_vals = z_vals_large_scale
  elif scales == 'small':  z_vals = z_vals_small_scale
  elif scales == 'middle': z_vals = z_vals_middle_scale
  elif scales == 'all':    z_vals = z_vals_all_scale
  else: 
    print( "ERROR: Scales = large,  small, middle or all ")
    # return


  nrows = 3
  ncols = 4


  if scales == 'middle':flags = np.zeros( (nrows, ncols ))

  # if scales == 'middle': nrows = 2


  fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(2*fig_width,fig_height*nrows))
  plt.subplots_adjust( hspace = 0.02, wspace=0.02)


  c_pchw18 = pylab.cm.viridis(.7)
  c_hm12 = pylab.cm.cool(.3)

  c_boss = pylab.cm.viridis(.3)
  c_walther = pylab.cm.viridis(.3)
  c_viel = 'C1'
  c_boera = pylab.cm.Purples(.7)

  text_color  = 'black'
  color_line = c_pchw18

  if scales == 'middle' or scales == 'all':
    c_walther = c_pchw18

  alpha_bar = 0.4

  for index, current_z in enumerate( z_vals ):



    indx_j = index % ncols
    indx_i = index//ncols

    if nrows > 1: ax = ax_l[indx_i][indx_j]
    else: ax = ax_l[indx_j]


    if scales == 'middle': flags[indx_i,  indx_j] = 1

    if ps_samples != []:
      diff = np.abs( z_samples - current_z )
      diff_min = diff.min()
      if diff_min < 0.05:
        index = np.where( diff == diff_min )[0][0]
        k_vals = ps_samples[index]['k_vals']
        delta_mean = ps_samples[index]['mean']
        delta_sigma = ps_samples[index]['sigma'] 
        delta_higher = ps_samples[index]['higher']
        delta_lower = ps_samples[index]['lower']
        delta_p = delta_mean + delta_sigma
        delta_m = delta_mean - delta_sigma 
        
        if 'Highest_Likelihood' in ps_samples[index]:
          print( 'Plotting Highest_Likelihood PS')
          delta_mean = ps_samples[index]['Highest_Likelihood']
        else:
          print( 'Plotting mean of distribution (Not Highest_Likelihood)')
        ax.plot( k_vals, delta_mean, linewidth=3, color=color_line, zorder=1   )
        # ax.fill_between( k_vals, delta_p, delta_m, facecolor=color_line, alpha=alpha_bar, zorder=1   )
        ax.fill_between( k_vals, delta_higher, delta_lower, facecolor=color_line, alpha=alpha_bar, zorder=1   )

    #       # ax.plot( k, delta, c=color_line, linewidth=3, label=sim_data['plot_label']  )
    # 

    ax.text(0.85, 0.95, r'$z={0:.1f}$'.format(current_z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size, color=text_color) 


    if scales == 'large' or scales == 'middle' or scales == 'all':

      # Add Boss data
      z_diff = np.abs( data_z_boss - current_z )
      diff_min = z_diff.min()
      if diff_min < 1e-1:
        data_index = np.where( z_diff == diff_min )[0][0]
        data_z_local = data_z_boss[data_index]

        data_k = data_boss[data_index]['k_vals']
        data_delta_power = data_boss[data_index]['delta_power']
        data_delta_power_error = data_boss[data_index]['delta_power_error']
        label_boss = 'eBOSS (2019)'
        d_boss = ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_boss, label=label_boss, zorder=2)

    if scales == 'small' or scales == 'middle' or scales == 'all':

      # Add Walther data
      z_diff = np.abs( data_z_w - current_z )
      diff_min = z_diff.min()
      if diff_min < 1e-1:
        data_index = np.where( z_diff == diff_min )[0][0]
        data_z_local = data_z_w[data_index]

        data_k = data_walther[data_index]['k_vals']
        data_delta_power = data_walther[data_index]['delta_power']
        data_delta_power_error = data_walther[data_index]['delta_power_error']
        label_walther ='Walther et al. (2018)' 
        d_walther = ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_walther, label=label_walther, zorder=2)


      # Add Boera data
      z_diff = np.abs( data_z_b - current_z )
      diff_min = z_diff.min()
      if diff_min < 1e-1:
        data_index = np.where( z_diff == diff_min )[0][0]
        data_z_local = data_z_b[data_index]

        data_k = data_boera[data_index]['k_vals']
        data_delta_power = data_boera[data_index]['delta_power']
        data_delta_power_error = data_boera[data_index]['delta_power_error']
        label_boera ='Boera et al. (2019)'
        d_boera = ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_boera, label=label_boera, zorder=2 )


      # Add Viel data
      z_diff = np.abs( data_z_v - current_z )
      diff_min = z_diff.min()
      if diff_min < 1e-1:
        data_index = np.where( z_diff == diff_min )[0][0]
        data_z_local = data_z_v[data_index]

        data_k = data_viel[data_index]['k_vals']
        data_delta_power = data_viel[data_index]['delta_power']
        data_delta_power_error = data_viel[data_index]['delta_power_error']
        label_viel = 'Viel et al. (2013)'
        d_viel = ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_viel, label=label_viel, zorder=2 )


    legend_loc = 3
    if indx_i == nrows-1 and nrows!=2: legend_loc = 2

    if scales == 'large': legend_loc = 2
    label_bars =  r'1$\sigma$ skewers $P\,(\Delta_F^2)$'

    add_legend = False
    if indx_j == 0: add_legend = True

    if scales == 'middle' and indx_i == nrows-1 and indx_j == ncols-1: add_legend = True


    if add_legend:
      # leg = ax.legend( loc=legend_loc, frameon=False, fontsize=12)
      leg = ax.legend(  loc=legend_loc, frameon=False, prop=prop    )

      for text in leg.get_texts():
          plt.setp(text, color = text_color)



    if scales == 'small':
      x_min, x_max = 4e-3, 2.5e-1
      if indx_i == 0: y_min, y_max = 1e-3, 9e-2
      if indx_i == 1: y_min, y_max = 5e-3, 2e-1
      if indx_i == 2: y_min, y_max = 5e-2, 3

    if scales == 'large':
      x_min, x_max = 2e-3, 2.3e-2
      if indx_i == 0: y_min, y_max = 1e-2, 1.2e-1
      if indx_i == 1: y_min, y_max = 2e-2, 2.5e-1
      if indx_i == 2: y_min, y_max = 5e-2, 7e-1

    if scales == 'middle':
      x_min, x_max = 5e-3, 1e-1
      if indx_i == 0: y_min, y_max = 4e-3, 9e-2
      if indx_i == 1: y_min, y_max = 1e-2, 5e-1
      
    if scales == 'all':
      x_min, x_max = 1e-3, 1e0
      if indx_i == 0: y_min, y_max = 8e-5, 9e-2
      if indx_i == 1: y_min, y_max = 8e-5, 2e-1
      if indx_i == 2: y_min, y_max = 5e-4, 6e-1


    ax.set_xlim( x_min, x_max )
    ax.set_ylim( y_min, y_max )
    ax.set_xscale('log')
    ax.set_yscale('log')


    [sp.set_linewidth(border_width) for sp in ax.spines.values()]

    if indx_j > 0:ax.set_yticklabels([])
    if indx_i != nrows-1 :ax.set_xticklabels([])

    ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in' )
    ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')

    if indx_j == 0: ax.set_ylabel( r' $\Delta_F^2(k)$', fontsize=label_size, color= text_color )
    if indx_i == nrows-1: ax.set_xlabel( r'$ k   \,\,\,  [\mathrm{s}\,\mathrm{km}^{-1}] $',  fontsize=label_size, color= text_color )



  if scales == 'middle':
    for i in range( nrows ):
      for j in range( ncols ):
        if not flags[i,j]:
          ax = ax_l[i][j].axis('off')


  fileName = output_dir + f'flux_ps_sampling_{scales}'
  if name != '': fileName += f'_{name}'
  fileName += '.png'
  # fileName += '.pdf'
  fig.savefig( fileName,  pad_inches=0.1, bbox_inches='tight', dpi=fig_dpi)
  print('Saved Image: ', fileName)