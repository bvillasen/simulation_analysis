import os, sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
import pymc 
import palettable
import pylab
from load_tabulated_data import *
from data_optical_depth import *
from data_optical_depth_HeII import data_tau_HeII_Worserc_2019
from data_thermal_history import *
analysis_dir = os.path.dirname(os.getcwd()) + '/'
sys.path.append(analysis_dir + 'phase_diagram')
sys.path.append(analysis_dir + 'lya_statistics')
sys.path.append(analysis_dir + 'tools')
from tools import *

color_map_0 = palettable.cmocean.sequential.Ice_20_r
color_map_1 = palettable.cmocean.sequential.Amp_20
color_map_2 = palettable.cmocean.sequential.Tempo_20
color_map_3 = palettable.cmocean.sequential.Dense_20
color_map_4 = palettable.cmocean.sequential.Algae_20
color_map_list = [ color_map_0, color_map_1, color_map_2, color_map_3, color_map_4 ]

c_boss = pylab.cm.viridis(.3)
c_walther = pylab.cm.viridis(.3)

c_0 = pylab.cm.viridis(.7)
c_1 = pylab.cm.cool(.3)
c_2 = 'C3'
c_3 = pylab.cm.Purples(.7)
c_4 = 'C1'

color_lines_list = [ c_0, c_1, c_2, c_3  ] 

use_color_from_colormap = False

def Plot_tau_HeII_Sampling( samples_tau_H, samples_tau_HeII, output_dir, system='Shamrock', label='', multiple=False ):
  
  if not multiple:
    labels_multiple = [label]
    samples_tau_H_multiple = {}
    samples_tau_H_multiple[0] = samples_tau_H
    samples_tau_HeII_multiple = {}
    samples_tau_HeII_multiple[0] = samples_tau_HeII
  else:
    labels_multiple = label
    samples_tau_H_multiple = samples_tau_H
    samples_tau_HeII_multiple = samples_tau_HeII
    
   
  from scipy import interpolate as interp 
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
  alpha = 0.2
  linewidth = 2

  c_pchw18 = pylab.cm.viridis(.7)
  c_hm12 = pylab.cm.cool(.3)

  c_boss = pylab.cm.viridis(.3)
  c_walther = pylab.cm.viridis(.3)
  c_viel = 'C1'
  c_boera = pylab.cm.Purples(.7)

  text_color  = 'black'
  color_line = c_pchw18
  color_becker = c_boss
  color_bosman = c_viel
  
  
  
  fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))


  ax = ax_l[0]
  for data_id in samples_tau_H_multiple:
    colormap = color_map_list[data_id]
    colors = colormap.mpl_colors
    n_colors = len( colors )
    if use_color_from_colormap: color_line = colors[n_colors//2]
    else:color_line = color_lines_list[data_id]
    label = labels_multiple[data_id]
    samples = samples_tau_H_multiple[data_id]
    z = samples['z']
    mean = samples['mean']
    high = samples['higher']
    low = samples['lower']
    if 'Highest_Likelihood' in samples:
      print( 'Plotting Highest_Likelihood T0')
      mean = samples['Highest_Likelihood']
    ax.plot( z, mean, color=color_line, zorder=1, label=label, linewidth=linewidth )
    ax.fill_between( z, high, low, color=color_line, alpha=alpha, zorder=1 )  
  

  data_set = data_optical_depth_Bosman_2018
  data_name = data_set['name']
  data_z = data_set['z']
  data_tau = data_set['tau'] 
  data_tau_sigma = data_set['tau_sigma'] 
  ax.errorbar( data_z, data_tau, yerr=data_tau_sigma, fmt='none',  alpha=0.8, ecolor= color_bosman, zorder=2)
  ax.scatter( data_z, data_tau, label=data_name, alpha=0.8, color= color_bosman, zorder=2) 

  
  data_set = data_optical_depth_Becker_2013
  data_name = data_set['name']
  data_z = data_set['z']
  z_analytical = np.linspace(2,5, 50)
  data_analytical = np.array([ Compute_analytical_TauEff_Becker(z) for z in z_analytical])
  data_tau = data_set['tau'] 
  data_tau_sigma = data_set['tau_sigma'] 
  ax.errorbar( data_z, data_tau, yerr=data_tau_sigma, fmt='none',  alpha=0.8, ecolor= color_becker, zorder=3)
  ax.scatter( data_z, data_tau, label=data_name, alpha=0.8, color= color_becker, zorder=3) 
  # ax.plot( z_analytical, data_analytical, '--', c=color_becker, zorder=4, label=data_name +  ' analytical fit')

  ax.tick_params(axis='both', which='major', direction='in', labelsize=label_size )
  ax.tick_params(axis='both', which='minor', direction='in' )
  ax.set_ylabel( r'$\tau_{eff} \,\, \mathrm{HI}$', fontsize=font_size  )
  ax.set_xlabel( r'$z$', fontsize=font_size )
  leg = ax.legend(loc=2, frameon=False, fontsize=font_size, prop=prop)
  ax.set_xlim( 2, 6 )
  ax.set_ylim( 0.1, 8)
  ax.set_yscale('log')
  

  ax = ax_l[1]
  
  for data_id in samples_tau_HeII_multiple:
    colormap = color_map_list[data_id]
    colors = colormap.mpl_colors
    n_colors = len( colors )
    if use_color_from_colormap: color_line = colors[n_colors//2]
    else:color_line = color_lines_list[data_id]
    label = labels_multiple[data_id]
    samples = samples_tau_HeII_multiple[data_id]
    z = samples['z']
    mean = samples['mean']
    high = samples['higher']
    low = samples['lower']
    if 'Highest_Likelihood' in samples:
      print( 'Plotting Highest_Likelihood T0')
      mean = samples['Highest_Likelihood']
    # ax.plot( z, mean, color=color_line, zorder=1, label=label )
    # ax.fill_between( z, high, low, color=color_line, alpha=alpha, zorder=1 )  
    sort_indices = np.argsort( z )
    z = z[sort_indices]
    mean = mean[sort_indices]
    high = high[sort_indices]
    low  = low[sort_indices]
    n_samples_intgerp = 10000
    z_interp = np.linspace( z[0], z[-1], n_samples_intgerp )  
    f_mean = interp.interp1d( z, mean, kind='cubic' )
    f_high = interp.interp1d( z, high, kind='cubic' )
    f_low  = interp.interp1d( z, low,  kind='cubic' )
    ax.plot( z_interp, f_mean(z_interp), color=color_line, zorder=1, label=label )
    ax.fill_between( z_interp, f_high(z_interp), f_low(z_interp), color=color_line, alpha=alpha, zorder=1 )  

    

  data_set = data_tau_HeII_Worserc_2019
  data_name = data_set['name']
  data_z = data_set['z']
  data_tau = data_set['tau'] 
  data_tau_sigma = data_set['tau_sigma'] 
  tau_p = data_set['tau_sigma_p']
  tau_m = data_set['tau_sigma_m']
  tau_error = [ data_tau - tau_m , tau_p - data_tau  ]
  ax.scatter( data_z, data_tau, label=data_name, alpha=0.8, color= color_becker, zorder=2) 
  ax.errorbar( data_z, data_tau, yerr=tau_error, fmt='none',  alpha=0.8, ecolor= color_becker, zorder=2)

  ax.tick_params(axis='both', which='major', direction='in', labelsize=label_size )
  ax.tick_params(axis='both', which='minor', direction='in' )
  ax.set_ylabel( r'$\tau_{eff} \,\, \mathrm{HeII}$', fontsize=font_size  )
  ax.set_xlabel( r'$z$', fontsize=font_size )
  leg = ax.legend(loc=2, frameon=False, fontsize=font_size, prop=prop)
  ax.set_xlim( 2, 3.2 )
  ax.set_ylim( 0., 8)
  # ax.set_yscale('log')


  
  figure_name = output_dir + f'fig_tau_HeII_sampling.png'
  fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
  print( f'Saved Figure: {figure_name}' )



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
  alpha = 0.4

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


def Plot_T0_Sampling( samples, output_dir, system='Shamrock', label='', plot_splines=False, multiple=False ):
   
  from scipy import interpolate as interp 
  import matplotlib
  import matplotlib.font_manager
  matplotlib.rcParams['mathtext.fontset'] = 'cm'
  matplotlib.rcParams['mathtext.rm'] = 'serif'

  if system == 'Lux':      prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/brvillas/fonts', "Helvetica.ttf"), size=12)
  if system == 'Shamrock': prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=12)

  nrows = 1
  ncols = 1
  
  if not multiple: 
    samples_multiple = {}
    samples_multiple[0] = samples
    labels_multiple = [ label ]
  else:
    samples_multiple = samples
    labels_multiple = label
  
  
  
  tick_size_major, tick_size_minor = 6, 4
  tick_label_size_major, tick_label_size_minor = 14, 12
  tick_width_major, tick_width_minor = 1.5, 1

  font_size = 18
  label_size = 16
  alpha = 0.2
  linewidth = 2

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
  
  for data_id in samples_multiple:
    colormap = color_map_list[data_id]
    colors = colormap.mpl_colors
    n_colors = len( colors )
    if use_color_from_colormap: color_line = colors[n_colors//2]
    else:color_line = color_lines_list[data_id]
    samples = samples_multiple[data_id]
    obs_name = 'T0'
    z = samples['z']
    mean = samples['mean']
    high = samples['higher']
    low = samples['lower']
    label = labels_multiple[data_id]
    if 'Highest_Likelihood' in samples:
      print( 'Plotting Highest_Likelihood T0')
      mean = samples['Highest_Likelihood']
    
    if plot_splines:
      print( '  Plotting Splines interpolation')
      n_samples_intgerp = 10000
      sort_indices = np.argsort( z )
      z = z[sort_indices]
      mean = mean[sort_indices]
      high = high[sort_indices]
      low  = low[sort_indices]
      z_interp = np.linspace( z[0], z[-1], n_samples_intgerp )  
      f_mean = interp.interp1d( z, mean, kind='cubic' )
      f_high = interp.interp1d( z, high, kind='cubic' )
      f_low  = interp.interp1d( z, low,  kind='cubic' )
      ax.plot( z_interp, f_mean(z_interp)/1e4, color=color_line, zorder=1, label=label )
      ax.fill_between( z_interp, f_high(z_interp)/1e4, f_low(z_interp)/1e4, color=color_line, alpha=alpha, zorder=1 )  
    else:
      ax.plot( z, mean/1e4, color=color_line, zorder=1, label=label, linewidth=linewidth )
      ax.fill_between( z, high/1e4, low/1e4, color=color_line, alpha=alpha, zorder=1 )  

  data_set = data_thermal_history_Gaikwad_2020a
  data_z = data_set['z']
  data_mean = data_set['T0'] 
  data_error = 0.5 * ( data_set['T0_sigma_plus'] + data_set['T0_sigma_minus'] )
  name = data_set['name']   
  ax.errorbar( data_z, data_mean/1e4, yerr=data_error/1e4, fmt='none',  alpha=0.8, ecolor= c_viel, zorder=2)
  ax.scatter( data_z, data_mean/1e4, label=name, alpha=0.8, color= c_viel, zorder=2) 
  
  data_set = data_thermal_history_Gaikwad_2020b
  data_z = data_set['z']
  data_mean = data_set['T0'] 
  data_error = 0.5 * ( data_set['T0_sigma_plus'] + data_set['T0_sigma_minus'] )
  name = data_set['name']   
  ax.errorbar( data_z, data_mean/1e4, yerr=data_error/1e4, fmt='none',  alpha=0.8, ecolor= c_boss, zorder=2)
  ax.scatter( data_z, data_mean/1e4, label=name, alpha=0.8, color= c_boss, zorder=2) 

  ax.tick_params(axis='both', which='major', direction='in', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major  )
  ax.tick_params(axis='both', which='minor', direction='in', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor  )
  ax.set_ylabel( r'$T_0   \,\,\,\, [10^4 \,\,\,\mathrm{K}\,]$', fontsize=font_size  )
  ax.set_xlabel( r'$z$', fontsize=font_size )
  leg = ax.legend(loc=1, frameon=False, fontsize=font_size, prop=prop)
  ax.set_xlim( 1.8, 8 )
  ax.set_ylim( 6000/1e4, 18000/1e4)

  figure_name = output_dir + f'fig_T0_sampling.png'
  fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
  print( f'Saved Figure: {figure_name}' )





def Plot_Comparable_Data( field, comparable_data, comparable_grid, output_dir, log_ps=False  ):

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


  if not log_ps: 
    ax.set_yscale('log')
  else:
    ax.set_ylim( -5, -1)
  ax.legend( frameon=False )

  figure_name = output_dir + 'data_for_fit.png'
  fig.savefig( figure_name, bbox_inches='tight', dpi=500 )
  print( f'Saved Figure: {figure_name}' )


def Plot_Corner( samples, data_label, labels, output_dir, n_bins_1D=20, n_bins_2D=30,  lower_mask_factor=50, multiple=False, system='Shamrock'  ):
  
  
  from scipy import interpolate as interp 
  import matplotlib
  matplotlib.rcParams['mathtext.fontset'] = 'cm'
  matplotlib.rcParams['mathtext.rm'] = 'serif'
  if system == 'Lux':      prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/brvillas/fonts', "Helvetica.ttf"), size=12)
  if system == 'Shamrock': prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=12)


  
  if not multiple:
    data_labels = [data_label]
    samples_multiple = {}
    samples_multiple[0] = samples
  else: 
    data_labels = data_label
    samples_multiple = samples
    
  samples = samples_multiple[0]  
  param_ids = samples.keys()
  n_param = len( param_ids )

      
  color = 'C0'
  data_color = 'C9'
  font_size = 16
  label_size = 26
  alpha = 0.6
  fig_size = 5
  space = 0.05
  
  n_tricks = 6
  
  tick_label_size = 14
  tick_length = 7
  tick_width = 2
  border_width = 2.0
  hist_1D_line_width = 2
  hist_1D_line_color = 'C0'
  hist_2D_colormap = palettable.cmocean.sequential.Ice_20_r.mpl_colormap
  
  color_map_0 = palettable.cmocean.sequential.Ice_20_r
  color_map_1 = palettable.cmocean.sequential.Amp_20
  color_map_2 = palettable.cmocean.sequential.Tempo_20
  color_map_3 = palettable.cmocean.sequential.Dense_20
  color_map_4 = palettable.cmocean.sequential.Algae_20
  color_map_list = [ color_map_0, color_map_1, color_map_2, color_map_3, color_map_4 ]
  

  fig, ax_l = plt.subplots(nrows=n_param, ncols=n_param, figsize=(fig_size*n_param,fig_size*n_param),  sharex='col' )
  fig.subplots_adjust( wspace=space, hspace=space )

  for j in range( n_param ):
    for i in range( n_param ):
      
      add_data_label = False
      if i== 0 and j == 0: add_data_label = True

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
      
      
      for data_id in samples_multiple:
        samples = samples_multiple[data_id]
        
        colormap = color_map_list[data_id]
        hist_2D_colormap = colormap.mpl_colormap
        contour_colormap = colormap.mpl_colors
        n_colors = len( contour_colormap )
        line_color = contour_colormap[n_colors//2]
        contour_colors = [contour_colormap[n_colors//2], contour_colormap[-1]]
        
        if plot_type == '1D':
          name  = samples[j]['name']
          trace = samples[j]['trace']
          hist, bin_edges = np.histogram( trace, bins=n_bins_1D ) 
          hist = hist / hist.sum()
          bin_centers = ( bin_edges[:-1] + bin_edges[1:] ) / 2.
          bin_width = bin_centers[0] - bin_centers[1]  
          bin_centers_interp = np.linspace( bin_centers[0], bin_centers[-1], 10000 )
          f_interp  = interp.interp1d( bin_centers, hist,  kind='cubic' )
          if add_data_label: data_label = data_labels[data_id] 
          else: data_label = ''
          ax.plot( bin_centers_interp, f_interp(bin_centers_interp),   color=line_color, linewidth=hist_1D_line_width, label=data_label  )
          # ax.plot( bin_centers, hist,   color=line_color, linewidth=hist_1D_line_width  ), 
          # ax.step( bin_centers, hist, where='mid',  color=line_color, linewidth=hist_1D_line_width  )

        if plot_type == '2D':
          trace_y = samples[j]['trace']
          trace_x = samples[i]['trace']
          hist, x_edges, y_edges = np.histogram2d( trace_x, trace_y, bins=[n_bins_2D, n_bins_2D] )
          hist = hist.astype( np.float ) 
          # hist = hist.T / hist.sum()
          hist = hist.T 
          hist_mean = hist.mean()
          hist_sigma = np.sqrt( ((hist-hist_mean)**2).mean() )
          extent = [ trace_x.min(), trace_x.max(), trace_y.min(), trace_y.max() ]
          hist_2D = hist
          lower = hist_2D.max() / lower_mask_factor
          hist_2D_masked = np.ma.masked_where( hist_2D < lower, hist_2D )
          ax.imshow( hist_2D_masked[::-1], cmap=hist_2D_colormap, extent=extent, aspect='auto', interpolation='bilinear' )
          ax.contour( hist, [hist_sigma, 2* hist_sigma], extent=extent, colors= contour_colors, linewidths=2 )
        
        if add_data_label:
          ax.legend( loc=0, frameon=False, fontsize=font_size, prop=prop )

      [sp.set_linewidth(border_width) for sp in ax.spines.values()]
      
      ax.xaxis.set_major_locator(plt.MaxNLocator(n_tricks))
      ax.yaxis.set_major_locator(plt.MaxNLocator(n_tricks))

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




def Plot_Power_Spectrum_Sampling( ps_samples, ps_data_dir, output_dir, scales='small', linewidth=2.5, system=None, name='', label='', plot_type='samples', rescaled_walther=False, rescale_walter_file=None, multiple=False ):

  
  if not multiple:
    labels_multiple = [ label ]
    ps_samples_multiple = {}
    ps_samples_multiple[0] = ps_samples
  else:
    labels_multiple = label
    ps_samples_multiple = ps_samples


  import matplotlib
  import pylab
  matplotlib.rcParams['mathtext.fontset'] = 'cm'
  matplotlib.rcParams['mathtext.rm'] = 'serif'
  if system == 'Lux':      prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/brvillas/fonts', "Helvetica.ttf"), size=12)
  if system == 'Shamrock': prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=12)


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
  
  if rescaled_walther:
    print(f" Loading Walther rescale values: {rescale_walter_file}")
    rescale_walter_alphas = Load_Pickle_Directory( rescale_walter_file)

  z_vals_small_scale  = [ 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 4.2, 4.6, 5.0, 5.4 ]
  z_vals_large_scale  = [ 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.6 ]
  z_vals_middle_scale = [ 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,   ]
  z_vals_all_scale    = z_vals_large_scale
  z_high = [ 5.0, 5.4 ]




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

  if scales == 'middle': nrows = 2


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
    c_walther = 'C1'

  alpha_bar = 0.4

  for index, current_z in enumerate( z_vals ):



    indx_j = index % ncols
    indx_i = index//ncols

    if nrows > 1: ax = ax_l[indx_i][indx_j]
    else: ax = ax_l[indx_j]


    if scales == 'middle': flags[indx_i,  indx_j] = 1

    
    if plot_type == 'multiple_lines':
      
      n_samples = len( ps_samples )
      for ps_index in range( n_samples ):
        ps_samples_local = ps_samples[ps_index]
        label = ps_samples_local['label']
        n_z_indices = len(ps_samples_local) - 1
        # print( f'N z Indices: {n_z_indices}')
        z_samples = np.array([ ps_samples_local[i]['z'] for i in range(n_z_indices) ])
        diff = np.abs( z_samples - current_z )
        diff_min = diff.min()
        if diff_min < 0.05:
          index = np.where( diff == diff_min )[0][0]
          k_vals = ps_samples_local[index]['k_vals']
          delta_mean = ps_samples_local[index]['mean']
          ax.plot( k_vals, delta_mean, linewidth=2,  zorder=1, label=label   )
          
      
    if plot_type == 'samples':
      
      for data_id in ps_samples_multiple:
        
        colormap = color_map_list[data_id]
        colors = colormap.mpl_colors
        n_colors = len( colors )
        if use_color_from_colormap: color_line = colors[n_colors//2]
        else:color_line = color_lines_list[data_id]
        
        ps_samples = ps_samples_multiple[data_id]
        label = labels_multiple[data_id]
        
        n_samples = len( ps_samples )
        z_samples = np.array([ ps_samples[i]['z'] for i in range(n_samples) ])
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
            ax.plot( k_vals, delta_mean,  linewidth=linewidth, color=color_line, zorder=1, label=label   )
            # ax.fill_between( k_vals, delta_p, delta_m, facecolor=color_line, alpha=alpha_bar, zorder=1   )
            ax.fill_between( k_vals, delta_higher, delta_lower, facecolor=color_line, alpha=alpha_bar, zorder=1   )

    
    
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
        
        if rescaled_walther and data_index in rescale_walter_alphas:
  
          rescale_z = rescale_walter_alphas[data_index]['z']
          rescale_alpha = rescale_walter_alphas[data_index]['alpha']
          print( f'  Rescaling z={rescale_z:.1f}    alpha={rescale_alpha:.3f} ')
          data_delta_power *= rescale_alpha
          d_walther = ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_walther, label=label_walther, zorder=2)

        
      # # Add Boera data
      # z_diff = np.abs( data_z_b - current_z )
      # diff_min = z_diff.min()
      # if diff_min < 1e-1:
      #   data_index = np.where( z_diff == diff_min )[0][0]
      #   data_z_local = data_z_b[data_index]
      # 
      #   data_k = data_boera[data_index]['k_vals']
      #   data_delta_power = data_boera[data_index]['delta_power']
      #   data_delta_power_error = data_boera[data_index]['delta_power_error']
      #   label_boera ='Boera et al. (2019)'
      #   d_boera = ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_boera, label=label_boera, zorder=2 )
      # 
      # 
      # # Add Viel data
      # z_diff = np.abs( data_z_v - current_z )
      # diff_min = z_diff.min()
      # if diff_min < 1e-1:
      #   data_index = np.where( z_diff == diff_min )[0][0]
      #   data_z_local = data_z_v[data_index]
      # 
      #   data_k = data_viel[data_index]['k_vals']
      #   data_delta_power = data_viel[data_index]['delta_power']
      #   data_delta_power_error = data_viel[data_index]['delta_power_error']
      #   label_viel = 'Viel et al. (2013)'
      #   d_viel = ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_viel, label=label_viel, zorder=2 )


    legend_loc = 3
    if indx_i == nrows-1 and nrows!=2: legend_loc = 2
    

    if scales == 'large': legend_loc = 2
    label_bars =  r'1$\sigma$ skewers $P\,(\Delta_F^2)$'

    add_legend = False
    if indx_j == 0: add_legend = True

    if scales == 'middle' and indx_i == nrows-1 and indx_j == ncols-1: add_legend = True

    if scales == 'middle':
      if indx_i == 0: legend_loc = 3
      if indx_i == 1: legend_loc = 3
      if indx_i == 2: legend_loc = 3
      
    if scales == 'all':
      if indx_i == 0: legend_loc = 3
      if indx_i == 1: legend_loc = 3
      if indx_i == 2: legend_loc = 3
  

    if add_legend:
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
      x_min, x_max = 2e-3, 1e-1
      if indx_i == 0: y_min, y_max = 4e-3, 9e-2
      if indx_i == 1: y_min, y_max = 1e-2, 2e-1
      
    if scales == 'all':
      x_min, x_max = 2e-3, 1e-1 
      if indx_i == 0: y_min, y_max = 4e-3, 9e-2
      if indx_i == 1: y_min, y_max = 1e-2, 2e-1
      if indx_i == 2: y_min, y_max = 3e-2, 6e-1


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


  fileName = output_dir + f'fig_flux_ps_{plot_type}_{scales}'
  if name != '': fileName += f'_{name}'
  fileName += '.png'
  # fileName += '.pdf'
  fig.savefig( fileName,  pad_inches=0.1, bbox_inches='tight', dpi=fig_dpi)
  print('Saved Image: ', fileName)