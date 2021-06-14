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
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from tools import *


def Plot_Grid_tau_HeII( sim_data_sets, output_dir, system=None, plot_splines=False ):

  from data_optical_depth_HeII import data_tau_HeII_Worserc_2019 
  
  from scipy import interpolate as interp 
  import matplotlib
  import matplotlib.font_manager
  matplotlib.rcParams['mathtext.fontset'] = 'cm'
  matplotlib.rcParams['mathtext.rm'] = 'serif'

  if system == 'Lux':      prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/brvillas/fonts', "Helvetica.ttf"), size=12)
  if system == 'Shamrock': prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=12)

  nrows = 1
  ncols = 1
  
  tick_size_major, tick_size_minor = 6, 4
  tick_label_size_major, tick_label_size_minor = 14, 12
  tick_width_major, tick_width_minor = 1.5, 1

  font_size = 18
  label_size = 16
  alpha = 0.7
  
  line_width = 0.6
  
  border_width = 1.5
    
  c_pchw18 = pylab.cm.viridis(.7)
  c_hm12 = pylab.cm.cool(.3)

  c_boss = pylab.cm.viridis(.3)
  c_walther = pylab.cm.viridis(.3)
  c_viel = 'C1'
  c_boera = pylab.cm.Purples(.7)

  text_color  = 'black'
  color_line = c_pchw18
  color_data = c_boss
  
  n_lines = len( sim_data_sets )
  colormap = palettable.cmocean.sequential.Matter_20.mpl_colormap
  colormap = palettable.cmocean.sequential.Tempo_20.mpl_colormap
  colormap = palettable.cmocean.sequential.Dense_20.mpl_colormap
  colormap = palettable.scientific.sequential.LaPaz_20_r.mpl_colormap
  colormap = palettable.colorbrewer.sequential.Blues_9_r.mpl_colormap
  colors = colormap( np.linspace(0,1,n_lines) )
  
  oranges = palettable.colorbrewer.sequential.YlOrBr_9.mpl_colors
  orange = oranges[4]
  
  black_background = True
  if black_background:
    text_color = 'white'
    color_data = orange

  fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,5*nrows))

  
  obs_name = 'tau_HeII'
  for sim_id, sim_data in enumerate(sim_data_sets):
    color_line = colors[sim_id]
    z = sim_data['z']
    mean = sim_data[obs_name]
    
    if plot_splines:
      n_samples_intgerp = 10000
      sort_indices = np.argsort( z )
      z = z[sort_indices]
      mean = mean[sort_indices]
      z_interp = np.linspace( z[0], z[-1], n_samples_intgerp )  
      f_mean = interp.interp1d( z, mean, kind='cubic' )
      
      if sim_id == n_lines-1: label = 'Simulation'
      else: label = ''
      ax.plot( z_interp, f_mean(z_interp), color=color_line, zorder=1, label=label, alpha=alpha, lw=line_width )
    else:
      ax.plot( z, vals , label=sim_data['plot_label'], color=color_line, zorder=1, alpha=alpha, lw=line_width )
      
    
  
  data_set = data_tau_HeII_Worserc_2019
  data_name = data_set['name']
  data_z = data_set['z']
  data_tau = data_set['tau'] 
  data_tau_sigma = data_set['tau_sigma'] 
  tau_p = data_set['tau_sigma_p']
  tau_m = data_set['tau_sigma_m']
  tau_error = [ data_tau - tau_m , tau_p - data_tau  ]
  ax.errorbar( data_z, data_tau, yerr=tau_error, fmt='o', color=color_data, label=data_name, zorder=2 )
  
  lower_lims = [  [3.16, 5.2] ]
  x_lenght = 0.025/4
  for lower_lim in lower_lims:
    lim_x, lim_y = lower_lim
    ax.plot( [lim_x-x_lenght, lim_x+x_lenght], [lim_y, lim_y], color=color_data,  zorder=2  )
    dx, dy = 0, 1
    ax.arrow( lim_x, lim_y, dx, dy,  color=color_data, head_width=0.01, head_length=0.08,  zorder=2   )
    

  ax.tick_params(axis='both', which='major', direction='in', color=text_color, labelcolor=text_color, labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major  )
  ax.tick_params(axis='both', which='minor', direction='in', color=text_color, labelcolor=text_color, labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor  )

  ax.set_ylabel( r'$\tau_{eff} \,\, \mathrm{HeII}$', fontsize=font_size, color=text_color  )
  ax.set_xlabel( r'$z$', fontsize=font_size, color=text_color )
  ax.set_xlim( 2, 3.4 )
  ax.set_ylim( 0, 8 )

  leg = ax.legend(loc=2, frameon=False, fontsize=22, prop=prop)
  for text in leg.get_texts():
    plt.setp(text, color = text_color)

  if black_background: 
    fig.patch.set_facecolor('black') 
    ax.set_facecolor('k')
    [ spine.set_edgecolor(text_color) for spine in list(ax.spines.values()) ]
      
  [sp.set_linewidth(border_width) for sp in ax.spines.values()]



  figure_name = output_dir + f'fig_sim_grid_tau_HeII.png'
  fig.savefig( figure_name, bbox_inches='tight', dpi=300, facecolor=fig.get_facecolor() )
  print( f'Saved Figure: {figure_name}' )




def Plot_Grid_TO_gamma( sim_data_sets, output_dir, plot_splines=True, system='Shamrock' ):
    
  from scipy import interpolate as interp 
  import matplotlib
  matplotlib.rcParams['mathtext.fontset'] = 'cm'
  matplotlib.rcParams['mathtext.rm'] = 'serif'

  if system == 'Lux':      prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/brvillas/fonts', "Helvetica.ttf"), size=12)
  if system == 'Shamrock': prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=12)

  nrows = 1
  ncols = 2

  tick_size_major, tick_size_minor = 6, 4
  tick_label_size_major, tick_label_size_minor = 14, 12
  tick_width_major, tick_width_minor = 1.5, 1

  font_size = 20
  label_size = 16
  alpha = 0.7
  
  border_width = 1.5

  c_pchw18 = pylab.cm.viridis(.7)
  c_hm12 = pylab.cm.cool(.3)

  c_boss = pylab.cm.viridis(.3)
  c_walther = pylab.cm.viridis(.3)
  c_viel = 'C1'
  c_boera = pylab.cm.Purples(.7)


  purples = palettable.colorbrewer.sequential.Purples_9.mpl_colors
  purple = purples[-1]

  oranges = palettable.colorbrewer.sequential.YlOrBr_9.mpl_colors
  orange = oranges[4]


  blues = palettable.colorbrewer.sequential.Blues_9.mpl_colors
  blue = blues[-2]

  greens = palettable.colorbrewer.sequential.BuGn_9.mpl_colors
  green = greens[-2]


  text_color  = 'black'
  color_line = c_pchw18
  color_data = c_boss

  black_background = True
  if black_background:
    text_color = 'white'
    color_data_0 = orange
    color_data_1 = 'C3'
    
    

    
  n_lines = len( sim_data_sets )
  colormap = palettable.cmocean.sequential.Matter_20.mpl_colormap
  colormap = palettable.cmocean.sequential.Tempo_20.mpl_colormap
  colormap = palettable.cmocean.sequential.Dense_20.mpl_colormap
  colormap = palettable.scientific.sequential.LaPaz_20_r.mpl_colormap
  colormap = palettable.colorbrewer.sequential.Blues_9_r.mpl_colormap
  colors = colormap( np.linspace(0,1,n_lines) )


  fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))
  plt.subplots_adjust( hspace = 0.1, wspace=0.1)
 
  ax = ax_l[0]
  for sim_id, sim_data in enumerate( sim_data_sets ):
    color_line = colors[sim_id]
    obs_name = 'T0'
    z = sim_data['z']
    mean = sim_data[obs_name]
    
    if plot_splines:
      n_samples_intgerp = 10000
      sort_indices = np.argsort( z )
      z = z[sort_indices]
      mean = mean[sort_indices]
      z_interp = np.linspace( z[0], z[-1], n_samples_intgerp )  
      f_mean = interp.interp1d( z, mean, kind='cubic' )
      
      if sim_id == n_lines-1: label = 'Simulation'
      else: label = ''
      ax.plot( z_interp, f_mean(z_interp)/1e4, color=color_line, zorder=1, label=label, alpha=alpha )
    else:
      ax.plot( z, mean/1e4, color=color_line, zorder=1, label=label )
    
  data_set = data_thermal_history_Gaikwad_2020a
  data_z = data_set['z']
  data_mean = data_set['T0'] 
  data_error = 0.5 * ( data_set['T0_sigma_plus'] + data_set['T0_sigma_minus'] )
  name = data_set['name']   
  ax.errorbar( data_z, data_mean/1e4, yerr=data_error/1e4, label=name, fmt='o', color= color_data_1, zorder=2)
  # ax.scatter( data_z, data_mean/1e4, label=name, alpha=0.8, color= color_data_1, zorder=2) 

  data_set = data_thermal_history_Gaikwad_2020b
  data_z = data_set['z']
  data_mean = data_set['T0'] 
  data_error = 0.5 * ( data_set['T0_sigma_plus'] + data_set['T0_sigma_minus'] )
  name = data_set['name']   
  ax.errorbar( data_z, data_mean/1e4, yerr=data_error/1e4, label=name, fmt='o', color= color_data_0, zorder=2)
  # ax.scatter( data_z, data_mean/1e4, label=name, alpha=0.8, color= color_data_0, zorder=2) 

  ax.tick_params(axis='both', which='major', direction='in', color=text_color, labelcolor=text_color, labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major  )
  ax.tick_params(axis='both', which='minor', direction='in', color=text_color, labelcolor=text_color, labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor  )
  ax.set_ylabel( r'$T_0   \,\,\,\, [10^4 \,\,\,\mathrm{K}\,]$', fontsize=font_size, color=text_color  )
  ax.set_xlabel( r'$z$', fontsize=font_size, color=text_color )
  ax.set_xlim( 1.8, 8 )
  ax.set_ylim( 6000/1e4, 18000/1e4)
  leg = ax.legend(loc=1, frameon=False, fontsize=22, prop=prop)
  for text in leg.get_texts():
    plt.setp(text, color = text_color)

  if black_background: 
    fig.patch.set_facecolor('black') 
    ax.set_facecolor('k')
    [ spine.set_edgecolor(text_color) for spine in list(ax.spines.values()) ]
      
  [sp.set_linewidth(border_width) for sp in ax.spines.values()]
    
    
  ax = ax_l[1]
  for sim_id, sim_data in enumerate( sim_data_sets ):
    color_line = colors[sim_id]
    obs_name = 'gamma'
    z = sim_data['z']
    mean = sim_data[obs_name] + 1
    
    if plot_splines:
      n_samples_intgerp = 10000
      sort_indices = np.argsort( z )
      z = z[sort_indices]
      mean = mean[sort_indices]
      z_interp = np.linspace( z[0], z[-1], n_samples_intgerp )  
      f_mean = interp.interp1d( z, mean, kind='cubic' )
      
      if sim_id == n_lines-1: label = 'Simulation'
      else: label = ''
      ax.plot( z_interp, f_mean(z_interp), color=color_line, zorder=1, label=label, alpha=alpha )
    else:
      ax.plot( z, mean, color=color_line, zorder=1, label=label )
    
  data_set = data_thermal_history_Gaikwad_2020a
  data_z = data_set['z']
  data_mean = data_set['gamma'] 
  data_error = 0.5 * ( data_set['gamma_sigma_plus'] + data_set['gamma_sigma_minus'] )
  name = data_set['name']   
  ax.errorbar( data_z, data_mean, yerr=data_error, label=name, fmt='o', color= color_data_1, zorder=2)
  # ax.scatter( data_z, data_mean/1e4, label=name, alpha=0.8, color= color_data_1, zorder=2) 

  data_set = data_thermal_history_Gaikwad_2020b
  data_z = data_set['z']
  data_mean = data_set['gamma'] 
  data_error = 0.5 * ( data_set['gamma_sigma_plus'] + data_set['gamma_sigma_minus'] )
  name = data_set['name']   
  ax.errorbar( data_z, data_mean, yerr=data_error, label=name, fmt='o', color= color_data_0, zorder=2)
  # ax.scatter( data_z, data_mean/1e4, label=name, alpha=0.8, color= color_data_0, zorder=2) 

  ax.tick_params(axis='both', which='major', direction='in', color=text_color, labelcolor=text_color, labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major  )
  ax.tick_params(axis='both', which='minor', direction='in', color=text_color, labelcolor=text_color, labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor  )
  ax.set_ylabel( r'$\gamma$', fontsize=font_size, color=text_color  )
  ax.set_xlabel( r'$z$', fontsize=font_size, color=text_color )
  ax.set_xlim( 1.8, 8 )
  # ax.set_ylim( 6000/1e4, 18000/1e4)
  leg = ax.legend(loc=1, frameon=False, fontsize=22, prop=prop)
  for text in leg.get_texts():
    plt.setp(text, color = text_color)

  if black_background: 
    fig.patch.set_facecolor('black') 
    ax.set_facecolor('k')
    [ spine.set_edgecolor(text_color) for spine in list(ax.spines.values()) ]
      
  [sp.set_linewidth(border_width) for sp in ax.spines.values()]
    
    
    
  figure_name = output_dir + f'fig_sim_grid_T0_gamma.png'
  fig.savefig( figure_name, bbox_inches='tight', dpi=300, facecolor=fig.get_facecolor() )
  print( f'Saved Figure: {figure_name}' )
  

def Plot_Grid_TO( sim_data_sets, output_dir, plot_splines=True, system='Shamrock' ):
    
  from scipy import interpolate as interp 
  import matplotlib
  matplotlib.rcParams['mathtext.fontset'] = 'cm'
  matplotlib.rcParams['mathtext.rm'] = 'serif'

  if system == 'Lux':      prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/brvillas/fonts', "Helvetica.ttf"), size=12)
  if system == 'Shamrock': prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=12)

  nrows = 1
  ncols = 1

  tick_size_major, tick_size_minor = 6, 4
  tick_label_size_major, tick_label_size_minor = 14, 12
  tick_width_major, tick_width_minor = 1.5, 1

  font_size = 20
  label_size = 16
  alpha = 0.7
  
  border_width = 1.5

  c_pchw18 = pylab.cm.viridis(.7)
  c_hm12 = pylab.cm.cool(.3)

  c_boss = pylab.cm.viridis(.3)
  c_walther = pylab.cm.viridis(.3)
  c_viel = 'C1'
  c_boera = pylab.cm.Purples(.7)


  purples = palettable.colorbrewer.sequential.Purples_9.mpl_colors
  purple = purples[-1]

  oranges = palettable.colorbrewer.sequential.YlOrBr_9.mpl_colors
  orange = oranges[4]


  blues = palettable.colorbrewer.sequential.Blues_9.mpl_colors
  blue = blues[-2]

  greens = palettable.colorbrewer.sequential.BuGn_9.mpl_colors
  green = greens[-2]


  text_color  = 'black'
  color_line = c_pchw18
  color_data = c_boss

  black_background = True
  if black_background:
    text_color = 'white'
    color_data_0 = orange
    color_data_1 = 'C3'
    
    

    
  n_lines = len( sim_data_sets )
  colormap = palettable.cmocean.sequential.Matter_20.mpl_colormap
  colormap = palettable.cmocean.sequential.Tempo_20.mpl_colormap
  colormap = palettable.cmocean.sequential.Dense_20.mpl_colormap
  colormap = palettable.scientific.sequential.LaPaz_20_r.mpl_colormap
  colormap = palettable.colorbrewer.sequential.Blues_9_r.mpl_colormap
  colors = colormap( np.linspace(0,1,n_lines) )


  fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))

  for sim_id, sim_data in enumerate( sim_data_sets ):

    color_line = colors[sim_id]
    obs_name = 'T0'
    z = sim_data['z']
    mean = sim_data[obs_name]
    
    if plot_splines:
      n_samples_intgerp = 10000
      sort_indices = np.argsort( z )
      z = z[sort_indices]
      mean = mean[sort_indices]
      z_interp = np.linspace( z[0], z[-1], n_samples_intgerp )  
      f_mean = interp.interp1d( z, mean, kind='cubic' )
      
      if sim_id == n_lines-1: label = 'Simulation'
      else: label = ''
      ax.plot( z_interp, f_mean(z_interp)/1e4, color=color_line, zorder=1, label=label, alpha=alpha )
    else:
      ax.plot( z, mean/1e4, color=color_line, zorder=1, label=label )
    
  data_set = data_thermal_history_Gaikwad_2020a
  data_z = data_set['z']
  data_mean = data_set['T0'] 
  data_error = 0.5 * ( data_set['T0_sigma_plus'] + data_set['T0_sigma_minus'] )
  name = data_set['name']   
  ax.errorbar( data_z, data_mean/1e4, yerr=data_error/1e4, label=name, fmt='o', color= color_data_1, zorder=2)
  # ax.scatter( data_z, data_mean/1e4, label=name, alpha=0.8, color= color_data_1, zorder=2) 

  data_set = data_thermal_history_Gaikwad_2020b
  data_z = data_set['z']
  data_mean = data_set['T0'] 
  data_error = 0.5 * ( data_set['T0_sigma_plus'] + data_set['T0_sigma_minus'] )
  name = data_set['name']   
  ax.errorbar( data_z, data_mean/1e4, yerr=data_error/1e4, label=name, fmt='o', color= color_data_0, zorder=2)
  # ax.scatter( data_z, data_mean/1e4, label=name, alpha=0.8, color= color_data_0, zorder=2) 

  ax.tick_params(axis='both', which='major', direction='in', color=text_color, labelcolor=text_color, labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major  )
  ax.tick_params(axis='both', which='minor', direction='in', color=text_color, labelcolor=text_color, labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor  )
  ax.set_ylabel( r'$T_0   \,\,\,\, [10^4 \,\,\,\mathrm{K}\,]$', fontsize=font_size, color=text_color  )
  ax.set_xlabel( r'$z$', fontsize=font_size, color=text_color )
  ax.set_xlim( 1.8, 8 )
  ax.set_ylim( 6000/1e4, 18000/1e4)
  leg = ax.legend(loc=1, frameon=False, fontsize=22, prop=prop)
  for text in leg.get_texts():
    plt.setp(text, color = text_color)

  if black_background: 
    fig.patch.set_facecolor('black') 
    ax.set_facecolor('k')
    [ spine.set_edgecolor(text_color) for spine in list(ax.spines.values()) ]
      
  [sp.set_linewidth(border_width) for sp in ax.spines.values()]
    
  figure_name = output_dir + f'fig_sim_grid_T0.png'
  fig.savefig( figure_name, bbox_inches='tight', dpi=300, facecolor=fig.get_facecolor() )
  print( f'Saved Figure: {figure_name}' )