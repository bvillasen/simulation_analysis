import sys, os
import numpy as np
import h5py as h5
import palettable
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import pylab
import pickle
from matplotlib.legend_handler import HandlerTuple
import os, sys
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from tools import *
from load_tabulated_data import load_power_spectrum_table, load_tabulated_data_boera, load_tabulated_data_viel, load_data_boss


import matplotlib
import matplotlib.font_manager
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'


# print(matplotlib.font_manager.get_cachedir())




def plot_power_spectrum_grid( ps_data_dir, output_dir, scales='large', sim_data_sets=None, system=None, high_z_only=False  ):
  
  if system == 'Lux' or system == 'Summit': matplotlib.use('Agg')
  import matplotlib.pyplot as plt


  fig_height = 5
  fig_width = 8
  fig_dpi = 300
  
  if high_z_only: fig_height = 8

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
  z_high = [ 5.0, 5.4 ]
  
  
  if scales == 'large': z_vals = z_vals_large_scale
  elif scales == 'small': z_vals = z_vals_small_scale
  elif scales == 'middle': z_vals = z_vals_middle_scale
  else: 
    print( "ERROR: Scales = large,  small of middle ")
    return
    
  if high_z_only: z_vals = z_high
    
  nrows = 3
  ncols = 4
  
  if high_z_only:    nrows, ncols = 1, 2
  
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
  
  if scales == 'middle':
    c_walther = 'C3'
  

  for index, current_z in enumerate( z_vals ):



    indx_j = index % ncols
    indx_i = index//ncols

    if nrows > 1: ax = ax_l[indx_i][indx_j]
    else: ax = ax_l[indx_j]
    
    
    if scales == 'middle': flags[indx_i,  indx_j] = 1

    if sim_data_sets:
      for sim_data in sim_data_sets:
        sim_z_vals = sim_data['z']
        diff = np.abs( sim_z_vals - current_z )
        diff_min = diff.min()
        index = np.where( diff == diff_min )[0][0]
        # print( index )
        if diff_min < 0.08:
          k = sim_data['ps_kvals'][index]
          ps = sim_data['ps_mean'][index]
          delta = ps * k / np.pi 
          ax.plot( k, delta, linewidth=3, label=sim_data['plot_label'], zorder=1  )
          # ax.plot( k, delta, c=color_line, linewidth=3, label=sim_data['plot_label']  )
          

    ax.text(0.85, 0.95, r'$z={0:.1f}$'.format(current_z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size, color=text_color) 

    
    if scales == 'large' or scales == 'middle':

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

    if scales == 'small' or scales == 'middle':
      
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
      
    if high_z_only: y_min, y_max = 5e-2, 3
      
      



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
        

  fileName = output_dir + f'flux_ps_grid_{scales}'
  if high_z_only: fileName += '_highZ'
  fileName += '.png'
  # fileName += '.pdf'
  fig.savefig( fileName,  pad_inches=0.1, bbox_inches='tight', dpi=fig_dpi)
  print('Saved Image: ', fileName)


