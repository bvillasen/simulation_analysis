import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import palettable
import pylab
import pickle
from matplotlib.legend_handler import HandlerTuple
import os, sys
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from tools import *
from load_tabulated_data import load_power_spectrum_table, load_tabulated_data_boera, load_tabulated_data_viel, load_data_boss, load_data_irsic
from plotting_tools import smooth_line


data_dir = '/raid/bruno/data/'
input_dir_0  = data_dir + f'cosmo_sims/rescaled_P19/1024_50Mpc/analysis_files/ps_statistics/'
input_dir_1  = data_dir + f'cosmo_sims/rescaled_P19/2048_100Mpc/analysis_files/ps_statistics/'
input_dir_2  = data_dir + f'cosmo_sims/rescaled_P19/2048_200Mpc/analysis_files/ps_statistics/'
output_dir = data_dir + f'cosmo_sims/rescaled_P19/figures/'
create_directory( output_dir )

input_dir_list = [ input_dir_0, input_dir_1, input_dir_2 ]
input_dir_list = [ input_dir_0, input_dir_1 ]

line_widths = [ 2, 2, 2 ]
line_styles = [ 'solid', 'dashed', 'dashed' ]

show_error_bar = [ True, True, False ]

labels = [ r'$50 \,\,\,\, h^{-1}\mathrm{Mpc} \,\,\, 1024^3 $', r'$100 \,h^{-1}\mathrm{Mpc} \,\,\, 2048^3 $', r'$200 \,h^{-1}\mathrm{Mpc} \,\,\, 2048^3 $']


file_ids = range(15, 56)
data_sets = []
for data_id, input_dir in  enumerate(input_dir_list) :
  data_set = {}
  z_vals = []
  for index,file_id in enumerate(file_ids):
    in_file_name = input_dir + f'{file_id}_stats.h5'
    print( f'Loading File: {in_file_name} ' )
    in_file = h5.File( in_file_name, 'r' )
    
    data_set[index] = {}
    for key in in_file.attrs.keys():
      data_set[index][key] = in_file.attrs[key][0]
    z_vals.append( in_file.attrs['current_z'][0] )
    k_vals = in_file['k_vals'][...]
    ps_mean = in_file['mean'][...]
    ps_max  = in_file['max'][...]
    ps_high = in_file['higher'][...]
    ps_low  = in_file['lower'][...]
    n_independent = in_file['n_independent'][...]
    in_file.close()
    data_set[index]['n_independent'] = n_independent 
    data_set[index]['k_vals'] = k_vals 
    data_set[index]['mean'] = k_vals / np.pi * ps_mean
    data_set[index]['max']  = k_vals / np.pi * ps_max
    data_set[index]['high'] = k_vals / np.pi * ps_high
    data_set[index]['low']  = k_vals / np.pi * ps_low
  z_vals = np.array( z_vals )
  data_set['z'] = z_vals
  data_set['label'] = labels[data_id]
  data_sets.append( data_set )


 
plot_ps_data = data_sets
scales = 'middle'
scales = 'all'
# scales = 'large'
  
high_redshift = True
# high_redshift = False


dir_boss = 'data/data_power_spectrum_boss/'
data_filename = dir_boss + 'data_table.py'
data_boss = load_data_boss( data_filename )
data_z_boss = data_boss['z_vals']

dir_irsic = 'data/data_power_spectrum_irsic_2017/'
data_filename = dir_irsic + 'data_table.py'
data_irsic = load_data_irsic( data_filename )
data_z_irsic = data_irsic['z_vals']

data_filename = 'data/data_power_spectrum_walther_2019/data_table.txt'
data_walther = load_power_spectrum_table( data_filename )
data_z_w = data_walther['z_vals']

dir_data_boera = 'data/data_power_spectrum_boera_2019/'
data_boera = load_tabulated_data_boera( dir_data_boera )
data_z_b = data_boera['z_vals']

data_dir_viel = 'data/data_power_spectrum_viel_2013/'
data_viel = load_tabulated_data_viel( data_dir_viel)
data_z_v = data_viel['z_vals']



z_vals_small_scale  = [ 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 4.2, 4.6, 5.0, 5.4 ]
z_vals_large_scale  = [ 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4 ]
z_vals_middle_scale = [ 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2   ]
z_vals_all_scale    = [ 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4 ]
z_high = [ 4.6, 4.8, 5.0, 5.4 ]
z_high = [ 4.2, 4.6, 5.0  ]

if scales == 'large':    z_vals = z_vals_large_scale
elif scales == 'small':  z_vals = z_vals_small_scale
elif scales == 'middle': z_vals = z_vals_middle_scale
elif scales == 'all':    z_vals = z_vals_all_scale
else: 
  print( "ERROR: Scales = large,  small, middle or all ")
  # return
if high_redshift: z_vals = z_high

plot_boss, plot_boera, plot_viel, plot_walther, plot_irsic = False, False, False, False, False 
if scales == 'large' : plot_boss = True

# plot_boss = plot_irsic = True
plot_boss = True 
plot_boera = True

fig_width = 4
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
text_color  = 'black'
linewidth = 2
alpha_bar = 0.5

# Colors
bright_green = pylab.cm.viridis(.7)
light_blue = pylab.cm.cool(.3)
dark_blue = pylab.cm.viridis(.3) 
purple = pylab.cm.Purples(.7)
blue = 'C0'
orange = 'C1'
green = 'C2'
red = 'C3'
purple_2 = 'C4'

c_boss = dark_blue
c_boera = purple_2
c_irsic = purple_2


colors = [ blue, orange, green ]

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'
prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=12)
  
nrows, ncols = 3, 4
if high_redshift: nrows, ncols = 1, 3
if scales == 'middle': nrows, ncols = 2, 4
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_width*ncols,5*nrows))
plt.subplots_adjust( hspace = 0.02, wspace=0.02)


n_data_sets = len( plot_ps_data )

for index, current_z in enumerate( z_vals ):



  indx_j = index % ncols
  indx_i = index//ncols
  if nrows > 1: ax = ax_l[indx_i][indx_j]
  else: ax = ax_l[indx_j]
  
  for data_id in range( n_data_sets ):
    data_set = plot_ps_data[data_id]
    z_vals = data_set['z']
    diff = np.abs( z_vals - current_z )
    diff_min = diff.min()
    color_line = colors[data_id]
    if diff_min < 0.05:
      index = np.where( diff == diff_min )[0][0]
      data_ps = data_set[index]
      k_vals = data_ps['k_vals']
      delta_ps = data_ps['mean']
      label = data_set['label']
      linewidth = line_widths[data_id]
      linestyle = line_styles[data_id]
      ax.plot( k_vals, delta_ps,  linewidth=linewidth, color=color_line, zorder=1, label=label, linestyle=linestyle   )
      if show_error_bar[data_id]:
        high = data_ps['high']
        low = data_ps['low']
        n_ind = data_ps['n_independent']
        sigma  = 0.5*(  high - low )
        high_ind = delta_ps + 1/np.sqrt(n_ind) * ( high - delta_ps ) 
        low_ind =  delta_ps - 1/np.sqrt(n_ind) * (  delta_ps - low ) 
        if data_id == 0: high_ind[4]*=1.1
        if data_id == 1: high_ind[4]*=1.1
        ax.fill_between( k_vals, high_ind, low_ind, facecolor=color_line, alpha=alpha_bar, zorder=1  )
        
    
  
  ax.text(0.85, 0.95, r'$z={0:.1f}$'.format(current_z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size, color=text_color) 

  # Add Boss data
  if plot_boss:
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

  # Add Irsic data
  if plot_irsic:
    z_diff = np.abs( data_z_irsic - current_z )
    diff_min = z_diff.min()
    if diff_min < 1e-1:
      data_index = np.where( z_diff == diff_min )[0][0]
      data_z_local = data_z_irsic[data_index]

      data_k = data_irsic[data_index]['k_vals']
      data_delta_power = data_irsic[data_index]['delta_power']
      data_delta_power_error = data_irsic[data_index]['delta_power_error']
      label_irsic = 'Irsic et al. (2017)'
      d_irsic = ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_irsic, label=label_irsic, zorder=2)

  # Add Walther data
  if plot_walther:
    z_diff = np.abs( data_z_w - current_z )
    diff_min = z_diff.min()
    if diff_min < 1e-1:
      data_index = np.where( z_diff == diff_min )[0][0]
      data_z_local = data_z_w[data_index]

      data_k = data_walther[data_index]['k_vals']
      data_delta_power = data_walther[data_index]['delta_power']
      data_delta_power_error = data_walther[data_index]['delta_power_error']
      label_walther ='Walther et al. (2018)' 


  # Add Boera data
  if plot_boera:
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
  if plot_viel:
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


  if scales == 'large': legend_loc = 2
  if scales == 'all':   legend_loc = 3
  if scales == 'middle':   legend_loc = 3
  
    
  if indx_j == 0:
    ax.legend( loc=legend_loc, frameon=False, prop=prop )
  
  if scales == 'small':
    x_min, x_max = 4e-3, 2.5e-1
    if indx_i == 0: y_min, y_max = 1e-3, 9e-2
    if indx_i == 1: y_min, y_max = 5e-3, 2e-1
    if indx_i == 2: y_min, y_max = 5e-2, 3

  if scales == 'large':
    x_min, x_max = 8e-4, 2.3e-2
    if indx_i == 0: y_min, y_max = 6e-3, 9e-2
    if indx_i == 1: y_min, y_max = 1.2e-2, 2.05e-1
    if indx_i == 2: y_min, y_max = 2.5e-2, 6e-1

  if scales == 'middle':
    x_min, x_max = 2e-3, 7e-2
    if indx_i == 0: y_min, y_max = 1e-2, 2e-1
    if indx_i == 1: y_min, y_max = 3e-2, 5e-1
    
  if scales == 'all':
    x_min, x_max = 2e-4, 1e-0 
    x_min, x_max = 1e-3, 3e-1 
    if indx_i == 0: y_min, y_max = 3e-7, 9e-2
    if indx_i == 1: y_min, y_max = 1e-6, 2e-1
    if indx_i == 2: y_min, y_max = 1e-5, 6e-1
    
  if high_redshift:
    if indx_i == 0: y_min, y_max = 5e-3, 8e-1
  
    


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


file_name = output_dir + f'flux_power_spectrum_grid_{scales}'
if high_redshift: file_name += '_highZ'
file_name += '.png'

fig.savefig( file_name,  pad_inches=0.1, bbox_inches='tight', dpi=fig_dpi)
print('Saved Image: ', file_name )


   