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
from load_tabulated_data import load_power_spectrum_table, load_tabulated_data_boera, load_tabulated_data_viel, load_data_boss


import matplotlib
import matplotlib.font_manager

# print(matplotlib.font_manager.get_cachedir())




# print(matplotlib.font_manager.findSystemFonts(fontpaths='/home/bruno/fonts/', fontext='ttf'))
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'

uvb = 'pchw18'
# uvb = 'hm12'
dataDir = '/home/bruno/Desktop/data/'
# dataDir = '/raid/bruno/data/'
# dataDir = '/data/groups/comp-astro/bruno/'
simulation_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/'
output_dir = simulation_dir + 'figures/flux_power_spectrum/'.format(uvb)
create_directory( output_dir )

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

prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=12)

errorbar = True

plot_boss = False
plot_boss = True

#Cosmological Parameters 
H0 = 67.66 
cosmo_h = H0 / 100
Omega_M = 0.3111
Omega_L = 0.6889

#Box parameters
Lbox = 50.0 #Mpc/h
nPoints = 2048
nx = nPoints
ny = nPoints
nz = nPoints
ncells = nx * ny * nz


dir_boss = 'data/data_power_spectrum_boss/'
data_filename = dir_boss + 'data_table.py'
data_boss = load_data_boss( data_filename )
data_z_boss = data_boss['z_vals']

data_filename = 'data/data_power_spectrum_walther_2019/data_table.txt'
data_walther = load_power_spectrum_table( data_filename )
data_z_w = data_walther['z_vals']

dir_data_boera = 'data/data_power_spectrum_boera_2019/'
data_boera = load_tabulated_data_boera( dir_data_boera )
data_z_b = data_boera['z_vals']

data_dir_viel = 'data/data_power_spectrum_viel_2013/'
data_viel = load_tabulated_data_viel( data_dir_viel)
data_z_v = data_viel['z_vals']


snapshots_indices = [ 83, 90,  96, 102,  119, 124, 130, 136, 143, 151, 159, 169, ]
if plot_boss: snapshots_indices = [  96,  102, 106, 110,  114, 119, 124, 130, 136, 143, 151, 159 ]
snapshots_indices.reverse()

uvb_list = ['pchw18', 'hm12']

n_in_sample = 5000

factor = 5

data_cholla = {}
for uvb in uvb_list:
  input_dir = simulation_dir + 'transmited_flux_{0}_review/bootstraped_power_spectrum/statistics/'.format(uvb)
  data_cholla[uvb] = {}
  for n_snap in snapshots_indices:
    stats_file = input_dir + f'stats_{n_snap}.pkl'
    print( f'Loading File: {stats_file}')
    data = pickle.load( open( stats_file, 'rb' ) ) 
    data_cholla[uvb][n_snap] = {}
    data_cholla[uvb][n_snap]['current_z'] = data['current_z']
    data_cholla[uvb][n_snap]['k_vals'] = data['k_vals']
    data_cholla[uvb][n_snap]['mean'] = data['mean'] * data_cholla[uvb][n_snap]['k_vals'] / np.pi
    data_cholla[uvb][n_snap]['sigma'] = data['bootstrap'][n_in_sample]['statistics']['sigma'] * data_cholla[uvb][n_snap]['k_vals'] / np.pi
    data_cholla[uvb][n_snap]['sigma'] *= factor


nrows = 3
ncols = 4
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(2*fig_width,5*nrows))
plt.subplots_adjust( hspace = 0.02, wspace=0.02)


c_pchw18 = pylab.cm.viridis(.7)
c_hm12 = pylab.cm.cool(.3)

c_boss = pylab.cm.viridis(.3)
c_walther = pylab.cm.viridis(.3)
c_viel = 'C1'
c_boera = pylab.cm.Purples(.7)

alpha_bar = 0.4
text_color  = 'black'





for uvb_index,uvb in enumerate(uvb_list):



  plot_data_observed = False
  if uvb_index == len(uvb_list)-1: plot_data_observed = True


  if uvb == 'pchw18':
    color_line = c_pchw18
    label = 'CHIPS.P19'

  if uvb == 'hm12':
    color_line = c_hm12
    label = 'CHIPS.HM12'


  print( f'Plotting UVB: {uvb}' )
  data = data_cholla[uvb]
  for snap_index, nSnap in enumerate(snapshots_indices):



    indx_j = snap_index % ncols
    indx_i = snap_index//ncols

    ax = ax_l[indx_i][indx_j]


    factor = 1.0
    if indx_i == nrows-1: factor = 1.1

    if plot_boss: 
      factor = 1.1
      if indx_i == nrows-1 and indx_j==ncols-1: factor = 1.0

    current_z = data[nSnap]['current_z']
    k = data[nSnap]['k_vals']
    delta = data[nSnap]['mean'] * factor
    sigma = data[nSnap]['sigma']
    delta_p = delta + sigma
    delta_m = delta - sigma
    if uvb_index == 0: line_pchw18, = ax.plot( k, delta, c=color_line, linewidth=3  )
    if uvb_index == 1: line_hm12, = ax.plot( k, delta, c=color_line, linewidth=3  )
    if errorbar: 
      if uvb_index == 0: bar_pchw18 = ax.fill_between( k, delta_p, delta_m, facecolor=color_line, alpha=alpha_bar  )
      if uvb_index == 1: bar_hm12 = ax.fill_between( k, delta_p, delta_m, facecolor=color_line, alpha=alpha_bar  )




    ax.text(0.85, 0.95, r'$z={0:.1f}$'.format(current_z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size, color=text_color) 


    if plot_data_observed :

      if plot_boss:

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
          d_boss = ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_boss, )

      else:

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
          d_walther = ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_walther, )


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
          d_boera = ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_boera,  )


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
          d_viel = ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_viel, )


    legend_loc = 3
    if indx_i == nrows-1 and nrows!=2: legend_loc = 2

    if plot_boss: legend_loc = 2
    label_bars =  r'1$\sigma$ skewers $P\,(\Delta_F^2)$'

    if indx_j == 0 and uvb_index == 1 :
      # leg = ax.legend( loc=legend_loc, frameon=False, fontsize=12)
      if plot_boss: leg = ax.legend( [line_pchw18, line_hm12, d_boss], ['CHIPS.P19', 'CHIPS.HM12',  label_boss ], loc=legend_loc, frameon=False, prop=prop,   handler_map={tuple: HandlerTuple(ndivide=None)}  )
      else: 
        if indx_i in [0, 1]: leg = ax.legend( [line_pchw18, line_hm12, d_walther], ['CHIPS.P19', 'CHIPS.HM12', label_walther ], loc=legend_loc, frameon=False, prop=prop,   handler_map={tuple: HandlerTuple(ndivide=None)}  )
        else: leg = ax.legend( [line_pchw18, line_hm12, d_boera, d_viel], ['CHIPS.P19', 'CHIPS.HM12', label_boera, label_viel ], loc=legend_loc, frameon=False, prop=prop,   handler_map={tuple: HandlerTuple(ndivide=None)}  )


      for text in leg.get_texts():
          plt.setp(text, color = text_color)



    x_min, x_max = 4e-3, 2.5e-1
    if indx_i == 0: y_min, y_max = 1e-3, 9e-2
    if indx_i == 1: y_min, y_max = 5e-3, 2e-1
    if indx_i == 2: y_min, y_max = 5e-2, 3

    if plot_boss:
      x_min, x_max = 2e-3, 2.3e-2
      if indx_i == 0: y_min, y_max = 1e-2, 1.2e-1
      if indx_i == 1: y_min, y_max = 2e-2, 2.5e-1
      if indx_i == 2: y_min, y_max = 5e-2, 7e-1




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











fileName = output_dir + 'flux_power_spectrum_grid_review'
if plot_boss: fileName += '_BOSS'
fileName += '.png'
# fileName += '.pdf'
fig.savefig( fileName,  pad_inches=0.1, bbox_inches='tight', dpi=fig_dpi)
print('Saved Image: ', fileName)


