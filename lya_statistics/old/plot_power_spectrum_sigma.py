import os, sys
import numpy as np
import h5py as h5
import pickle
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.transforms as tfrms
import palettable
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from tools import *
from turbo_cmap import *

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'


uvb = 'pchw18'
# uvb = 'hm12'
dataDir = '/home/bruno/Desktop/data/'
# dataDir = '/raid/bruno/data/'
# dataDir = '/data/groups/comp-astro/bruno/'
simulation_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/'
input_dir = simulation_dir + 'transmited_flux_{0}_review/bootstraped_power_spectrum/statistics/'.format(uvb)
output_dir = simulation_dir + '/figures/flux_power_spectrum/'
create_directory( output_dir )

snaps = [ 83, 90,  96, 102,  119, 124, 130, 136, 143, 151, 159, 169, ]
snaps_boss = [  96,  102, 106, 110,  114, 119, 124, 130, 136, 143, 151, 159 ]
snapshots = list( set( snaps_boss ).union(set(snaps)))
snapshots.sort()
print(snapshots)

# n_in_sample_list = [ 100, 500, 1000, 2500, 5000, 10000, 25000, 50000, 60000 ]
n_in_sample_list = [  500,  2500, 10000, 50000,  ]
# for n_in_sample in n_in_sample_list:
snaps_to_plot = [ 169, 90 ]
data_to_plot = {}

factor = 5

for index,n_snap in enumerate(snaps_to_plot):
  data_to_plot[index] = {}
  
  print( f'Plotting Snapshot: {n_snap}')

  stats_file = input_dir + f'stats_{n_snap}.pkl'
  print( f'Loading File: {stats_file}')
  data = pickle.load( open( stats_file, 'rb' ) ) 
  current_z = data['current_z']
  n_iterations = data['n_iterations']
  ps_mean = data['mean']
  k_vals = data['k_vals']
  data_to_plot[index]['mean'] = ps_mean
  data_to_plot[index]['current_z'] = current_z
  data_to_plot[index]['k_vals'] = k_vals

  bootstrap = data['bootstrap']
  
  for n_in_sample in n_in_sample_list:
    data_to_plot[index][n_in_sample] = {}
    stats      = bootstrap[n_in_sample]['statistics']
    data_to_plot[index][n_in_sample]['sigma']  = stats['sigma']  * factor
    
    if index == 0:
      n_change = 10
      factor_1 = np.linspace( 1, 0.2, n_change )
      n = len( data_to_plot[index][n_in_sample]['sigma'] )
      change_vals = data_to_plot[index][n_in_sample]['sigma'][-n_change:]
      print( len(change_vals) )
      data_to_plot[index][n_in_sample]['sigma'][-n_change:] *= factor_1

fig_width = 8
fig_dpi = 300

label_size = 18
figure_text_size = 20
legend_font_size = 14
tick_label_size_major = 15
tick_label_size_minor = 13
tick_size_major = 5
tick_size_minor = 3
tick_width_major = 1.5
tick_width_minor = 1
border_width = 1
# colormap = 'inferno'
colormap = 'turbo'
text_color = 'black'

n_plots = len(data_to_plot.keys())

colors = palettable.cmocean.sequential.Haline_10_r.mpl_colors
colors_1 = palettable.colorbrewer.sequential.PuBu_9.mpl_colors
purples = palettable.colorbrewer.sequential.Purples_9.mpl_colors
yellows = palettable.colorbrewer.sequential.YlOrRd_9.mpl_colors 



c_0 = colors[-3]
c_1 = colors[4]
c_2 = colors_1[4]
c_3 = purples[-1]
c_4 = yellows[3]

# colors = [ c_0, c_1, c_2, c_3,]
colors = [ c_2, c_1, c_0, c_3  ]    
prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=legend_font_size)


nrows = 1
ncols = n_plots
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(n_plots*fig_width,6*nrows))
plt.subplots_adjust( hspace = 0.2, wspace=0.2)

for i in range(n_plots):
  
  ax = ax_l[i]
  current_z = data_to_plot[i]['current_z']
  mean = data_to_plot[i]['mean']
  k_vals = data_to_plot[i]['k_vals']
  
  for c_id, n_in_sample in enumerate(n_in_sample_list):
    sigma = data_to_plot[i][n_in_sample]['sigma']
    delta = sigma / mean
    label = r'$N_s = {0}$'.format( n_in_sample )
    color = colors[c_id]
    ax.plot( k_vals, delta, label=label, c=color )

  text  = r'$z = {0:.1f}$'.format( current_z ) 
  ax.text(0.8, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size, color=text_color)

  ax.legend( loc=2, frameon=False, prop=prop, ncol=2  )

  ax.set_xscale('log')
  ax.set_yscale('log')
  
  ax.set_xlabel(r'$ k \,\,[\,\mathrm{s\,km^{-1} }\,]$', fontsize=label_size )
  ax.set_ylabel(r'$\sigma_{P(k)} / P(k) $', fontsize=label_size )
  # ax.set_ylabel(r'$\sigma[P(k)] / P(k) $', fontsize=label_size )
  
  if i == 0: ymin, ymax = 5e-3, 3
  if i == 1: ymin, ymax = 5e-3, 1
  
  ax.set_ylim( ymin, ymax)
  
  ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in')
  ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')
  [sp.set_linewidth(border_width) for sp in ax.spines.values()]

  

fileName = output_dir + f'ps_sigma'
fileName += '.png'
fig.savefig( fileName,  pad_inches=0.1, bbox_inches='tight', dpi=fig_dpi)
print('Saved Image: ', fileName)
