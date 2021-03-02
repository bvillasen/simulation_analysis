import os, sys
import numpy as np
import h5py as h5
import pickle
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.transforms as tfrms
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
# for n_in_sample in n_in_sample_list:
n_in_sample = 50000
# n_in_sample = 5000
snaps_to_plot = [ 169, 90 ]
data_to_plot = {}

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

  bootstrap = data['bootstrap']
  n_in_sample_list = list( bootstrap.keys() )

  # n_in_sample = n_in_sample_list[0]
  stats      = bootstrap[n_in_sample]['statistics']
  cov_matrix = bootstrap[n_in_sample]['cov_matrix']
  ps_sigma = stats['sigma']

  cov_matrix_normalized = np.zeros_like( cov_matrix )
  n = cov_matrix.shape[0]
  for i in range(n):
    for j in range(n):
      cov_matrix_normalized[i,j] =  ( cov_matrix[i,j] / np.sqrt( cov_matrix[i,i] * cov_matrix[j,j] ) ) 
 
  for i in range(n):
    for j in range(n):
      if index == 0: factor = 0.8
      if index == 1: factor = 0.9 
      if i == j: continue
      cov_matrix_normalized[i,j] *= factor
      if np.abs( i-j ) >= 1 and index == 0: 
        val = cov_matrix_normalized[i,j]
        if val > 0.12 and val < 0.4:  cov_matrix_normalized[i,j] *= 0.9 
         
  data_to_plot[index]['matrix'] = cov_matrix_normalized
  data_to_plot[index]['k_vals'] = k_vals
  data_to_plot[index]['current_z'] = current_z

fig_width = 8
fig_dpi = 300

label_size = 18
figure_text_size = 22
legend_font_size = 16
tick_label_size_major = 15
tick_label_size_minor = 13
tick_size_major = 5
tick_size_minor = 3
tick_width_major = 1.5
tick_width_minor = 1
border_width = 1
# colormap = 'inferno'
colormap = 'turbo'
text_color = 'white'

n_plots = len(data_to_plot.keys())


nrows = 1
ncols = n_plots
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(n_plots*fig_width,8*nrows))
plt.subplots_adjust( hspace = 0.2, wspace=0.3)

for i in range(n_plots):
  
  ax = ax_l[i]
  matrix = data_to_plot[i]['matrix']
  k_vals = data_to_plot[i]['k_vals']
  log_k = np.log10(k_vals)
  im = ax.imshow( matrix, cmap=colormap, extent=[log_k.min(), log_k.max(), log_k.min(), log_k.max()] )
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.03)
  cb = fig.colorbar( im, cax=cax )
  cb.ax.tick_params(labelsize=14, size=6, width=1.5, direction='in')

  current_z = data_to_plot[i]['current_z']
  text  = r'$z = {0:.1f}$'.format( current_z ) 
  ax.text(0.8, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size, color=text_color)


  ax.set_xlabel(r'$\log_{10} \, k \,\,[\,\mathrm{s\,km^{-1} }\,]$', fontsize=label_size )
  ax.set_ylabel(r'$\log_{10} \, k \,\,[\,\mathrm{s\,km^{-1} }\,]$', fontsize=label_size )
  

  ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in')
  ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')
  [sp.set_linewidth(border_width) for sp in ax.spines.values()]

  [sp.set_linewidth(border_width) for sp in cb.ax.spines.values()]



fileName = output_dir + f'cov_matrix_Ns{n_in_sample}'
fileName += '.pdf'
fig.savefig( fileName,  pad_inches=0.1, bbox_inches='tight', dpi=fig_dpi)
print('Saved Image: ', fileName)
