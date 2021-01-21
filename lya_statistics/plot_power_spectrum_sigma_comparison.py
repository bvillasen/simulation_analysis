import sys, os, time
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import pickle
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from flux_power_spectrum import get_skewer_flux_power_spectrum
from tools import *
from stats_functions import compute_distribution, get_highest_probability_interval

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'

# dataDir = '/gpfs/alpine/proj-shared/ast149/'
# dataDir = '/data/groups/comp-astro/bruno/'
dataDir = '/raid/bruno/data/'

n_points = 2048
nx = n_points
ny = n_points
nz = n_points
ncells = nx * ny * nz

uvb = 'pchw18'
# uvb = 'hm12'

# cosmo_name = 'cosmo_3'
cosmo_name = ''


input_dir   = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/flux_power_spectrum_grid_{1}/'.format(n_points, uvb, cosmo_name )
figures_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/figures/flux_ps_diagnostics/'.format( n_points )



snaps = [ 83, 90,  96, 102,  119, 124, 130, 136, 143, 151, 159, 169, ]
snaps_boss = [  96,  102, 106, 110,  114, 119, 124, 130, 136, 143, 151, 159 ]
snapshots = list( set( snaps_boss ).union(set(snaps)))
snapshots.sort()
# print(snapshots)

n_snapshot = snapshots[-2]
snapshot_dir = figures_dir + f'snap_{n_snapshot:03}/'
create_directory( snapshot_dir )

axis = 'x'
  

file_name = input_dir + f'ps_grid_stats_{axis}_{n_snapshot}.pkl'
print ( f'Loading File: {file_name }' )
ps_stats = pickle.load( open( file_name, 'rb' ) ) 



# for k_index in ps_stats:
#   data = ps_stats[k_index]
#   k_val = data['k_val']
#   delta_mean = data['delta_mean']
#   print( k_val, delta_mean )




fig_width = 8
fig_dpi = 300

label_size = 16
figure_text_size = 14

k_index = 3
for k_index in range(1, 18):

  data = ps_stats[k_index]

  nrows = 1
  ncols = 2
  fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols*fig_width,6*nrows))
  plt.subplots_adjust( hspace = 0.2, wspace=0.3)



  k_val = data['k_val']
  delta_mean = data['delta_mean']
  delta_sigma = data['delta_sigma']
  bin_centers = data['bin_centers']
  distribution = data['distribution']
  sigma_l = data['sigma_left']
  sigma_r = data['sigma_right']
  sigma_c = data['sigma_max']
  mean_samples = data['mean_samples']
  n_independent = data['n_independent']
  n_groups = len( mean_samples )

  n_bins_for_distribution = int(np.sqrt(n_groups)) 
  mean_distribution, mean_bin_centers = compute_distribution( mean_samples, n_bins_for_distribution, log=True )

  mean_samples_mean = mean_samples.mean()
  mean_samples_sigma = mean_samples.std()
  mean_samples_sigma_l = mean_samples_mean - mean_samples_sigma
  mean_samples_sigma_r = mean_samples_mean + mean_samples_sigma

  ids  = np.where( distribution > 0.001 )[0]
  left  = ids.min()
  right = ids.max()


  ax = ax_l[0]

  ax.plot( bin_centers, distribution )
  ax.axvline( x=delta_mean, c='C1'  )
  ax.axvline( x=sigma_l, linestyle='--', c='C1'  )
  ax.axvline( x=sigma_r, linestyle='--', c='C1'  )

  text =  r'$k={0:.2e}$ '.format(k_val) + r'$\,\,  \mathrm{s}\,\mathrm{km}^{-1}$' 
  ax.text(0.2, 0.95, text, horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size ) 

  x_text = 0.03
  y_text = 0.85
  delta_text = 0.05

  text = r'$avr: {0:.1e}$ '.format( delta_mean, ) 
  ax.text(x_text, y_text, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size ) 

  text = r'$\sigma_p:\,\, {0:.1e}$'.format(sigma_r - delta_mean)
  ax.text(x_text, y_text - delta_text, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size ) 

  text = r'$\sigma_m: {0:.1e}$'.format( delta_mean - sigma_l)
  ax.text(x_text, y_text - 2*delta_text, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size ) 

  text = r'$\sigma: {0:.1e}$'.format( delta_sigma)
  ax.text(x_text, y_text - 3*delta_text, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size ) 

  text = r'$\sigma / \sqrt{N_{ind}}:$' +' ${0:.1e}$'.format( delta_sigma / np.sqrt( n_independent ))
  ax.text(x_text, y_text - 4*delta_text, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size ) 


  ax.set_xscale('log')

  ax.set_ylim( 0, distribution.max()* 1.2 )
  ax.set_xlim( bin_centers[left], bin_centers[right])

  ax.set_xlabel( r' $\Delta_F^2$', fontsize=label_size )
  ax.set_ylabel( r' $P\,[\Delta_F^2]$', fontsize=label_size )


  ax = ax_l[1]

  ax.plot( mean_bin_centers, mean_distribution )
  ax.axvline( x=mean_samples_mean, c='C2'  )
  ax.axvline( x=mean_samples_sigma_l, linestyle='--', c='C2'  )
  ax.axvline( x=mean_samples_sigma_r, linestyle='--', c='C2'  )

  text =  r'$k={0:.2e}$ '.format(k_val) + r'$\,\,  \mathrm{s}\,\mathrm{km}^{-1}$' 
  ax.text(0.2, 0.95, text, horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size ) 



  text = r'$avr: {0:.1e}$ '.format( mean_samples_mean ) 
  ax.text(x_text, y_text, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size ) 

  text = r'$\sigma:\,\, {0:.1e}$'.format(mean_samples_sigma)
  ax.text(x_text, y_text - delta_text, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size ) 

  text = r'$N_{ind}:$' +' ${0:}$'.format( n_independent)
  ax.text(x_text, y_text - 2*delta_text, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size ) 

  ax.set_xscale('log')

  ax.set_ylim( 0, mean_distribution.max()* 1.2 )

  ax.set_xlabel( r' $\Delta_F^2$', fontsize=label_size )
  ax.set_ylabel( r' $P\,[\Delta_F^2]$', fontsize=label_size )




  fileName = snapshot_dir + f'delta_PS_distribution_k{k_index:03}.png'
  fig.savefig( fileName,  pad_inches=0.1, bbox_inches='tight', dpi=fig_dpi)
  print('Saved Image: ', fileName)


