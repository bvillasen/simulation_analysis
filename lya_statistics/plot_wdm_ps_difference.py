import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import palettable
import pylab
import pickle
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from tools import *
from plot_flux_power_spectrum import plot_power_spectrum_grid


data_dir = '/raid/bruno/data/'
input_dir_0  = data_dir + f'cosmo_sims/rescaled_P19/1024_50Mpc/analysis_files/'
input_dir_1  = data_dir + f'cosmo_sims/rescaled_P19/wdm/1024_50Mpc_wdm_m0.5kev/analysis_files/'
input_dir_2  = data_dir + f'cosmo_sims/rescaled_P19/wdm/1024_50Mpc_wdm_m1.0kev/analysis_files/'
input_dir_3  = data_dir + f'cosmo_sims/rescaled_P19/wdm/1024_50Mpc_wdm_m3.0kev/analysis_files/'
output_dir = data_dir + f'cosmo_sims/rescaled_P19/wdm/figures/'
create_directory( output_dir )

input_dir_list = [ input_dir_0, input_dir_1, input_dir_2, input_dir_3 ]
labels = ['CDM', 'WDM m = 0.5 keV', 'WDM m = 1.0 keV', 'WDM m = 3.0 keV' ]

n_sims = len(input_dir_list)
n_files = 56

data_all = {}
for sim_id in range(n_sims):
  input_dir = input_dir_list[sim_id]
  data_sim = {}
  for n_file in range(n_files):
    file_name = input_dir + f'{n_file}_analysis.h5'
    infile = h5.File( file_name, 'r' )
    current_z = infile.attrs['current_z'][0]
    ps_data = infile['lya_statistics']['power_spectrum']
    k_vals = ps_data['k_vals'][...]
    ps_mean = ps_data['p(k)'][...]
    indices = np.where( ps_mean > 0 )
    ps_mean = ps_mean[indices]
    k_vals = k_vals[indices]
    data_sim[n_file] = { 'z':current_z, 'k_vals':k_vals, 'ps_mean':ps_mean  }
  sim_z_vals = np.array([ data_sim[i]['z'] for i in range(n_files) ])
  data_sim['label'] = labels[sim_id]
  data_sim['z_vals'] = sim_z_vals
  data_all[sim_id] = data_sim
  


z_vals = [ 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5  ]


ps_data_dir = os.getcwd() + '/data/'
plot_power_spectrum_grid( ps_data_dir, output_dir, ps_data=data_all, scales='small_highz', system='Shamrock',  )
  



# nrows, ncols = 4, 2
# fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(7.5*ncols,5*nrows), sharey='row', sharex='col')
# plt.subplots_adjust( hspace = 0.02, wspace=0.02)
# 
# import matplotlib
# matplotlib.rcParams['mathtext.fontset'] = 'cm'
# matplotlib.rcParams['mathtext.rm'] = 'serif'
# prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=12)
# 
# 
# colors = [ 'C3', 'C0', 'C1', 'C2']
# 
# fig_dpi = 300
# label_size = 18
# figure_text_size = 18
# legend_font_size = 16
# tick_label_size_major = 15
# tick_label_size_minor = 13
# tick_size_major = 5
# tick_size_minor = 3
# tick_width_major = 1.5
# tick_width_minor = 1
# border_width = 1
# text_color  = 'black'
# linewidth = 2
# alpha_bar = 0.5
# 
# for index, z in enumerate(z_vals):
# 
#   indx_j = index % ncols
#   indx_i = index//ncols
#   ax = ax_l[indx_i][indx_j]
# 
#   z_diff = np.abs( sim_z_vals - z)
#   z_indx = np.where( z_diff == z_diff.min() )[0][0]
# 
#   data_sim = data_all[0][z_indx]
#   k_vals_0 = data_sim['k_vals']
#   ps_mean_0 = data_sim['ps_mean'] 
# 
#   for sim_id in [ 1, 2, 3]:
#     data_sim = data_all[sim_id][z_indx]
#     k_vals = data_sim['k_vals']
#     ps_mean = data_sim['ps_mean'] 
#     ps_diff = ( ps_mean - ps_mean_0 ) / ps_mean_0
#     if (k_vals != k_vals_0).any(): print( 'ERROR: K vals mismatch')
#     label = data_sim['label']
#     ax.plot( k_vals, ps_diff, lw=2, c=colors[sim_id], label=label )
# 
#   ax.  axhline( y=0, ls='--', color='C3')
#   ax.text(0.9, 0.93, r'$z={0:.1f}$'.format(z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size, color=text_color) 
# 
#   if indx_j == 0 and indx_i == 0:
#     legend_loc = 3
#     ax.legend( loc=legend_loc, frameon=False, prop=prop )
# 
#   ax.set_xscale( 'log' )  
#   ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in' )
#   ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')
# 
#   [sp.set_linewidth(border_width) for sp in ax.spines.values()]
#   if indx_j == 0: ax.set_ylabel( r' $\Delta P\,(k) / P\,(k)$', fontsize=label_size, color= text_color )
#   if indx_i == nrows-1: ax.set_xlabel( r'$ k   \,\,\,  [\mathrm{s}\,\mathrm{km}^{-1}] $',  fontsize=label_size, color= text_color )
# 
# figure_name = output_dir + f'fig_wdm_ps_difference.png'
# fig.savefig( figure_name, bbox_inches='tight', dpi=500 )
# print( f'Saved Figure: {figure_name}' )