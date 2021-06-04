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
  data_sim = { 'z':[], 'tau_HI':[], 'tau_HeII':[]}
  for n_file in range(n_files):
    file_name = input_dir + f'{n_file}_analysis.h5'
    infile = h5.File( file_name, 'r' )
    current_z = infile.attrs['current_z'][0]
    lya_data = infile['lya_statistics']
    F_mean_HI   = lya_data.attrs['Flux_mean_HI']
    F_mean_HeII = lya_data.attrs['Flux_mean_HeII']
    F_mean_HI[F_mean_HI<1e-20] = 1e-20
    F_mean_HeII[F_mean_HeII<1e-20] = 1e-20
    tau_HI   = -np.log( F_mean_HI )
    tau_HeII = -np.log( F_mean_HeII )
    data_sim['z'].append( current_z )
    data_sim['tau_HI'].append( tau_HI )
    data_sim['tau_HeII'].append( tau_HeII )
  data_sim = { key:np.array(data_sim[key]) for key in data_sim }
  data_sim['label'] = labels[sim_id]
  data_all[sim_id] = data_sim



nrows, ncols = 1, 2
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))
plt.subplots_adjust( hspace = 0.1, wspace=0.1)

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'
prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=12)


font_size = 18
label_size = 16
alpha = 0.5

colors = [ 'C3', 'C0', 'C1', 'C2']





ax = ax_l[0]

for sim_id in range(n_sims):
  data_sim = data_all[sim_id]
  z_vals = data_sim['z']
  tau = data_sim['tau_HI']
  label = data_sim['label']
  ax.plot( z_vals, tau, lw=2, c=colors[sim_id], label=label )

ax.tick_params(axis='both', which='major', direction='in', labelsize=label_size )
ax.tick_params(axis='both', which='minor', direction='in' )
ax.set_ylabel( r'$\tau_{eff} \,\, \mathrm{HI}$', fontsize=font_size  )
ax.set_xlabel( r'$z$', fontsize=font_size )
leg = ax.legend(loc=2, frameon=False, fontsize=font_size, prop=prop)
ax.set_xlim( 2, 6 )
ax.set_ylim( 0.1, 9)
ax.set_yscale('log')


ax = ax_l[1]

for sim_id in range(n_sims):
  data_sim = data_all[sim_id]
  z_vals = data_sim['z']
  tau = data_sim['tau_HeII']
  label = data_sim['label']
  ax.plot( z_vals, tau, lw=2, c=colors[sim_id], label=label )



ax.tick_params(axis='both', which='major', direction='in', labelsize=label_size )
ax.tick_params(axis='both', which='minor', direction='in' )
ax.set_ylabel( r'$\tau_{eff} \,\, \mathrm{HeII}$', fontsize=font_size  )
ax.set_xlabel( r'$z$', fontsize=font_size )
leg = ax.legend(loc=2, frameon=False, fontsize=font_size, prop=prop)
ax.set_xlim( 2, 3.2 )
ax.set_ylim( 0., 7)
# ax.set_yscale('log')



figure_name = output_dir + f'fig_wdm_tau_difference.png'
fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
print( f'Saved Figure: {figure_name}' )