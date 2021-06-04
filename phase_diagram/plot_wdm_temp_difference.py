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
  data_sim = { 'z':[], 'T0':[], 'gamma':[]}
  for n_file in range(n_files):
    file_name = input_dir + f'{n_file}_analysis.h5'
    infile = h5.File( file_name, 'r' )
    current_z = infile.attrs['current_z'][0]
    infile.close()
    file_name = input_dir + f'fit_mcmc/fit_{n_file}.pkl'
    data = Load_Pickle_Directory( file_name )
    T0 = 10**data['T0']['mean'] 
    gamma =  1 + data['gamma']['mean']
    data_sim['z'].append( current_z )
    data_sim['T0'].append( T0 )
    data_sim['gamma'].append( gamma )
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



from scipy import interpolate as interp 

ax = ax_l[0]

for sim_id in range(n_sims):
  data_sim = data_all[sim_id]
  z_vals = data_sim['z']
  T0 = data_sim['T0'] / 1e4 
  label = data_sim['label']
  n_samples_intgerp = 10000
  z_interp = np.linspace( z_vals[0], z_vals[-1], n_samples_intgerp )  
  f = interp.interp1d( z_vals, T0, kind='cubic' )
  ax.plot( z_interp, f(z_interp), lw=2, c=colors[sim_id], label=label )

ax.tick_params(axis='both', which='major', direction='in', labelsize=label_size )
ax.tick_params(axis='both', which='minor', direction='in' )
ax.set_ylabel( r'$T_0 \,\,\, [\,10^4\, \mathrm{K}\,]$', fontsize=font_size  )
ax.set_xlabel( r'$z$', fontsize=font_size )
leg = ax.legend(loc=1, frameon=False, fontsize=font_size, prop=prop)
ax.set_xlim( 2, 8 )
ax.set_ylim( 0.6, 1.8)
# ax.set_yscale('log')


ax = ax_l[1]

for sim_id in range(n_sims):
  data_sim = data_all[sim_id]
  z_vals = data_sim['z']
  gamma = data_sim['gamma']
  label = data_sim['label']
  # z_interp = np.linspace( z_vals[0], z_vals[-1], n_samples_intgerp )  
  # f = interp.interp1d( z_vals, gamma, kind='cubic' )
  ax.plot( z_vals, gamma, lw=2, c=colors[sim_id], label=label )



ax.tick_params(axis='both', which='major', direction='in', labelsize=label_size )
ax.tick_params(axis='both', which='minor', direction='in' )
ax.set_ylabel( r'$\gamma$', fontsize=font_size  )
ax.set_xlabel( r'$z$', fontsize=font_size )
leg = ax.legend(loc=1, frameon=False, fontsize=font_size, prop=prop)
ax.set_xlim( 2, 8 )
# ax.set_ylim( 0., 7)
# ax.set_yscale('log')



figure_name = output_dir + f'fig_wdm_temp_difference.png'
fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
print( f'Saved Figure: {figure_name}' )