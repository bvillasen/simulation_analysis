import os, sys
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.font_manager
import pylab
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from tools import *
from data_optical_depth_HeII import data_tau_HeII


uvb = 'pchw18'
# dataDir = '/home/bruno/Desktop/ssd_0/data/'
dataDir = '/raid/bruno/data/'
# dataDir = '/data/groups/comp-astro/bruno/'
simulation_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/'
input_dir = simulation_dir + 'skewers_{0}_HeII/los_F/'.format(uvb)
output_dir = simulation_dir + 'figures/'.format(uvb)
create_directory( output_dir )




chem_type = 'HeII'

data = {}
z, tau_vals, tau_sigma_vals = [], [], []

snapshots = range( 100, 170 )
for index, n_snap in enumerate( snapshots ):

  file_name = input_dir + f'los_transmitted_flux_{n_snap}_{chem_type}.h5'
  file = h5.File( file_name, 'r')
  current_z = file.attrs['current_z']
  F_los = file['los_F'][...]
  file.close()

  F_vals = F_los.mean( axis=1 )
  F_mean  = F_vals.mean() 
  F_sigma = F_vals.std()

  tau = -np.log( F_mean )
  tau_sigma = 1/F_mean * F_sigma
  tau_sigma = min( tau_sigma, .30*tau)
  z.append( current_z )
  tau_vals.append(tau)
  tau_sigma_vals.append( tau_sigma )

z = np.array( z )
tau = np.array( tau_vals )
tau_sigma = np.array( tau_sigma_vals )




matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'

system = 'Shamrock'
if system == 'Lux':      prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/brvillas/fonts', "Helvetica.ttf"), size=12)
if system == 'Shamrock': prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=12)

nrows = 1
ncols = 1

font_size = 18
label_size = 16
alpha = 0.8

c_pchw18 = pylab.cm.viridis(.7)
c_hm12 = pylab.cm.cool(.3)

c_boss = pylab.cm.viridis(.3)
c_walther = pylab.cm.viridis(.3)
c_viel = 'C1'
c_boera = pylab.cm.Purples(.7)

text_color  = 'black'
color_line = c_pchw18
color_data = c_boss

fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,5*nrows))


tau_p, tau_m = tau+tau_sigma, tau-tau_sigma

ax.plot( z, tau , label='CHIPS.P19', c=c_pchw18 )
ax.fill_between( z, tau_p, tau_m, facecolor=c_pchw18, alpha=0.3,  )


data_z = data_tau_HeII['z']
data_tau = data_tau_HeII['tau']


ax.scatter( data_z, data_tau, color=color_data, s=4, label='Worseck et al. 2019' )

# 
# data_set = comparable_data[obs_name]
# data_z = data_set['z']
# data_mean = data_set['mean'] 
# data_error = data_set['sigma'] 
# ax.errorbar( data_z, data_mean, yerr=data_error, fmt='none',  alpha=0.8, ecolor= color_data, zorder=2)
# ax.scatter( data_z, data_mean, label='Data for MCMC Fit', alpha=0.8, color= color_data, zorder=2) 

ax.tick_params(axis='both', which='major', direction='in', labelsize=label_size )
ax.tick_params(axis='both', which='minor', direction='in' )
ax.set_ylabel( r'$\tau_{eff} \,\, \mathrm{HeII}$', fontsize=font_size  )
ax.set_xlabel( r'$z$', fontsize=font_size )
leg = ax.legend(loc=1, frameon=False, fontsize=font_size, prop=prop)
ax.set_xlim( 2, 3.4 )
ax.set_ylim( 0, 8 )



figure_name = output_dir + f'fig_tau_HeII_CHIPS_fixed.png'
fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
print( f'Saved Figure: {figure_name}' )

