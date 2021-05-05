import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import palettable
import pylab

cosmo_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
sys.path.append( cosmo_dir + 'lya_statistics/data' )
from tools import *
from data_optical_depth import *


data_dir = '/raid/bruno/data/'
input_dir_0  = data_dir + f'cosmo_sims/rescaled_P19/1024_50Mpc/analysis_files/'
input_dir_1  = data_dir + f'cosmo_sims/rescaled_P19/2048_100Mpc/analysis_files/'
input_dir_2  = data_dir + f'cosmo_sims/rescaled_P19/2048_200Mpc/analysis_files/'
output_dir = data_dir + f'cosmo_sims/rescaled_P19/figures/'

input_dir_list = [ input_dir_0, input_dir_1, input_dir_2 ]

line_widths = [ 2, 2, 2 ]
line_styles = [ 'solid', 'dashed', 'dashed' ]

labels = [ r'$50 \,\,\,\, h^{-1}\mathrm{Mpc} \,\,\, 1024^3 $', r'$100 \,h^{-1}\mathrm{Mpc} \,\,\, 2048^3 $', r'$200 \,h^{-1}\mathrm{Mpc} \,\,\, 2048^3 $']



n_files = range( 14, 56)

data_sets = []

for data_id, input_dir in  enumerate(input_dir_list) :

  z_vals, tau_HI_vals, tau_HeII_vals = [], [], []
  for n_file in n_files:
    file_name = input_dir + f'{n_file}_analysis.h5'
    print( f'Loading File: {file_name} ' )
    file = h5.File( file_name, 'r' )
    current_z = file.attrs['current_z'][0]
    F_mean_HI   = file['lya_statistics'].attrs['Flux_mean_HI'][0]
    F_mean_HeII = file['lya_statistics'].attrs['Flux_mean_HeII'][0]
    tau_HI   = -np.log(F_mean_HI)
    tau_HeII = -np.log(F_mean_HeII)
    z_vals.append( current_z )
    tau_HI_vals.append( tau_HI )
    tau_HeII_vals.append( tau_HeII )
  data_set = {}
  data_set['z'] = z_vals
  data_set['tau_HI']   = np.array(tau_HI_vals)
  data_set['tau_HeII'] = np.array(tau_HeII_vals)
  data_sets.append( data_set )


n_data_sets = len(data_sets)

system = 'Shamrock'

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'

if system == 'Lux':      prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/brvillas/fonts', "Helvetica.ttf"), size=12)
if system == 'Shamrock': prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=12)

nrows = 1
ncols = 2

font_size = 18
label_size = 16
alpha = 0.5

c_pchw18 = pylab.cm.viridis(.7)
c_hm12 = pylab.cm.cool(.3)

c_boss = pylab.cm.viridis(.3)
c_walther = pylab.cm.viridis(.3)
c_viel = 'C1'
c_boera = pylab.cm.Purples(.7)

text_color  = 'black'
color_line = c_pchw18
color_becker = c_boss
color_bosman = c_viel

blue = 'C0'
orange = 'C1'
green = 'C2'


colors = [ blue, orange, green ]


fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))
plt.subplots_adjust( hspace = 0.1, wspace=0.1)

ax = ax_l[0]

for data_id in range( n_data_sets ):
  data_set = data_sets[data_id]
  z = data_set['z']
  tau = data_set['tau_HI']
  linewidth = line_widths[data_id]
  linestyle = line_styles[data_id]
  color_line = colors[data_id]
  label  = labels[data_id]
  ax.plot( z, tau,  linewidth=linewidth, color=color_line, zorder=1, label=label, linestyle=linestyle   )
  


ax.tick_params(axis='both', which='major', direction='in', labelsize=label_size )
ax.tick_params(axis='both', which='minor', direction='in' )
ax.set_ylabel( r'$\tau_{eff} \,\, \mathrm{HI}$', fontsize=font_size  )
ax.set_xlabel( r'$z$', fontsize=font_size )
leg = ax.legend(loc=2, frameon=False, fontsize=font_size, prop=prop)
ax.set_xlim( 2, 6 )
ax.set_ylim( 0.1, 6)
ax.set_yscale('log')


ax = ax_l[1]
for data_id in range( n_data_sets ):
  data_set = data_sets[data_id]
  z = data_set['z']
  tau = data_set['tau_HeII']
  linewidth = line_widths[data_id]
  linestyle = line_styles[data_id]
  color_line = colors[data_id]
  label  = labels[data_id]
  ax.plot( z, tau,  linewidth=linewidth, color=color_line, zorder=1, label=label, linestyle=linestyle   )
  


ax.tick_params(axis='both', which='major', direction='in', labelsize=label_size )
ax.tick_params(axis='both', which='minor', direction='in' )
ax.set_ylabel( r'$\tau_{eff} \,\, \mathrm{HeII}$', fontsize=font_size  )
ax.set_xlabel( r'$z$', fontsize=font_size )
leg = ax.legend(loc=2, frameon=False, fontsize=font_size, prop=prop)
ax.set_xlim( 2, 3.2 )
ax.set_ylim( 0., 7)
# ax.set_yscale('log')



figure_name = output_dir + f'fig_tau_res.png'
fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
print( f'Saved Figure: {figure_name}' )
