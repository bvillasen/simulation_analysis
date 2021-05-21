import os, sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
import pylab
import palettable
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from simulation_parameters import *


ps_data_dir = 'lya_statistics/data/'
input_dir = root_dir + 'interpolated_observables/'
output_dir = root_dir + 'figures/interpolated_ps/'
create_directory( output_dir )

in_file_name = input_dir + 'interpolated_observables.pkl'
fields_data = Load_Pickle_Directory( in_file_name )['fields']

params = fields_data['params']
temp_data =fields_data['T0']
z = temp_data['z']


n_to_plot = 4
indices_to_plot = { 0:[1,2,3,4], 1:[0,1,2,3,], 2:[0,1,2,3,], 3:[1,2,3,4] }
vary_params = [ 'scale_He', 'scale_H', 'deltaZ_He', 'deltaZ_H' ]

labels = {'scale_He':' \\beta_{\mathrm{He}}', 'scale_H':' \\beta_{\mathrm{H}}', 'deltaZ_He':'\Delta z _{\mathrm{He}}', 'deltaZ_H':'\Delta z _{\mathrm{H}}' }


fig_width = 6
fig_height = 2
fig_dpi = 500
label_size = 14
figure_text_size = 14
legend_font_size = 10
tick_label_size_major = 10
tick_label_size_minor = 13
tick_size_major = 4
tick_size_minor = 3
tick_width_major = 1.5
tick_width_minor = 1
border_width = 1
text_color  = 'black'
linewidth = 2
alpha_bar = 0.5

colors = palettable.cmocean.sequential.Haline_10_r.mpl_colors
colors_1 = palettable.colorbrewer.sequential.PuBu_9.mpl_colors
purples = palettable.colorbrewer.sequential.Purples_9.mpl_colors
yellows = palettable.colorbrewer.sequential.YlOrRd_9.mpl_colors 



c_0 = colors[-3]
c_1 = colors[4]
c_2 = colors_1[4]
c_3 = purples[-1]
c_4 = yellows[3]
    
colors = [ c_2, c_1, c_0, c_3  ] 

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'
prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=11)
  
nrows, ncols = 4, 1
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_width*ncols,fig_height*nrows))
plt.subplots_adjust( hspace = 0.0, wspace=0.04)

for indx_i, vary_param in enumerate(vary_params):
  
    id = indx_i
    ax = ax_l[indx_i]
    
    temp_param = temp_data[vary_param]
    
    for i in range(n_to_plot):
      temp = temp_param[i]
      label_param = labels[vary_param]
      color = colors[i]
      label = r'${0} \, = \, {1:.1f}$'.format(label_param, p_val)
      ax.plot( k_vals, ps, linewidth=1, label=label, color=colors[i] )
    # 
    # 
    # legend_loc = 3
    # leg = ax.legend(  loc=legend_loc, frameon=False, prop=prop, fontsize=legend_font_size    )
    
    # ax.text(0.85, 0.95, r'$z={0:.1f}$'.format(z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size, color=text_color) 
    
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    # 
    # if indx_j > 0:ax.set_yticklabels([])
    # if indx_i != nrows-1 :ax.set_xticklabels([])
    # 
    ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in' )
    ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')
    
    ax.set_ylabel( r' $T_0$', fontsize=label_size, color= text_color )
    if indx_i == nrows-1: ax.set_xlabel( r'$z$',  fontsize=label_size, color= text_color )
    
    



  


# file_name = output_dir + f'flux_power_spectrum_interp_{z_index}.png'
file_name = output_dir + f'temp_interp.png'
fig.savefig( file_name,  pad_inches=0.1, bbox_inches='tight', dpi=fig_dpi)
print('Saved Image: ', file_name )


