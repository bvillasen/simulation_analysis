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
ps_data = Load_Pickle_Directory( in_file_name )['power_spectrum']

z_vals = ps_data['z_vals']
params = ps_data['params']

k_max = 2e-1
n_to_plot = 4
indices_to_plot = { 0:[1,2,3,4], 1:[0,1,2,3,], 2:[0,1,2,3,], 3:[1,2,3,4] }
vary_params = [ 'scale_He', 'scale_H', 'deltaZ_He', 'deltaZ_H' ]

z_vals = np.array([ 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.6, 4.4,  5.0,  5.4,   ])

# for z in z_vals:
z = 3.8
z_index = np.where( z_vals == z )[0][0]


ps_plot = {}
for param_indx, vary_param in enumerate(vary_params):
  
  ps_plot[param_indx] = { 'param_name': vary_param }
  ps_param = ps_data[vary_param]
  for i,indx in enumerate(indices_to_plot[param_indx]):
    ps_plot[param_indx][i] = { 'k_vals':ps_param[indx]['ps'][z_index]['k_vals'], 'ps':ps_param[indx]['ps'][z_index]['mean'], 'p_val':ps_param[indx]['params'][param_indx] }
    

labels = {'scale_He':' \\beta_{\mathrm{He}}', 'scale_H':' \\beta_{\mathrm{H}}', 'deltaZ_He':'\Delta z _{\mathrm{He}}', 'deltaZ_H':'\Delta z _{\mathrm{H}}' }


fig_width = 4
fig_width = 5
fig_dpi = 300
label_size = 14
figure_text_size = 14
legend_font_size = 10
tick_label_size_major = 12
tick_label_size_minor = 13
tick_size_major = 5
tick_size_minor = 3
tick_width_major = 1.5
tick_width_minor = 1
border_width = 1
text_color  = 'black'
linewidth = 2
alpha_bar = 0.5


# Colors
bright_green = pylab.cm.viridis(.7)
light_blue = pylab.cm.cool(.3)
dark_blue = pylab.cm.viridis(.3) 
purple = pylab.cm.Purples(.7)
blue = 'C0'
orange = 'C1'
green = 'C2'
red = 'C3'
purple_2 = 'C4'

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
  
nrows, ncols = 2, 2
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_width*ncols,5*nrows))
plt.subplots_adjust( hspace = 0.04, wspace=0.04)

for indx_i in range(ncols):
  for indx_j in range(ncols):
    
    id = indx_i*ncols + indx_j
    ax = ax_l[indx_i][indx_j]
    
    ps_param = ps_plot[id]
    param_name = ps_param['param_name']
    for i in range(n_to_plot):
      ps_d = ps_param[i]
      k_vals = ps_d['k_vals']
      indices = k_vals < k_max
      k_vals = k_vals[indices]
      ps = ps_d['ps'][indices]
      p_val = ps_d['p_val']
      label_param = labels[param_name]
      color = colors[i]
      label = r'${0} \, = \, {1:.1f}$'.format(label_param, p_val)
      ax.plot( k_vals, ps, linewidth=1, label=label, color=colors[i] )
    

    legend_loc = 3
    leg = ax.legend(  loc=legend_loc, frameon=False, prop=prop, fontsize=legend_font_size    )
    
    ax.text(0.85, 0.95, r'$z={0:.1f}$'.format(z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size, color=text_color) 
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    if indx_j > 0:ax.set_yticklabels([])
    if indx_i != nrows-1 :ax.set_xticklabels([])

    ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in' )
    ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')

    if indx_j == 0: ax.set_ylabel( r' $\Delta_F^2(k)$', fontsize=label_size, color= text_color )
    if indx_i == nrows-1: ax.set_xlabel( r'$ k   \,\,\,  [\mathrm{s}\,\mathrm{km}^{-1}] $',  fontsize=label_size, color= text_color )



  




# file_name = output_dir + f'flux_power_spectrum_interp_{z_index}.png'
file_name = output_dir + f'flux_power_spectrum_interp.png'
fig.savefig( file_name,  pad_inches=0.1, bbox_inches='tight', dpi=fig_dpi)
print('Saved Image: ', file_name )



