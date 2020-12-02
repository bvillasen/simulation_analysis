import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
from mpl_toolkits.axes_grid1 import ImageGrid
import palettable

root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from tools import *
from load_data_enzo import load_snapshot_enzo
from load_data import load_cholla_snapshot_file
from phase_diagram_functions import *
from turbo_cmap import *

import matplotlib
matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'


# data_dir = '/home/bruno/Desktop/ssd_0/data/'
data_dir = '/raid/bruno/data/'
enzo_dir = data_dir + 'cosmo_sims/enzo/256_cool_uv_50Mpc_HLLC_grav4/h5_files/'
cholla_dir = data_dir + 'cosmo_sims/256_cool_uv_50Mpc/snapshots/'

# output_dir = '/home/bruno/Desktop/'
output_dir = data_dir + 'cosmo_sims/figures/'
create_directory( output_dir )

n_points = 256
n_cells = n_points**3

data = {}
types = [ 'enzo', 'cholla' ]
types = [ 'enzo' ]
data['enzo'] = load_snapshot_enzo( 33, enzo_dir, hydro=True, temp=True )
data['cholla'] = load_cholla_snapshot_file( 33, cholla_dir, dm=False )

current_z = data['enzo']['current_z']
  

#Get Bin Egdes for the histogram
dens_start, dens_end = -3, 5
temp_start, temp_end = 2, 8
nbins = 800
bins_dens = np.logspace( dens_start, dens_end, nbins, base=10 )
bins_temp = np.logspace( temp_start, temp_end, nbins, base=10 )

#Get the phase diagram
pd_data = {}
# for type in types:
type = 'enzo'
print(" Generating Phase Diagram ,   n_bins:{0}".format(nbins))
centers_dens, centers_temp, phase = get_phase_diagram_bins( data[type]['gas']['density'], data[type]['gas']['temperature'], bins_dens, bins_temp  )
temp_points, dens_points = np.meshgrid( centers_temp, centers_dens )
phase /= n_cells
pd_data['data'] = phase
pd_data['dens_min'], pd_data['dens_max'], pd_data['n_dens'] = 10**centers_dens.min(), 10**centers_dens.max(), len( centers_dens)
pd_data['temp_min'], pd_data['temp_max'], pd_data['n_temp'] = 10**centers_temp.min(), 10**centers_temp.max(), len( centers_temp)  
values_to_fit = get_density_tyemperature_values_to_fit( pd_data, delta_min=-0.3, delta_max=1.2, n_samples_line=50, fraction_enclosed=0.70 )
# fit_values = fit_thermal_parameters_mcmc( None, values_to_fit, None, save_file=False )
temp_points = temp_points.flatten()
dens_points = dens_points.flatten()
phase = phase.flatten() 
indices = np.where(phase > 0 )
phase = phase[indices]
dens_points = dens_points[indices]
temp_points = temp_points[indices]
pd_data[type] = { 'dens_points':dens_points, 'temp_points':temp_points, 'phase':phase }
  
  
# 
# 
# fig_width = 8
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
# 
# text_color = 'black'
# 
# n_data = 2
# n_rows = 1
# n_cols = n_data
# 
# title_all = [ 'Enzo', 'Cholla' ]
# 
# # Set up figure and image grid
# fig = plt.figure(0, figsize=(fig_width*n_cols,10*n_rows),  )
# grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
#                  nrows_ncols=(n_rows,n_cols),
#                  axes_pad=0.2,
#                  share_all=True,
#                  cbar_location="right",
#                  cbar_mode="single",
#                  cbar_size="5%",
#                  cbar_pad=0.1,
#                  )
# 
# colormap =  'turbo'
# alpha = 0.6
# 
# x_min, x_max = -1.5, 5.1
# y_min, y_max =  1.8, 8.5
# 
# for i, type in enumerate( types ):
# 
#   ax = grid[i]
# 
#   dens_points = pd_data[type]['dens_points']
#   temp_points = pd_data[type]['temp_points']
#   phase = np.log10(pd_data[type]['phase'])    
#   v_min, v_max = phase.min(), phase.max()
# 
#   im = ax.scatter( dens_points, temp_points, c=phase, s=0.1, vmin=v_min, vmax=v_max, alpha=alpha, cmap=colormap  )
# 
# 
#   cb = ax.cax.colorbar(im,   )
#   cb.ax.tick_params(labelsize=tick_label_size_major, size=tick_size_major, color=text_color, width=tick_width_major, length=tick_size_major, labelcolor=text_color, direction='in' )
#   ax.cax.toggle_label(True)
#   [sp.set_linewidth(border_width) for sp in cb.ax.spines.values()]
# 
# 
# 
#   font = {'fontname': 'Helvetica',
#       'color':  text_color,
#       'weight': 'normal',
#       'size': label_size,
#       'ha':'center'
#       }
#   cb.set_label_text( r'$\log_{10}  \,\, P\,(\Delta, T\,) $', fontdict=font )
#   ax.set_ylabel(r'$\log_{10} \, T \,\,[\,\mathrm{K}\,]$', fontsize=label_size , color=text_color)
#   ax.set_xlabel(r'$\log_{10} \, \Delta$ ', fontsize=label_size , color=text_color )
# 
# 
#   text  = r'$z = {0:.2f}$'.format( current_z ) 
#   if i == 0: ax.text(0.05, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size, color=text_color)
# 
#   text = title_all[i]
#   ax.text(0.95, 0.95, text, horizontalalignment='right',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size, color=text_color)
# 
# 
#   ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in')
#   ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')
#   [sp.set_linewidth(border_width) for sp in ax.spines.values()]
#   ax.set_xlim( x_min, x_max )
#   ax.set_ylim( y_min, y_max )
# 
# out_fileName = output_dir + 'phase_diagram_256.png'
# fig.savefig( out_fileName,  pad_inches=0.1,  facecolor=fig.get_facecolor(),  bbox_inches='tight', dpi=300)
# print(( 'Saved Image: ' + out_fileName ))
# 
# 
