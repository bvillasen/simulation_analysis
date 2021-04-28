import os, sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from parameters_UVB_rates import param_UVB_Rates
from simulation_grid import Simulation_Grid
from simulation_parameters import *
from plot_flux_power_spectrum import plot_power_spectrum_grid
from plot_T0_tau import plot_T0_and_tau, plot_tau_HeII
from mcmc_data_functions import Get_Comparable_Composite
from grid_ploting_functions import Plot_Grid_TO


ps_data_dir = 'lya_statistics/data/'
output_dir = root_dir + 'figures/'
create_directory( output_dir )

SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
sim_ids = list(SG.sim_ids)
# sim_ids = range(10)

rescaled_ps = True
ps_norm = {'normalization':'Becker', 'type':'tau_eff_global'}
SG.Load_Grid_Analysis_Data( sim_ids=sim_ids, load_fit=True, load_normalized_ps=rescaled_ps, ps_norm=ps_norm )


field = 'T0+tau'
z_min = 2.0
z_max = 5.0 
tau_extras = {'factor_sigma_becker':6.0, 'factor_sigma_keating':4.0}
comparable_data = Get_Comparable_Composite( field, z_min, z_max, tau_extras=tau_extras  )

sim_data_sets = [ ]
# sim_ids_to_plot = [ 1 ]
sim_ids_to_plot = sim_ids
for sim_id in sim_ids_to_plot:
  sim_data  = SG.Grid[sim_id]['analysis']
  sim_param = SG.Grid[sim_id]['parameters']
  beta_He = sim_param['scale_He'] 
  beta_H  = sim_param['scale_H'] 
  deltaZ_He = sim_param['deltaZ_He'] 
  # label = r'$\beta_{\mathrm{He}}:$' + f'{beta_He:.1f}' + ' ' + r'$\beta_{\mathrm{H}}:$' + f'{beta_H:.1f}' + ' ' + r'$\Delta_z:$' + f'{deltaZ_He:.1f}' 
  # label = 'P19 Modified'
  label = ''
  sim_data['plot_label'] = label
  sim_data_sets.append(sim_data)
  
  
Plot_Grid_TO( sim_data_sets, output_dir, plot_splines=True, system=system )


z_val = 2.6

black_background = True

from load_tabulated_data import *

current_z = z_val

rescale_walter_file = ps_data_dir + 'rescale_walther_to_boss.pkl' 
rescaled_walther = True


import palettable
import pylab
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'
if system == 'Lux':      prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/brvillas/fonts', "Helvetica.ttf"), size=12)
if system == 'Shamrock': prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=12)


if system == 'Lux': matplotlib.use('Agg')


fig_height = 5
fig_width = 8
fig_dpi = 300

alpha = 0.6

label_size = 18
figure_text_size = 18
legend_font_size = 16
tick_label_size_major = 15
tick_label_size_minor = 13
tick_size_major = 5
tick_size_minor = 3
tick_width_major = 1.5
tick_width_minor = 1
border_width = 1


dir_boss = ps_data_dir + 'data_power_spectrum_boss/'
data_filename = dir_boss + 'data_table.py'
data_boss = load_data_boss( data_filename )
data_z_boss = data_boss['z_vals']

data_filename = ps_data_dir + 'data_power_spectrum_walther_2019/data_table.txt'
data_walther = load_power_spectrum_table( data_filename )
data_z_w = data_walther['z_vals']

dir_data_boera = ps_data_dir + 'data_power_spectrum_boera_2019/'
data_boera = load_tabulated_data_boera( dir_data_boera )
data_z_b = data_boera['z_vals']

data_dir_viel = ps_data_dir + 'data_power_spectrum_viel_2013/'
data_viel = load_tabulated_data_viel( data_dir_viel)
data_z_v = data_viel['z_vals']

if rescaled_walther:
  print(f" Loading Walther rescale values: {rescale_walter_file}")
  rescale_walter_alphas = Load_Pickle_Directory( rescale_walter_file)



c_pchw18 = pylab.cm.viridis(.7)
c_hm12 = pylab.cm.cool(.3)

c_boss = pylab.cm.viridis(.3)
c_walther = pylab.cm.viridis(.3)
c_viel = 'C1'
c_boera = pylab.cm.Purples(.7)

text_color  = 'black'
color_line = c_pchw18
purples = palettable.colorbrewer.sequential.Purples_9.mpl_colors
purple = purples[-1]

oranges = palettable.colorbrewer.sequential.YlOrBr_9.mpl_colors
orange = oranges[4]


blues = palettable.colorbrewer.sequential.Blues_9.mpl_colors
blue = blues[-2]

greens = palettable.colorbrewer.sequential.BuGn_9.mpl_colors
green = greens[-2]


sim_data_sets.reverse()
  
n_lines = len( sim_data_sets )
colormap = palettable.cmocean.sequential.Matter_20.mpl_colormap
colormap = palettable.cmocean.sequential.Tempo_20.mpl_colormap
colormap = palettable.cmocean.sequential.Dense_20.mpl_colormap
colormap = palettable.scientific.sequential.LaPaz_20_r.mpl_colormap
colormap = palettable.colorbrewer.sequential.Blues_9_r.mpl_colormap
colors = colormap( np.linspace(0,1,n_lines) )


text_color  = 'black'
color_line = c_pchw18
color_data = c_boss

black_background = True
if black_background:
  text_color = 'white'
  color_data_0 = orange
  color_data_1 = 'C3'
    
alpha_bar = 0.4

nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(6,7))
plt.subplots_adjust( hspace = 0.02, wspace=0.02)



for sim_id, sim_data in enumerate( sim_data_sets ):
  
  if sim_id %2 == 0: continue
  if sim_id %4 == 0: continue
  
  color_line = colors[sim_id]
  sim_ps = sim_data_sets[sim_id]['power_spectrum']
  sim_z = sim_ps['z']

  # print( f'N z Indices: {n_z_indices}')
  diff = np.abs( sim_z - current_z )
  diff_min = diff.min()
  if sim_id == n_lines - 1: label = 'Simulation'
  else: label = '' 
  if diff_min < 0.05:
    index = np.where( diff == diff_min )[0][0]
    k_vals = sim_ps['k_vals'][index]
    ps_mean = sim_ps['ps_mean'][index]
    delta_mean = ps_mean * k_vals / np.pi
    ax.plot( k_vals, delta_mean, linewidth=1, color=color_line,  zorder=1, label=label, alpha=alpha  )


ax.text(0.85, 0.95, r'$z={0:.1f}$'.format(current_z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size, color=text_color) 



# Add Boss data
z_diff = np.abs( data_z_boss - current_z )
diff_min = z_diff.min()
if diff_min < 1e-1:
  data_index = np.where( z_diff == diff_min )[0][0]
  data_z_local = data_z_boss[data_index]

  data_k = data_boss[data_index]['k_vals']
  data_delta_power = data_boss[data_index]['delta_power']
  data_delta_power_error = data_boss[data_index]['delta_power_error']
  label_boss = 'eBOSS (2019)'
  d_boss = ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=color_data_0, label=label_boss, zorder=2)


# Add Walther data
z_diff = np.abs( data_z_w - current_z )
diff_min = z_diff.min()
if diff_min < 1e-1:
  data_index = np.where( z_diff == diff_min )[0][0]
  data_z_local = data_z_w[data_index]

  data_k = data_walther[data_index]['k_vals']
  data_delta_power = data_walther[data_index]['delta_power']
  data_delta_power_error = data_walther[data_index]['delta_power_error']
  label_walther ='Walther et al. (2018)' 
  
  if rescaled_walther and data_index in rescale_walter_alphas:

    rescale_z = rescale_walter_alphas[data_index]['z']
    rescale_alpha = rescale_walter_alphas[data_index]['alpha']
    print( f'  Rescaling z={rescale_z:.1f}    alpha={rescale_alpha:.3f} ')
    data_delta_power *= rescale_alpha
    d_walther = ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=color_data_1, label=label_walther, zorder=2)

    
  # # Add Boera data
  # z_diff = np.abs( data_z_b - current_z )
  # diff_min = z_diff.min()
  # if diff_min < 1e-1:
  #   data_index = np.where( z_diff == diff_min )[0][0]
  #   data_z_local = data_z_b[data_index]
  # 
  #   data_k = data_boera[data_index]['k_vals']
  #   data_delta_power = data_boera[data_index]['delta_power']
  #   data_delta_power_error = data_boera[data_index]['delta_power_error']
  #   label_boera ='Boera et al. (2019)'
  #   d_boera = ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_boera, label=label_boera, zorder=2 )
  # 
  # 
  # # Add Viel data
  # z_diff = np.abs( data_z_v - current_z )
  # diff_min = z_diff.min()
  # if diff_min < 1e-1:
  #   data_index = np.where( z_diff == diff_min )[0][0]
  #   data_z_local = data_z_v[data_index]
  # 
  #   data_k = data_viel[data_index]['k_vals']
  #   data_delta_power = data_viel[data_index]['delta_power']
  #   data_delta_power_error = data_viel[data_index]['delta_power_error']
  #   label_viel = 'Viel et al. (2013)'
  #   d_viel = ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_viel, label=label_viel, zorder=2 )


legend_loc = 3


label_bars =  r'1$\sigma$ skewers $P\,(\Delta_F^2)$'


leg = ax.legend(  loc=legend_loc, frameon=False, prop=prop    )
for text in leg.get_texts():
  plt.setp(text, color = text_color)

  
x_min, x_max = 1.5e-3, 2e-1 
y_min, y_max = 1e-3, 8e-2


ax.set_xlim( x_min, x_max )
ax.set_ylim( y_min, y_max )
ax.set_xscale('log')
ax.set_yscale('log')


if black_background: 
  fig.patch.set_facecolor('black') 
  ax.set_facecolor('k')
  [ spine.set_edgecolor(text_color) for spine in list(ax.spines.values()) ]
    
[sp.set_linewidth(border_width) for sp in ax.spines.values()]

ax.tick_params(axis='both', which='major', color=text_color, labelcolor=text_color, labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in' )
ax.tick_params(axis='both', which='minor', color=text_color, labelcolor=text_color, labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')

ax.set_ylabel( r' $\Delta_F^2(k)$', fontsize=label_size, color= text_color )
ax.set_xlabel( r'$ k   \,\,\,  [\mathrm{s}\,\mathrm{km}^{-1}] $',  fontsize=label_size, color= text_color )

fileName = output_dir + f'fig_grid_flux_ps'
fileName += '.png'
fig.savefig( fileName,  pad_inches=0.1, bbox_inches='tight', dpi=fig_dpi, facecolor=fig.get_facecolor() )
print('Saved Image: ', fileName)    



# plot_power_spectrum_grid( ps_data_dir, output_dir, scales='small_walther', sim_data_sets=sim_data_sets, system=system, plot_ps_normalized=rescaled_ps )
# plot_power_spectrum_grid( ps_data_dir, output_dir, scales='large', sim_data_sets=sim_data_sets, system=system )
# # plot_power_spectrum_grid( ps_data_dir, output_dir, scales='middle', sim_data_sets=sim_data_sets, system=system )
# 
# 
# 
# 
# plot_T0_and_tau( comparable_data, output_dir, sim_data_sets=sim_data_sets, system=system )
# 
# # plot_tau_HeII( output_dir, sim_data_sets=sim_data_sets, system=system )
