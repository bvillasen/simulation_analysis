import os, sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from parameters_UVB_rates import param_UVB_Rates
from simulation_grid import Simulation_Grid
from simulation_parameters import *
from mcmc_functions import *
from mcmc_data_functions import *
from data_thermal_history import *
from mcmc_plotting_functions import *

# data_sets = [ 'Boss', 'Walther', 'Boera', 'Viel' ]
# data_sets = [ 'Boss' ]
data_sets = [ 'Walther' ]
# data_sets = [ 'Boera' ]
# data_sets = [ 'Boss', 'Walther' ]
# data_sets = [ 'Walther', 'Boera' ]
# data_sets = [ 'Walther', 'Viel' ]


name = ''
for data_set in data_sets:
  name += data_set + '_'
name = name[:-1] 

field = 'P(k)'
ps_data_dir = 'lya_statistics/data/'
output_dir = root_dir + f'fit_results_{field}_{name}/'
create_directory( output_dir )

load_mcmc_stats = False
# load_mcmc_stats = True


SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Load_Grid_Analysis_Data()
ps_range = SG.Get_Power_Spectrum_Range( kmax=0.01 )
sim_ids = SG.sim_ids

z_min = 2.0
z_max = 5.0 
comparable_data = Get_Comparable_Power_Spectrum(  ps_data_dir, z_min, z_max, data_sets, ps_range )
comparable_grid = Get_Comparable_Power_Spectrum_from_Grid( comparable_data['separate'], SG )


stats_file = output_dir + 'fit_mcmc.pkl'

fields = [ 'P(k)' ]
z_vals = [ 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.6, 5.0,  ]
data_grid = Get_Data_Grid_Power_spectrum( z_vals, SG )

params = SG.parameters

if load_mcmc_stats:
  print( f'Loading File: {stats_file}')
  stats = pickle.load( open( stats_file, 'rb' ) )
  param_stats = {}
  for p_id in params.keys():
    p_name = params[p_id]['name']
    p_stats = stats[p_name]
    params[p_id]['mean'] = p_stats['mean']
    params[p_id]['sigma'] = p_stats['standard deviation']

else:
  nIter = 200000
  # nIter = 100000
  nBurn = nIter / 5
  nThin = 1
  model, params_mcmc = mcmc_model_3D( comparable_data, comparable_grid, field, 'mean', SG )
  MDL = pymc.MCMC( model )  
  MDL.sample( iter=nIter, burn=nBurn, thin=nThin )
  stats = MDL.stats()
  param_stats = {}
  for p_id in params.keys():
    p_name = params[p_id]['name']
    p_stats = stats[p_name]
    params[p_id]['mean'] = p_stats['mean']
    params[p_id]['sigma'] = p_stats['standard deviation']
  Plot_MCMC_Stats( stats, MDL, params_mcmc,  stats_file, output_dir, plot_corner=False )

  samples = {} 
  for p_id in params_mcmc.keys():
    param = params_mcmc[p_id]
    samples[p_id] = {}
    samples[p_id]['name'] = param['name']
    samples[p_id]['trace'] = param['sampler'].trace() 
  
  labels = { 'scale_He':r'$\beta_{\mathrm{He}}$', 'scale_H':r'$\beta_{\mathrm{H}}$', 'deltaZ_He':r'$\Delta z_{\mathrm{He}}$'    }
  Plot_Corner( samples, labels, output_dir  )




for p_id in params.keys():
  p_name = params[p_id]['name']
  p_stats = stats[p_name]
  params[p_id]['mean'] = p_stats['mean']
  params[p_id]['sigma'] = p_stats['standard deviation']


# Obtain distribution of the power spectrum
n_samples = 1000
ps_samples = Sample_Power_Spectrum( n_samples, params, data_grid, SG  )
# chi2_vals = Get_Chi2( observables, params, comparable_grid, comparable_data, SG )
# Plot_Observables( observables_samples, comparable_data, params, SG, 'sampling', output_dir, chi2=chi2_vals)
# Plot_Observables( observables_samples, comparable_data, params, SG, 'grid', output_dir, chi2=None)
# 


scales = 'small'


import matplotlib
import pylab

matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'


if system == 'Lux': matplotlib.use('Agg')

fig_height = 5
fig_width = 8
fig_dpi = 300


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

if system == 'Lux':      prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/brvillas/fonts', "Helvetica.ttf"), size=12)
if system == 'Shamrock': prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=12)

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

z_vals_small_scale  = [ 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 4.2, 4.6, 5.0, 5.4 ]
z_vals_large_scale  = [ 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.6 ]
z_vals_middle_scale = [ 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 4.2, 4.6 ]
z_high = [ 5.0, 5.4 ]



n_samples = len( ps_samples )
z_samples = np.array([ ps_samples[i]['z'] for i in range(n_samples) ])

if scales == 'large': z_vals = z_vals_large_scale
elif scales == 'small': z_vals = z_vals_small_scale
elif scales == 'middle': z_vals = z_vals_middle_scale
else: 
  print( "ERROR: Scales = large,  small of middle ")
  # return
  
  
nrows = 3
ncols = 4


if scales == 'middle':flags = np.zeros( (nrows, ncols ))

# if scales == 'middle': nrows = 2


fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(2*fig_width,fig_height*nrows))
plt.subplots_adjust( hspace = 0.02, wspace=0.02)


c_pchw18 = pylab.cm.viridis(.7)
c_hm12 = pylab.cm.cool(.3)

c_boss = pylab.cm.viridis(.3)
c_walther = pylab.cm.viridis(.3)
c_viel = 'C1'
c_boera = pylab.cm.Purples(.7)

text_color  = 'black'
color_line = c_pchw18

if scales == 'middle':
  c_walther = 'C3'

alpha_bar = 0.4

for index, current_z in enumerate( z_vals ):



  indx_j = index % ncols
  indx_i = index//ncols

  if nrows > 1: ax = ax_l[indx_i][indx_j]
  else: ax = ax_l[indx_j]
  
  
  if scales == 'middle': flags[indx_i,  indx_j] = 1

  
  diff = np.abs( z_samples - current_z )
  diff_min = diff.min()
  if diff_min < 0.05:
    index = np.where( diff == diff_min )[0][0]
    k_vals = ps_samples[index]['k_vals']
    delta_mean = ps_samples[index]['mean']
    delta_sigma = ps_samples[index]['sigma'] 
    delta_p = delta_mean + delta_sigma
    delta_m = delta_mean - delta_sigma 
    ax.plot( k_vals, delta_mean, linewidth=3, color=color_line, zorder=1   )
    ax.fill_between( k_vals, delta_p, delta_m, facecolor=color_line, alpha=alpha_bar, zorder=1   )
  
  #       # ax.plot( k, delta, c=color_line, linewidth=3, label=sim_data['plot_label']  )
  # 

  ax.text(0.85, 0.95, r'$z={0:.1f}$'.format(current_z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size, color=text_color) 

  
  if scales == 'large' or scales == 'middle':
  
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
      d_boss = ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_boss, label=label_boss, zorder=2)
  
  if scales == 'small' or scales == 'middle':
  
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
      d_walther = ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_walther, label=label_walther, zorder=2)
  
  
    # Add Boera data
    z_diff = np.abs( data_z_b - current_z )
    diff_min = z_diff.min()
    if diff_min < 1e-1:
      data_index = np.where( z_diff == diff_min )[0][0]
      data_z_local = data_z_b[data_index]
  
      data_k = data_boera[data_index]['k_vals']
      data_delta_power = data_boera[data_index]['delta_power']
      data_delta_power_error = data_boera[data_index]['delta_power_error']
      label_boera ='Boera et al. (2019)'
      d_boera = ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_boera, label=label_boera, zorder=2 )
  
  
    # Add Viel data
    z_diff = np.abs( data_z_v - current_z )
    diff_min = z_diff.min()
    if diff_min < 1e-1:
      data_index = np.where( z_diff == diff_min )[0][0]
      data_z_local = data_z_v[data_index]
  
      data_k = data_viel[data_index]['k_vals']
      data_delta_power = data_viel[data_index]['delta_power']
      data_delta_power_error = data_viel[data_index]['delta_power_error']
      label_viel = 'Viel et al. (2013)'
      d_viel = ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_viel, label=label_viel, zorder=2 )
  

  legend_loc = 3
  if indx_i == nrows-1 and nrows!=2: legend_loc = 2

  if scales == 'large': legend_loc = 2
  label_bars =  r'1$\sigma$ skewers $P\,(\Delta_F^2)$'

  add_legend = False
  if indx_j == 0: add_legend = True
  
  if scales == 'middle' and indx_i == nrows-1 and indx_j == ncols-1: add_legend = True
    
    
  if add_legend:
    # leg = ax.legend( loc=legend_loc, frameon=False, fontsize=12)
    leg = ax.legend(  loc=legend_loc, frameon=False, prop=prop    )
    
    for text in leg.get_texts():
        plt.setp(text, color = text_color)
        



  x_min, x_max = 4e-3, 2.5e-1
  if indx_i == 0: y_min, y_max = 1e-3, 9e-2
  if indx_i == 1: y_min, y_max = 5e-3, 2e-1
  if indx_i == 2: y_min, y_max = 5e-2, 3

  if scales == 'large':
    x_min, x_max = 2e-3, 2.3e-2
    if indx_i == 0: y_min, y_max = 1e-2, 1.2e-1
    if indx_i == 1: y_min, y_max = 2e-2, 2.5e-1
    if indx_i == 2: y_min, y_max = 5e-2, 7e-1

  if scales == 'middle':
    x_min, x_max = 5e-3, 1e-1
    if indx_i == 0: y_min, y_max = 4e-3, 9e-2
    if indx_i == 1: y_min, y_max = 1e-2, 5e-1
    


  ax.set_xlim( x_min, x_max )
  ax.set_ylim( y_min, y_max )
  ax.set_xscale('log')
  ax.set_yscale('log')


  [sp.set_linewidth(border_width) for sp in ax.spines.values()]

  if indx_j > 0:ax.set_yticklabels([])
  if indx_i != nrows-1 :ax.set_xticklabels([])

  ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in' )
  ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')

  if indx_j == 0: ax.set_ylabel( r' $\Delta_F^2(k)$', fontsize=label_size, color= text_color )
  if indx_i == nrows-1: ax.set_xlabel( r'$ k   \,\,\,  [\mathrm{s}\,\mathrm{km}^{-1}] $',  fontsize=label_size, color= text_color )



if scales == 'middle':
  for i in range( nrows ):
    for j in range( ncols ):
      if not flags[i,j]:
        ax = ax_l[i][j].axis('off')
      

fileName = output_dir + f'flux_ps_grid_{scales}'
fileName += '.png'
# fileName += '.pdf'
fig.savefig( fileName,  pad_inches=0.1, bbox_inches='tight', dpi=fig_dpi)
print('Saved Image: ', fileName)








# nrows, ncols = 1, 1
# fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20*ncols,5*nrows))
# 
# data_mean = comparable_data['P(k)']['mean']
# data_sigma = comparable_data['P(k)']['sigma']
# n_points = len( data_mean )
# x = np.arange( 0, n_points, 1)
# 
# sim_id = 0
# sim_mean = comparable_grid[sim_id]['P(k)']['mean']
# 
# ax.errorbar( x, data_mean, yerr=data_sigma, fmt='o', c='C0', label='Data', ms=1)
# ax.scatter(x, sim_mean, color='C1', s=1 )
# 
# ax.set_yscale('log')
# 
# figure_name = output_dir + 'ps_data_sim.png'
# fig.savefig( figure_name, bbox_inches='tight', dpi=500 )
# print( f'Saved Figure: {figure_name}' )