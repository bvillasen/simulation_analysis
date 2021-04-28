import os, sys
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from parameters_UVB_rates import param_UVB_Rates
from simulation_grid import Simulation_Grid
from simulation_parameters import *
from load_tabulated_data import load_power_spectrum_table, load_data_boss


ps_data_dir = 'lya_statistics/data/'


SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
root_dir = SG.root_dir

output_dir = root_dir + 'figures/rescaled_ps_walther/'
create_directory( output_dir )



sim_ids = list(SG.sim_ids)
# sim_ids = range(10)
SG.Load_Grid_Analysis_Data( sim_ids=sim_ids, load_fit=True )

dir_boss = ps_data_dir + 'data_power_spectrum_boss/'
data_filename = dir_boss + 'data_table.py'
data_boss = load_data_boss( data_filename )
data_z_boss = data_boss['z_vals']



data_filename = ps_data_dir + 'data_power_spectrum_walther_2019/data_table.txt'
data_walther = load_power_spectrum_table( data_filename )
data_walther_z = data_walther['z_vals']

# Select only the redshifts available on the BOSS dataset
z_min = 2.2
z_ids = np.where( data_walther_z >= z_min )[0]

k_max = 0.1

rescale_values = {}

# z_ids = [8]
for indx, z_id in enumerate(z_ids):

  ps_walther = data_walther[z_id]
  walther_z = ps_walther['z']
  walther_kvals = ps_walther['k_vals']
  k_indices = np.where( walther_kvals <= k_max )
  walther_kvals = walther_kvals[k_indices]
  walther_delta_ps = ps_walther['delta_power'][k_indices]
  walther_delta_ps_sigma = ps_walther['delta_power_error'][k_indices]

  k_max_compare = 0.02
  compare_indices = np.where( walther_kvals <= k_max_compare )
  walther_kvals_comapre = walther_kvals[compare_indices]
  walther_delta_ps_compare = walther_delta_ps[compare_indices]
  walther_delta_ps_sigma_compare = walther_delta_ps_sigma[compare_indices]



  z_diff_boss = np.abs( data_z_boss - walther_z )
  z_id_boss = np.where( z_diff_boss == z_diff_boss.min() )[0][0]
  if z_diff_boss.min() > 0.05: print( f'WARNING: Z  Boss doesnt match: {z_diff_boss.min()} ' )
  boss_z = data_z_boss[z_id_boss]
  boss_kvals = data_boss[z_id_boss]['k_vals']
  boss_delta_ps = data_boss[z_id_boss]['delta_power']
  boss_delta_ps_sigma = data_boss[z_id_boss]['delta_power_error']






  grid_data = {}

  sim_id = 0
  for sim_id in sim_ids:
    sim_ps = SG.Grid[sim_id]['analysis']['power_spectrum']
    sim_z_vals = sim_ps['z']
    z_diff = np.abs( sim_z_vals - walther_z ) 
    if z_diff.min() > 0.01: print( f'WARNING: Large z difference: {z_diff.min()}')
    sim_z_id = np.where( z_diff == z_diff.min()  )[0][0]
    sim_z = sim_ps['z'][sim_z_id]
    sim_power = sim_ps['ps_mean'][sim_z_id]
    sim_kvals = sim_ps['k_vals'][sim_z_id]
    sim_delta_ps = sim_power * sim_kvals / np.pi

    # Interpolate the sim ps to the walther ps k_vals
    sim_delta_ps_interp = np.interp( walther_kvals_comapre, sim_kvals, sim_delta_ps )
    # Get Chi Squared
    chi2 =  np.sum( ( ( sim_delta_ps_interp - walther_delta_ps_compare ) / walther_delta_ps_sigma_compare )**2  )
    grid_data[sim_id] = {}
    grid_data[sim_id]['chi2'] = chi2
    grid_data[sim_id]['k_vals'] = sim_kvals
    grid_data[sim_id]['delta_ps'] = sim_delta_ps
    grid_data[sim_id]['delta_ps_interp'] = sim_delta_ps_interp


  chi2_vals = np.array([ grid_data[sim_id]['chi2'] for sim_id in sim_ids ] )
  min_sim_id = np.where( chi2_vals == chi2_vals.min() )[0][0]

  sim_data = grid_data[min_sim_id]

  sim_kvals    = sim_data['k_vals']
  sim_delta_ps = sim_data['delta_ps']

  # Rescale the simulation P(k) to minimize chi2 with respect to BOSS P(k)
  compare_indices = np.where( boss_kvals <= k_max_compare )
  boss_kvals_comapre          = boss_kvals[compare_indices]
  boss_delta_ps_compare       = boss_delta_ps[compare_indices]
  boss_delta_ps_sigma_compare = boss_delta_ps_sigma[compare_indices]

  def get_Chi2_rescaled( alpha, ps, ps_goal, ps_goal_sigma ):
    ps_rescaled = alpha * ps
    chi2 =  np.sum( ( ( ps_rescaled - ps_goal ) / ps_goal_sigma )**2  )
    return chi2



  # Interpolate the sim ps to the walther ps k_vals
  sim_delta_ps_interp = np.interp( boss_kvals_comapre, sim_kvals, sim_delta_ps )

  from scipy import optimize

  min_result = optimize.minimize( get_Chi2_rescaled, x0=1, args=(sim_delta_ps_interp, boss_delta_ps_compare, boss_delta_ps_sigma_compare) )
  rescale_alpha = min_result.x[0]
  
  rescale_values[z_id] = {}
  rescale_values[z_id]['z'] = walther_z
  rescale_values[z_id]['alpha'] = rescale_alpha

  import pylab
  import matplotlib
  matplotlib.rcParams['mathtext.fontset'] = 'cm'
  matplotlib.rcParams['mathtext.rm'] = 'serif'


  if system == 'Lux':      prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/brvillas/fonts', "Helvetica.ttf"), size=12)
  if system == 'Shamrock': prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=12)

  fig_height = 7
  fig_width = 8

  label_size = 18
  figure_text_size = 18
  legend_font_size = 16
  tick_label_size_major = 15
  tick_label_size_minor = 13
  tick_size_major = 5
  tick_size_minor = 3
  tick_width_major = 1.5
  tick_width_minor = 1


  c_walther = pylab.cm.viridis(.3)
  # c_boss = pylab.cm.Purples(.7)
  c_boss = 'C1'
  
  x_min, x_max = 2e-3, 0.11
  if z_id == 2 : y_min, y_max = 5e-3, 5e-2 
  if z_id == 3 : y_min, y_max = 6e-3, 6e-2 
  if z_id == 4 : y_min, y_max = 5e-3, 7e-2 
  if z_id == 5 : y_min, y_max = 7e-3, 8e-2 
  if z_id == 6 : y_min, y_max = 1e-2, 1e-1 
  if z_id == 7 : y_min, y_max = 1e-2, 1.5e-1 
  if z_id == 8 : y_min, y_max = 1.5e-2, 2e-1 

  ncols, nrows = 2, 1
  fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols*fig_width,nrows*fig_height))
  plt.subplots_adjust( hspace = 0.02, wspace=0.02)

  walther_k = walther_kvals
  walther_delta_power = walther_delta_ps
  walther_delta_power_error = walther_delta_ps_sigma
  label_walther ='Walther et al. (2018)' 
  # 
  # 
  # ax = ax_l[0]
  # 
  # ax.errorbar( walther_k, walther_delta_power, yerr=walther_delta_power_error, fmt='o', c=c_walther, label=label_walther, zorder=2)
  # 
  # sim_k = sim_data['k_vals']
  # sim_delta_ps = sim_data['delta_ps']
  # ax.plot( sim_k, sim_delta_ps )
  # 
  # 
  # x_fill = [ k_max_compare, 10 ]
  # y_fill_max = [ 1, 1 ]
  # y_fill_min = [ 0, 0 ]
  # ax.fill_between( x_fill, y_fill_min, y_fill_max, color='gray', alpha=0.3 )
  # 
  # ax.text(0.9, 0.95, r'$z={0:.1f}$'.format(walther_z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size ) 
  # 
  # ax.set_ylabel( r' $\Delta_F^2(k)$', fontsize=label_size )
  # ax.set_xlabel( r'$ k   \,\,\,  [\mathrm{s}\,\mathrm{km}^{-1}] $',  fontsize=label_size  )
  # 
  # 
  # legend_loc = 3
  # ax.legend(  loc=legend_loc, frameon=False, prop=prop    )
  # 
  # ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in' )
  # ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')
  # 
  # ax.set_xscale('log')
  # ax.set_yscale('log')
  # 
  # 
  # ax.set_xlim( x_min, x_max )
  # ax.set_ylim( y_min, y_max )
  # 
  # 




  ax = ax_l[0]

  boss_k = boss_kvals
  boss_delta_power = boss_delta_ps
  boss_delta_power_error = boss_delta_ps_sigma
  label_boss = 'eBOSS (2019)'
  ax.errorbar( boss_k, boss_delta_power, yerr=boss_delta_power_error, fmt='o', c=c_boss, label=label_boss, zorder=2)

  ax.errorbar( walther_k, walther_delta_power, yerr=walther_delta_power_error, fmt='o', c=c_walther, label=label_walther, zorder=2)

  sim_k = sim_data['k_vals']
  sim_delta_ps = sim_data['delta_ps']
  ax.plot( sim_k, sim_delta_ps )


  x_fill = [ k_max_compare, 10 ]
  y_fill_max = [ 1, 1 ]
  y_fill_min = [ 0, 0 ]
  ax.fill_between( x_fill, y_fill_min, y_fill_max, color='gray', alpha=0.3 )

  ax.text(0.9, 0.95, r'$z={0:.1f}$'.format(boss_z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size ) 

  ax.set_ylabel( r' $\Delta_F^2(k)$', fontsize=label_size )
  ax.set_xlabel( r'$ k   \,\,\,  [\mathrm{s}\,\mathrm{km}^{-1}] $',  fontsize=label_size  )

  legend_loc = 3
  ax.legend(  loc=legend_loc, frameon=False, prop=prop    )

  ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in' )
  ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')

  # ax.tick_params(axis='y', which='both', labelsize=0, labelcolor='white' )


  ax.set_xscale('log')
  ax.set_yscale('log')

  ax.set_xlim( x_min, x_max )
  ax.set_ylim( y_min, y_max )





  ax = ax_l[1]

  boss_k = boss_kvals
  boss_delta_power = boss_delta_ps
  boss_delta_power_error = boss_delta_ps_sigma
  label_boss = 'eBOSS (2019)'
  ax.errorbar( boss_k, boss_delta_power, yerr=boss_delta_power_error, fmt='o', c=c_boss, label=label_boss, zorder=2)

  label_walther = label_walther + ' Rescaled'
  ax.errorbar( walther_k, rescale_alpha* walther_delta_power, yerr=walther_delta_power_error, fmt='o', c=c_walther, label=label_walther, zorder=2)

  sim_k = sim_data['k_vals']
  sim_delta_ps = sim_data['delta_ps']
  ax.plot( sim_k, rescale_alpha * sim_delta_ps )


  x_fill = [ k_max_compare, 10 ]
  y_fill_max = [ 1, 1 ]
  y_fill_min = [ 0, 0 ]
  ax.fill_between( x_fill, y_fill_min, y_fill_max, color='gray', alpha=0.3 )

  ax.text(0.9, 0.95, r'$z={0:.1f}$'.format(boss_z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size ) 

  ax.text(0.1, 0.95, r'$\alpha={0:.2f}$'.format(rescale_alpha), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size ) 

  # ax.set_ylabel( r' $\Delta_F^2(k)$', fontsize=label_size )
  ax.set_xlabel( r'$ k   \,\,\,  [\mathrm{s}\,\mathrm{km}^{-1}] $',  fontsize=label_size  )

  legend_loc = 3
  ax.legend(  loc=legend_loc, frameon=False, prop=prop    )

  ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in' )
  ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')

  ax.tick_params(axis='y', which='both', labelsize=0, labelcolor='white' )


  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.set_xlim( x_min, x_max )
  ax.set_ylim( y_min, y_max )


  file_name =  output_dir + f'rescaled_ps_walther_{z_id}.png'
  fig.savefig( file_name,  pad_inches=0.1, bbox_inches='tight', dpi=300)
  print('Saved Image: ', file_name)



out_file_name = ps_data_dir + 'rescale_walther_to_boss.pkl'
Write_Pickle_Directory( rescale_values,  out_file_name )