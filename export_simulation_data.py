import os, sys
import numpy as np
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from parameters_UVB_rates import param_UVB_Rates
from simulation_grid import Simulation_Grid
from simulation_parameters import *
from load_tabulated_data import *
from mcmc_data_functions import Get_Comparable_T0_Gaikwad



SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Load_Grid_Analysis_Data( )
sim_ids = SG.sim_ids
output_dir = SG.root_dir + 'simulation_grid_data/'
create_directory( output_dir )

grid_param_values = np.array([ SG.Grid[sim_id]['parameter_values'] for sim_id in sim_ids ])
header = f'Row_i = Values of [ beta_He, beta_H, delta_z_He, delta_z_H ] for simulation i \nn_rows, n_cols = {grid_param_values.shape[0]}, {grid_param_values.shape[1]}  '
file_name = output_dir + f'grid_parameters.txt'
np.savetxt( file_name, grid_param_values, header=header, fmt='%.2f', )
print( f'Saved File: {file_name}' )


z_vals = SG.Grid[0]['analysis']['z']

fields = [ 'T0', 'tau', 'tau_HeII' ]
data_all = {}
for field in fields:
  field_vals_all = []
  for sim_id in sim_ids:
    field_vals_all.append( SG.Grid[sim_id]['analysis'][field]  ) 
  data_all[field] = np.array( field_vals_all )


for field in fields:
  data_field = data_all[field]
  data_out = np.concatenate([ z_vals[:,None].T, data_field ], axis=0 )
  header = f'Row_0 = Redshift values     Row_i = Evolution of {field} for simulation i-1 \nn_rows, n_cols = {data_out.shape[0]}, {data_out.shape[1]} '
  file_name = output_dir + f'grid_{field}_evolution.txt'
  np.savetxt( file_name, data_out, header=header )
  print( f'Saved File: {file_name}' )



z_vals = [ 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 5.0  ]


for z in z_vals:

  ps_mean_data = []
  for sim_id in sim_ids:
    sim_ps_data = SG.Grid[sim_id]['analysis']['power_spectrum']
    sim_z = sim_ps_data['z']

    z_diff = np.abs( sim_z - z )
    z_diff_min = z_diff.min()
    index = np.where( z_diff == z_diff_min )[0]
    if len( index ) > 1: print("# WARNING:  more than one redshift index " )
    if z_diff_min > 0.05: print("# WARNING:  Large Redshift difference" )
    k_vals  = sim_ps_data['k_vals'][index[0]]
    ps_mean = sim_ps_data['ps_mean'][index[0]] 
    ps_mean_data.append( ps_mean )
  ps_mean_data = np.array( ps_mean_data )
  data_ps = np.concatenate([ k_vals[:,None].T, ps_mean_data ])
  file_name = output_dir + f'grid_power_spectrum_z={z:.1f}.txt'
  header = f'Row_0 = k values [ s/km ]     Row_i = P(k) for simulation  i-1  \nn_rows, n_cols = {data_ps.shape[0]}, {data_ps.shape[1]}'
  np.savetxt( file_name, data_ps, header=header )
  print( f'Saved File: {file_name} ')



output_data_dir = output_dir + 'data/'
create_directory( output_data_dir )



data_T0 = Get_Comparable_T0_Gaikwad()
data_out = np.array([ data_T0['z'], data_T0['mean'], data_T0['sigma'] ])
header = f'Row_0 = Redshift values      Row_1 = T0      Row_2 = sigma_T0   \nn_rows, n_cols = {data_out.shape[0]}, {data_out.shape[1]}'
file_name = output_data_dir + f'data_T0_Gaikwad.txt'
np.savetxt( file_name, data_out, header=header )
print( f'Saved File: {file_name} ')



# 
# 
# ps_data_dir = 'lya_statistics/data/'
# output_ps_data_dir = output_data_dir + 'power_spectrum'
# 
# z_vals = [ 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 5.0  ]
# 
# dir_boss = ps_data_dir + 'data_power_spectrum_boss/'
# data_filename = dir_boss + 'data_table.py'
# data_boss = load_data_boss( data_filename )
# 
# data_filename = ps_data_dir + 'data_power_spectrum_walther_2019/data_table.txt'
# data_walther = load_power_spectrum_table( data_filename )
# 
# # dir_data_boera = ps_data_dir + 'data_power_spectrum_boera_2019/'
# # data_boera = load_tabulated_data_boera( dir_data_boera )
# # data_z_b = data_boera['z_vals']
# # 
# # data_dir_viel = ps_data_dir + 'data_power_spectrum_viel_2013/'
# # data_viel = load_tabulated_data_viel( data_dir_viel)
# # data_z_v = data_viel['z_vals']
# 
# 
# 
# 
# data_power_spectrum_all = { 'Boss':data_boss, 'Walther':data_walther }
# 
# 
# 
# for data_name in data_power_spectrum_all:
# 
#   print( f'Exporting data power spectrum: {data_name}')
#   out_dir = f'{output_ps_data_dir}_{data_name}/' 
#   create_directory( out_dir )
#   data_ps = data_power_spectrum_all[data_name]
# 
#   data_z = data_ps['z_vals']
# 
#   for z in z_vals:
#     z_diff = np.abs( data_z - z )
#     z_diff_min = z_diff.min()
#     index = np.where( z_diff == z_diff_min )[0]
#     if len( index ) > 1: print("# WARNING:  more than one redshift index " )
#     if z_diff_min > 0.05: 
#       print( f' Skipping this redshift: {z}')
#       continue
#     index = index[0]
#     k_vals = data_ps[index]['k_vals']
#     power_spectrum = data_ps[index]['power_spectrum']
#     sigma_power_spectrum = data_ps[index]['sigma_power_spectrum']
#     data_out = np.array([ k_vals, power_spectrum, sigma_power_spectrum ])
#     file_name =  out_dir + f'data_z={z:.1f}.txt'
#     header = f'Row_0 = k values [ s/km ]     Row_1 = P(k)     Row_2 = sigma_P(k)   \nn_rows, n_cols = {data_out.shape[0]}, {data_out.shape[1]}'
#     np.savetxt( file_name, data_out, header=header )
#     print( f'Saved File: {file_name} ')
# 
# 
