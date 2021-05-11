import os, sys
import numpy as np
import h5py as h5
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from flux_power_spectrum import get_skewer_flux_power_spectrum
from tools import *
from stats_functions import compute_distribution, get_highest_probability_interval


n_points = 1024
L_MPC = 50

# n_points = 2048
# L_MPC = 100

data_dir = '/raid/bruno/data/'
input_dir  = data_dir + f'cosmo_sims/rescaled_P19/{n_points}_{L_MPC}Mpc/analysis_files/'
output_dir = data_dir + f'cosmo_sims/rescaled_P19/{n_points}_{L_MPC}Mpc/analysis_files/ps_statistics/'
create_directory( output_dir )

files_list = get_files_names( input_dir )
files_list.sort()

Lbox = L_MPC * 1e3 #kpc/h

n_files = len( files_list )

for n_file in range(n_files):
  # n_file = int(file_name.split('_')[0])
  file_name = f'{n_file}_analysis.h5'  

  print( f'\nn_file: {n_file}/{n_files}')

  in_file_name = input_dir + file_name
  print( f'Loading File: {in_file_name} ' )
  in_file = h5.File( in_file_name, 'r' ) 
  current_z = in_file.attrs['current_z'][0]
  print( f'current_z: {current_z}' )

  if current_z > 6.0: continue


  lya_data = in_file['lya_statistics']
  ps_data = lya_data['power_spectrum']

  skewers_data = lya_data['skewers_x']
  vel_Hubble_x = skewers_data['vel_Hubble'][...]
  flux_HI_x    = skewers_data['los_transmitted_flux_HI'][...]

  skewers_data = lya_data['skewers_y']
  vel_Hubble_y = skewers_data['vel_Hubble'][...]
  flux_HI_y    = skewers_data['los_transmitted_flux_HI'][...]

  skewers_data = lya_data['skewers_z']
  vel_Hubble_z = skewers_data['vel_Hubble'][...]
  flux_HI_z    = skewers_data['los_transmitted_flux_HI'][...]

  n_skewers = lya_data.attrs['n_skewers'][0]
  F_mean_HI = lya_data.attrs['Flux_mean_HI'][0]

  k_vals_sim  = ps_data['k_vals'][...]
  ps_mean_sim = ps_data['p(k)'][...]
  indices = ps_mean_sim > 0
  k_vals_sim = k_vals_sim[indices]
  ps_mean_sim = ps_mean_sim[indices]

  vel_Hubble = vel_Hubble_x
  flux_all = np.concatenate([ flux_HI_x, flux_HI_y, flux_HI_z ])

  d_log_k = 0.1
  F_min = 1e-100

  flux_ps_all = []
  for F_los in flux_all:
    F_los[ F_los < F_min ] = F_min 
    delta_F = ( F_los - F_mean_HI ) / F_mean_HI
    k_vals, flux_power_spectrum = get_skewer_flux_power_spectrum( vel_Hubble, delta_F, d_log_k=d_log_k )
    flux_ps_all.append( flux_power_spectrum )

  flux_ps_all = np.array( flux_ps_all )  
  ps_mean = flux_ps_all.mean( axis=0 )

  fill_sum = 0.67
  sigma_vals, lower_vals, higher_vals, max_vals, mean_vals = [], [], [], [], []
  for ps_slice in flux_ps_all.T:
    v_mean = ps_slice.mean()
    sigma = ps_slice.std()
    distribution, bin_centers = compute_distribution( ps_slice, 40, log=True )
    v_l, v_r, v_max, sum = get_highest_probability_interval( bin_centers, distribution, fill_sum, log=False, n_interpolate=500)
    sigma_vals.append( sigma )
    lower_vals.append( v_l )
    higher_vals.append( v_r )
    max_vals.append( v_max )
    mean_vals.append( v_mean )
  sigma_vals = np.array( sigma_vals )
  lower_vals = np.array( lower_vals )
  higher_vals = np.array( higher_vals )  
  max_vals = np.array( max_vals )  
  mean_vals = np.array( mean_vals )
    
  n_points = len( vel_Hubble )
  dv = vel_Hubble[1] - vel_Hubble[0]
  v_max = n_points * dv

  n_ind_vals = []
  for k in k_vals:
    v = 2*np.pi / k
    stride =  v / dv 
    n_ind =  ( n_points / stride )**2 #number of independent skewers 
    n_ind_vals.append( n_ind )
  n_independent = np.array( n_ind_vals )

  out_file_name = output_dir + f'{n_file}_stats.h5'
  out_file = h5.File( out_file_name, 'w' )
  for key in in_file.attrs.keys():
    out_file.attrs[key] = in_file.attrs[key]

  out_file.create_dataset( 'vel_hubble', data=vel_Hubble )  
  out_file.create_dataset( 'n_independent', data=n_independent )  
  out_file.create_dataset( 'k_vals', data=k_vals )
  out_file.create_dataset( 'mean', data=ps_mean )
  out_file.create_dataset( 'sigma', data=sigma_vals )
  out_file.create_dataset( 'higher', data=higher_vals )
  out_file.create_dataset( 'lower', data=lower_vals )  
  out_file.create_dataset( 'max', data=max_vals )

  out_file.close()
  print( f'Saved File: {out_file_name}' )
