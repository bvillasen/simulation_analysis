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
from normalize_power_spectrum  import Normaliz_Flux_Power_Spectrum
from data_optical_depth import Compute_analytical_TauEff_Becker
from flux_power_spectrum import get_skewer_flux_power_spectrum






def ReScale_Optical_Depth_To_F_Mean_Diff( alpha, F_mean, tau_los  ):
  # print(alpha)
  tau_los_rescaled = tau_los * alpha
  F_los_rescaled = np.exp( - tau_los_rescaled )
  F_mean_rescaled = F_los_rescaled.mean()
  diff = F_mean_rescaled - F_mean
  return diff


def ReScale_Optical_Depth_To_F_Mean( tau_los, F_mean ):
  from scipy import optimize
  guess = tau_eff / tau_los.mean()
  alpha = optimize.newton(ReScale_Optical_Depth_To_F_Mean_Diff, guess, args=(F_mean, tau_los ) ) 
  tau_los_rescaled = alpha * tau_los
  F_los_rescaled = np.exp( - tau_los_rescaled )
  F_mean_rescaled = F_los_rescaled.mean()
  diff = np.abs( F_mean_rescaled - F_mean ) / F_mean
  if diff > 1e-6: print( 'WARNING: Rescaled F_mean mismatch: {F_mean_rescaled}   {f_mean}')
  return  tau_los_rescaled




use_mpi = True
if use_mpi:
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  n_procs = comm.Get_size()
else:
  rank = 0
  n_procs = 1
  

# sim_ids = [0]
sim_ids = range(400)

sim_ids_proc = split_indices( sim_ids, rank, n_procs, adjacent=False )
print( f'Id: {rank} sim_ids:{sim_ids_proc}' )  

# sim_id = 0
for sim_id in sim_ids_proc:

  SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
  SG.Load_Grid_Analysis_Data( sim_ids = [sim_id], load_fit=False )
  ps_out_dir = SG.root_dir + 'flux_power_spectrum_files/'

  sim_name = SG.Grid[sim_id]['key']
  output_dir = ps_out_dir + sim_name + '/'
  create_directory( output_dir )

  data_sim  = SG.Grid[sim_id]['analysis']
  data_ps  = SG.Grid[sim_id]['analysis']['power_spectrum']
  available_indices = data_sim['ps_available_indices'] 
  # available_indices = [30]

  sim_dir = SG.Get_Simulation_Directory( sim_id )
  analysis_dir = sim_dir + 'analysis_files/'

  # n_file = 55
  for n_file in available_indices:

    output_file_name = output_dir + f'flux_ps_{n_file}.h5'
    file_exists = check_if_file_exists( output_file_name )
    outfile = h5.File( output_file_name, 'w' )

    file_name = analysis_dir + f'{n_file}_analysis.h5'
    print( f'Loading File: {file_name}' )
    file = h5.File( file_name, 'r' )
    current_z = file.attrs['current_z'][0]
    print( f' current_z: {current_z}' )
    outfile.attrs['current_z'] = current_z

    lya_data = file['lya_statistics']
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

    n_skewers = lya_data.attrs['n_skewers']
    F_mean_HI = lya_data.attrs['Flux_mean_HI']

    k_vals_0 = ps_data['k_vals'][...]
    ps_mean_0 = ps_data['p(k)'][...]
    indices = ps_mean_0 > 1e-10
    k_vals_sim  = k_vals_0[indices]
    ps_mean_sim = ps_mean_0[indices]

    vel_Hubble = vel_Hubble_x
    F_los_all = np.concatenate([ flux_HI_x, flux_HI_y, flux_HI_z ])

    tau_eff_simulation = -np.log(F_mean_HI[0])
    tau_eff_Becker =  Compute_analytical_TauEff_Becker( current_z )

    # print(f' Tau Effective Simulation: {tau_eff_simulation}' )
    # print(f' Tau Effective Becker:     {tau_eff_Becker}' )
    tau_eff = tau_eff_simulation
    F_mean = np.exp( -tau_eff )

    # Substitute the F=0 pixels 
    F_min = 1e-200
    F_los_all[ F_los_all < F_min ] = F_min

    tau_eff_rescale = tau_eff_Becker
    F_mean_rescale = np.exp( -tau_eff_rescale )

    tau_los_all = -np.log( F_los_all )
    tau_los_all_rescaled = ReScale_Optical_Depth_To_F_Mean( tau_los_all, F_mean_rescale )
    F_los_all_rescaled = np.exp( -tau_los_all_rescaled )
    F_mean_rescaled = F_los_all_rescaled.mean()

    print( f'F_mean Becker:       {F_mean_rescale }' )
    print( f'Rescaling to F_mean: {F_mean_rescaled }' )

    normalization = 'Becker'
    type = 'tau_eff_global'
    print( f'Normalization: {normalization}   {type}' )

    d_log_k = 0.1
    # Compute the power spectrum from the rescaled 
    flux_ps_all_rescaled = []
    for F_los in F_los_all_rescaled:
      delta_F = ( F_los - F_mean_rescaled ) / F_mean_rescaled
      k_vals, flux_power_spectrum = get_skewer_flux_power_spectrum( vel_Hubble, delta_F, d_log_k=d_log_k )
      flux_ps_all_rescaled.append( flux_power_spectrum )
      
    flux_ps_all_rescaled = np.array( flux_ps_all_rescaled )
    ps_mean_rescaled = flux_ps_all_rescaled.mean( axis=0 )

    group = outfile.create_group( normalization )
    group.attrs['F_mean'] = F_mean_rescale
    group.attrs['tau_eff'] = tau_eff_rescale
    group.create_dataset( type, data=ps_mean_rescaled )
    print( f'Saved Normalization: {normalization}   {type}' )

    outfile.create_dataset( 'k_vals', data=k_vals )
    outfile.close()
    print( f'Saved File: {output_file_name}' )
