import os, sys
import numpy as np
import h5py as h5
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






use_mpi = False
if use_mpi:
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  n_procs = comm.Get_size()
else:
  rank = 0
  n_procs = 1
  

sim_ids = range(400)
sim_ids = [0]

sim_ids_proc = split_indices( sim_ids, rank, n_procs, adjacent=False )

print( f'Id: {rank} sim_ids:{sim_ids_proc}' )  


# sim_id = 0
for sim_id in sim_ids_proc:  


  SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
  SG.Load_Grid_Analysis_Data( sim_ids = [sim_id], load_fit=False )
  ps_out_dir = SG.root_dir + 'flux_power_spectrum_files/'
  # if rank == 0: create_directory( ps_out_dir )

  # Normaliz_Flux_Power_Spectrum( sim_id, ps_out_dir, SG )

  sim_name = SG.Grid[sim_id]['key']
  output_dir = ps_out_dir + sim_name + '/'
  create_directory( output_dir )


  data_sim  = SG.Grid[sim_id]['analysis']
  data_ps  = SG.Grid[sim_id]['analysis']['power_spectrum']
  available_indices = data_sim['ps_available_indices'] 





  sim_dir = SG.Get_Simulation_Directory( sim_id )
  analysis_dir = sim_dir + 'analysis_files/'

  # index =  0
  # n_file = available_indices[index]
  
  available_indices = [55]
  for n_file in available_indices:

    output_file_name = output_dir + f'flux_ps_{n_file}.h5'
    file_exists = check_if_file_exists( output_file_name )
    # if file_exists:
      # print ( f'Skipping File: {output_file_name} ' )
      # continue
    
    outfile = h5.File( output_file_name, 'w' )

    file_name = analysis_dir + f'{n_file}_analysis.h5'
    print( f'Loading File: {file_name}' )
    file = h5.File( file_name, 'r' )
    current_z = file.attrs['current_z'][0]
    print( f' current_z: {current_z}' )

    outfile.attrs['current_z'] = current_z

    lya_data = file['lya_statistics']
    ps_data = lya_data['power_spectrum']


    skewers_key = 'skewers_x'

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
    flux_all = np.concatenate([ flux_HI_x, flux_HI_y, flux_HI_z ])


    tau_eff_simulation = -np.log(F_mean_HI[0])
    tau_eff_Becker =  Compute_analytical_TauEff_Becker( current_z )

    print(f' Tau Effective Simulation: {tau_eff_simulation}' )
    print(f' Tau Effective Becker:     {tau_eff_Becker}' )



    type = 'normalize_F_mean'
    type = 'normalize_tau_eff'
    # types = [ 'normalize_F_mean', 'normalize_tau_eff' ] 
    # types = [  'normalize_tau_eff' ]
    # types = [  'normalize_F_mean' ] 


    types = [ 'tau_eff_local', 'tau_eff_global',  ]
    
    type = 'tau_eff_local' 

    tau_eff_vals = { 'Simulation': tau_eff_simulation,  'Becker': tau_eff_Becker }

    for normalization in tau_eff_vals:

      print( f'Normalization: {normalization}' )

      group = outfile.create_group( normalization )

      tau_eff = tau_eff_vals[normalization]
      group.attrs['tau_eff'] = tau_eff


      d_log_k = 0.1
      F_min = 1e-100
      F_mean_all = []
      for type in types:
        print( f' Type: {type} ')
        flux_ps_all = []
        for F_los in flux_all:
          if type == 'tau_eff_local' or type == 'tau_eff_global':
            F_los[ F_los < F_min ] = F_min 
            tau_los  = -np.log( F_los )
            if type == 'tau_eff_local':  tau_los =  tau_los / tau_los.mean() * tau_eff
            if type == 'tau_eff_global': tau_los =  tau_los / tau_eff_simulation * tau_eff
            F_los = np.exp( - tau_los )
            
          F_mean_all.append( F_los.mean() )
          F_mean = np.exp( -tau_eff )
          # print( F_los.mean() )
          delta_F = ( F_los - F_mean ) / F_mean
          k_vals, flux_power_spectrum = get_skewer_flux_power_spectrum( vel_Hubble, delta_F, d_log_k=d_log_k )
          flux_ps_all.append( flux_power_spectrum )

        F_mean_all = np.array( F_mean_all )
        flux_ps_all = np.array( flux_ps_all )  
        ps_mean = flux_ps_all.mean( axis=0 )  

        diff_ps = np.abs( ps_mean - ps_mean_sim ) / ps_mean_sim
        print( diff_ps )
        
        
        print( f'F_mean  comparison: {F_mean}   {F_mean_all.mean()}' )
        print( f'tau_eff comparison: {tau_eff}   {-np.log(F_mean_all.mean() ) }' )
        
        
        group.create_dataset( type, data=ps_mean )
        print( f'Saved normalization: {normalization}')

    outfile.create_dataset( 'k_vals', data=k_vals )
    outfile.close()
    print( f'Saved File: {output_file_name}' )
    # comm.Barrier()
    
    






