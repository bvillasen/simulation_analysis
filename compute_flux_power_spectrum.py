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
from flux_power_spectrum import get_skewer_flux_power_spectrum



SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Load_Grid_Analysis_Data( sim_ids = [0], load_fit=False )

sim_id = 0
 
data_sim  = SG.Grid[sim_id]['analysis']
data_ps  = SG.Grid[sim_id]['analysis']['power_spectrum']
available_indices = data_sim['ps_available_indices'] 





sim_dir = SG.Get_Simulation_Directory( sim_id )
analysis_dir = sim_dir + 'analysis_files/'


index = 0
n_file = available_indices[index]


file_name = analysis_dir + f'{n_file}_analysis.h5'
print( f'Loading File: {file_name}' )
file = h5.File( file_name, 'r' )
current_z = file.attrs['current_z'][0]
print( f' current_z: {current_z}' )

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
k_vals_0  = k_vals_0[indices]
ps_mean_0 = ps_mean_0[indices]


vel_Hubble = vel_Hubble_x
flux_all = np.concatenate([ flux_HI_x, flux_HI_y, flux_HI_z ])


d_log_k = 0.1
flux_ps_all = []
for F_los in flux_all:
  delta_F = F_los / F_mean_HI
  k_vals, flux_power_spectrum = get_skewer_flux_power_spectrum( vel_Hubble, delta_F, d_log_k=d_log_k )
  flux_ps_all.append( flux_power_spectrum )

flux_ps_all = np.array( flux_ps_all )  
ps_mean = flux_ps_all.mean( axis=0 )  

diff_ps = np.abs( ps_mean - ps_mean_0 ) / ps_mean_0