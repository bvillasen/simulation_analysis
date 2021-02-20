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


ps_data_dir = 'lya_statistics/data/'
output_dir = root_dir + 'figures/'
create_directory( output_dir )

SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
sim_ids = SG.sim_ids
SG.Load_Grid_Analysis_Data( sim_ids=sim_ids, load_fit=True )

sim_data_sets = [ ]
# sim_ids_to_plot = [ 1 ]
sim_ids_to_plot = sim_ids
for sim_id in sim_ids_to_plot:
  sim_data  = SG.Grid[sim_id]['analysis']
  sim_param = SG.Grid[sim_id]['parameters']
  # beta_He = sim_param['scale_He'] 
  # beta_H  = sim_param['scale_H'] 
  # deltaZ_He = sim_param['deltaZ_He'] 
  # label = r'$\beta_{\mathrm{He}}:$' + f'{beta_He:.1f}' + ' ' + r'$\beta_{\mathrm{H}}:$' + f'{beta_H:.1f}' + ' ' + r'$\Delta_z:$' + f'{deltaZ_He:.1f}' 
  label = 'P19 Modified'
  sim_data['plot_label'] = label
  sim_data_sets.append(sim_data)

plot_power_spectrum_grid( ps_data_dir, output_dir, scales='small', sim_data_sets=sim_data_sets, system=system )
plot_power_spectrum_grid( ps_data_dir, output_dir, scales='large', sim_data_sets=sim_data_sets, system=system )
plot_power_spectrum_grid( ps_data_dir, output_dir, scales='middle', sim_data_sets=sim_data_sets, system=system )




plot_T0_and_tau( output_dir, sim_data_sets=sim_data_sets, system=system )

plot_tau_HeII( output_dir, sim_data_sets=sim_data_sets, system=system )
