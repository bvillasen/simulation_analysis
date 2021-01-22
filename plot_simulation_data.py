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
from plot_T0_tau import plot_T0_and_tau


ps_data_dir = 'lya_statistics/data/'
output_dir = root_dir + 'figures/'
create_directory( output_dir )

SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
sim_ids = SG.sim_ids
SG.Load_Grid_Analysis_Data( sim_ids=sim_ids, load_fit=True )

sim_data_sets = [ ]
for sim_id in sim_ids:
  sim_data  = SG.Grid[sim_id]['analysis']
  sim_param = SG.Grid[sim_id]['parameters']
  deltaZ_H = sim_param['deltaZ_H'] 
  sim_data['plot_label'] = 'P19 Mod ' + r'$\Delta z_{\mathrm{H}}$' + ' = {0:.1f}'.format( deltaZ_H )
  sim_data_sets.append(sim_data)

plot_power_spectrum_grid( ps_data_dir, output_dir, scales='small', sim_data_sets=sim_data_sets, system=system, high_z_only=True )


# plot_power_spectrum_grid( ps_data_dir, output_dir, scales='large', sim_data_sets=sim_data_sets, system=system )
# plot_power_spectrum_grid( ps_data_dir, output_dir, scales='middle', sim_data_sets=sim_data_sets, system=system )


plot_T0_and_tau( output_dir, sim_data_sets=sim_data_sets, system=system )
