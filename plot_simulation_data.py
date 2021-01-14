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
from plot_power_spectrum_functions import plot_power_spectrum_grid


output_dir = root_dir + 'figures/'
create_directory( output_dir )

SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
sim_ids = SG.sim_ids
SG.Load_Grid_Analysis_Data( sim_ids=sim_ids, load_fit=True )

sim_id = 0
sim_data  = SG.Grid[sim_id]['analysis']
sim_param = SG.Grid[sim_id]['parameters'] 


ps_data_dir = 'lya_statistics/data/'
plot_power_spectrum_grid( ps_data_dir, output_dir, scales='small', sim_data=sim_data )
plot_power_spectrum_grid( ps_data_dir, output_dir, scales='large', sim_data=sim_data )

