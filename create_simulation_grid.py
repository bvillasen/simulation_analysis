import os, sys
import numpy as np
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from parameters_UVB_rates import param_UVB_Rates
from simulation_grid import Simulation_Grid
from simulation_parameters import *
from plot_UVB_Rates import Plot_Grid_UVB_Rates


SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Create_Grid_Directory_Structure()
SG.Create_Directories_for_Simulations()
# SG.Create_All_Submit_Job_Scripts()
SG.Create_All_Parameter_Files()
SG.Create_UVB_Rates_Files()


output_dir = root_dir + 'figures/'
create_directory( output_dir ) 

# SG.Load_Grid_UVB_Rates()
# Plot_Grid_UVB_Rates( SG, output_dir )

