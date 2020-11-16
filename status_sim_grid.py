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


create_directory( root_dir )
create_directory( figures_dir )


SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Get_Grid_Status( check_queue=True )


# SG.Load_Grid_Analysis_Data( )

# SG.Load_Grid_UVB_Rates()
# Plot_Grid_UVB_Rates( SG, figures_dir )