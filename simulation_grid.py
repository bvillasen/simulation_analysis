import os, sys
import numpy as np
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from parameters_UVB_rates import param_UVB_Rates
from simulation_grid import Simulation_Grid
from simulation_parameters import *

create_directory( root_dir )


SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )

# SG.Submit_Simulation_Job( 0 )

# SG.Get_Grid_Status()


# SG.Fit_Grid_Phase_Diagram( )
