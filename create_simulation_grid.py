import os, sys
import numpy as np
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from parameters_UVB_rates import param_UVB_Rates
from simulation_grid import Simulation_Grid
from simulation_parameters import sim_params, job_params

# cholla_dir = '/home/bruno/cholla/'    
# root_dir = '/raid/bruno/data/cosmo_sims/sim_grid/test/'

root_dir   = '/data/groups/comp-astro/bruno/cosmo_sims/sim_grid/'
cholla_dir = '/home/brvillas/cholla/'    



create_directory( root_dir )



SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
# SG.Create_Grid_Directory_Structure()
# SG.Create_Directories_for_Simulations()
# SG.Create_All_Submit_Job_Scripts()
# SG.Create_All_Parameter_Files()
# SG.Create_UVB_Rates_Files()

# SG.Submit_Simulation_Job( 0 )

SG.Fit_Grid_Phase_Diagram( )
