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
# SG.Create_Grid_Directory_Structure()
# SG.Create_Directories_for_Simulations()
# SG.Create_All_Submit_Job_Scripts()
# SG.Create_All_Parameter_Files()
# SG.Create_UVB_Rates_Files()

# SG.Delete_Grid_Output_files()


# SG.Submit_Simulation_Job( 0 )

sim_ids = SG.sim_ids
for sim_id in sim_ids:
  if sim_id == 0: continue
  SG.Submit_Simulation_Job( sim_id )
# 
# SG.Get_Grid_Status()
