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

check_queue = True
if system == 'Shamrock': check_queue = False

SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
# SG.Get_Grid_Status( check_queue=check_queue )


params = { 'scale_He':None, 'deltaZ_He':None, 'scale_H':0.86, 'deltaZ_H':0.0 }
tolerance = 5e-3

sim_ids = SG.sim_ids
vals_to_find = { key:params[key] for key in params if params[key] != None }

selected_sims = []
for sim_id in sim_ids: 
  sim_data = SG.Grid[sim_id]
  sim_parameters = sim_data['parameters']
  sim_vals = { key:sim_parameters[key] for key in vals_to_find }
  same_vals = True
  for key in sim_vals:
    if np.abs( sim_vals[key] - vals_to_find[key] ) > tolerance: same_vals = False
  if same_vals: selected_sims.append( sim_id )


# SG.Split_Grid_Hydro_Particles_Sanpshots( )

# SG.Load_Grid_Analysis_Data( )

# SG.Load_Grid_UVB_Rates()
# Plot_Grid_UVB_Rates( SG, figures_dir )