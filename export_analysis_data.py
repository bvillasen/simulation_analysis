import os, sys
import numpy as np
import pickle
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
SG.Get_Grid_Status( check_queue=False )


SG.Load_Grid_Analysis_Data( )



sim_ids = SG.sim_ids
data_out = {}
for sim_id in sim_ids:
  data_out[sim_id] = {}
  data_out[sim_id]['z']  = SG.Grid[sim_id]['analysis']['z']
  data_out[sim_id]['T0'] = SG.Grid[sim_id]['analysis']['T0']
  data_out[sim_id]['F_mean'] = SG.Grid[sim_id]['analysis']['F_mean']
  data_out[sim_id]['parameters'] = SG.Grid[sim_id]['parameters']
  
out_file_name = root_dir + 'scale_H_T0.pkl'
f = open( out_file_name, "wb")
pickle.dump( data_out, f)
print ( f'Saved File: {out_file_name}' )
