import os, sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from parameters_UVB_rates import param_UVB_Rates
from simulation_grid import Simulation_Grid
from simulation_parameters import *
from mcmc_functions import *
from mcmc_data_functions import *
from data_thermal_history import *



SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Load_Grid_Analysis_Data()
sim_ids = SG.sim_ids

comparable_data = Get_Comparable_Composite_T0_tau()
comparable_grid = Get_Comparable_Composite_T0_tau_from_Grid( comparable_data, SG )


parameters = SG.parameters
param_id = 0
param_name = parameters[param_id]['name']

param_value = 0.45
param_vals = parameters[param_id]['values']
param_max = max(param_vals)
param_min = min(param_vals) 

#find closest simulations 
sim_ids = SG.sim_ids
sim_id_l, sim_id_r = 0, 0
diff_l, diff_r = -np.inf, np.inf
for sim_id in sim_ids:
  sim_params = SG.Grid[sim_id]
  sim_param_val = sim_params['parameters'][param_name]
  diff = sim_param_val - param_value
  if diff > 0 and diff < diff_r: sim_id_r, diff_r = sim_id, diff
  if diff < 0 and diff > diff_l: sim_id_l, diff_l = sim_id, diff  