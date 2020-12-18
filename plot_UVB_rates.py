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
from mcmc_plotting_functions import *
from plot_UVB_Rates import Plot_Grid_UVB_Rates


field = 'T0+tau'
output_dir = root_dir + f'fit_results_{field}/'
create_directory( output_dir )



SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Load_Grid_Analysis_Data()
sim_ids = SG.sim_ids

comparable_data = Get_Comparable_Composite_T0_tau( factor_sigma_tau_becker=6.0, factor_sigma_tau_keating=4.0,  )
comparable_grid = Get_Comparable_Composite_T0_tau_from_Grid( comparable_data, SG )

fields = [ 'T0', 'tau' ]
data_grid = Get_Data_Grid( fields, SG )

params = SG.parameters

deltaZ_H_vals = [ -0.45,  -0.15, 0.15, 0.45 ]
ids_closest = []
for deltaZ_H in deltaZ_H_vals:
  params_search = np.array([ 0.41, 0.82, 0.13, deltaZ_H ]) 
  id_closest = SG.Find_Closest_Simulation( params_search )
  p_vals_closest = SG.Grid[id_closest]['parameter_values']
  ids_closest.append( id_closest )



output_dir = root_dir + 'figures/'
create_directory( output_dir ) 
SG.Load_Grid_UVB_Rates()
Plot_Grid_UVB_Rates( SG, output_dir, ids_to_plot=ids_closest, plot_label='deltaZ_H' )



