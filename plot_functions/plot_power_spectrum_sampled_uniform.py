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


ps_data_dir = 'lya_statistics/data/'
output_dir = root_dir + f'figures/'
create_directory( output_dir )


SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Load_Grid_Analysis_Data()
sim_ids = SG.sim_ids


z_vals = [ 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.6, 5.0, 5.4 ]
data_grid = Get_Data_Grid_Power_spectrum( z_vals, SG )

params = {}
for p_id in SG.parameters.keys():
  p_vals = SG.parameters[p_id]['values']
  params[p_id] = { 'min':min(p_vals), 'max':max(p_vals) }
  if p_id == 1: params[p_id] = { 'min':0.81, 'max':0.82 }
  params[p_id]['name'] = SG.parameters[p_id]['name'] 

print( params )
  

# Obtain distribution of the power spectrum
n_samples = 1000
ps_samples = Sample_Power_Spectrum( n_samples, params, data_grid, SG, hpi_sum=0.97, sampling='uniform'  )


Plot_Power_Spectrum_Sampling( ps_samples, ps_data_dir, output_dir, scales='small', system=system, name='fixed_beta_H' )
Plot_Power_Spectrum_Sampling( ps_samples, ps_data_dir, output_dir, scales='large', system=system, name='fixed_beta_H' )




