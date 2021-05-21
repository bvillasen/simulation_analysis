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
from mcmc_sampling_functions import *


ps_data_dir = 'lya_statistics/data/'
output_dir = root_dir + 'interpolated_observables/'
create_directory( output_dir )



SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Load_Grid_Analysis_Data()
ps_range = SG.Get_Power_Spectrum_Range( kmax=1 )
sim_ids = SG.sim_ids

z_vals = np.array([ 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.6, 4.4,  5.0,  5.4,   ])
data_grid, data_grid_power_spectrum = Get_Data_Grid_Composite( ['P(k)', 'T0', 'tau', 'tau_HeII'], SG, z_vals=z_vals )

params = SG.parameters
scale_He_values  = params[0]['values']
scale_H_values   = params[1]['values']
deltaZ_He_values = params[2]['values']
deltaZ_H_values  = params[3]['values']


scale_He_HL  = 0.44
scale_H_HL   = 0.70
deltaZ_He_HL = 0.27
deltaZ_H_HL  = 0.05

interpolated_ps = {}
interpolated_ps['z_vals'] = z_vals
interpolated_ps['params'] = params

fields = [ 'T0', 'tau', 'tau_HeII' ]

interpolated_fields = { 'params':params }
for field in fields:
  interpolated_fields[field] = {'z':data_grid[0][field]['z']} 

vary_params = [ 'scale_He', 'scale_H', 'deltaZ_He', 'deltaZ_H' ]

p_indices = { 'scale_He':0, 'scale_H':1, 'deltaZ_He':2, 'deltaZ_H':3 }

for vary_param in vary_params:

  print( f'Vary Param: {vary_param} ' )
  interpolated_ps[vary_param] = {}
  for field in fields:
    interpolated_fields[field][vary_param] = {}

  if vary_param == 'scale_He':  vary_param_vals = scale_He_values
  if vary_param == 'scale_H':   vary_param_vals = scale_H_values
  if vary_param == 'deltaZ_He': vary_param_vals = deltaZ_He_values
  if vary_param == 'deltaZ_H':  vary_param_vals = deltaZ_H_values

  for i,param_val in enumerate(vary_param_vals):
    scale_He, scale_H, deltaZ_He, deltaZ_H = scale_He_HL, scale_H_HL, deltaZ_He_HL, deltaZ_H_HL
    if vary_param == 'scale_He':  scale_He  = param_val  
    if vary_param == 'scale_H':   scale_H   = param_val  
    if vary_param == 'deltaZ_He': deltaZ_He = param_val  
    if vary_param == 'deltaZ_H':  deltaZ_H  = param_val  
    p_vals = [ scale_He, scale_H, deltaZ_He, deltaZ_H ]
    ps_interp = Interpolate_Power_Spectrum( p_vals, data_grid_power_spectrum, SG )
    interpolated_ps[vary_param][i] = { 'ps':ps_interp, 'params':p_vals }  
  
    for field in fields:
      Interpolate_Field( p_vals, field, data_grid, SG )
      field_interp = Interpolate_Field( p_vals, field, data_grid, SG )
      p_val = p_vals[p_indices['vary_param']]
      interpolated_fields[field][vary_param][i] = {'p_val': p_val, 'field_interp':field_interp }
    

interpolated_observables = {}
interpolated_observables['power_spectrum'] = interpolated_ps
interpolated_observables['fields'] = interpolated_fields


out_file_name = output_dir + 'interpolated_observables.pkl'
Write_Pickle_Directory( interpolated_observables, out_file_name )








# 
# out_file = h5.File( out_file_name, 'w' )
# 
# observable = 'power_spectrum'
# ps_group = out_file.create_group( observable )
# ps_group.create_dataset( 'z_vals', data=z_vals )
# 
# 
# vary_param = 'scale_He'
# vary_group = ps_group.create_group( vary_param )
# ps_data = interpolated_ps[vary_param]
# for index in ps_data.keys():
#   params = ps_data[index]['params']
#   index_group = vary_group.create_group( str( )
#   index_group.atrrs['params'] = params
# 
# 
# 
# 
# 
# 
# 
# 
# 
# out_file.close()

