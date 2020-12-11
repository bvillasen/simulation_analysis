
# Parameters for changing the H and He Photoionization and Photoheating Rates

param_UVB_Rates = {}

param_UVB_Rates[0] = {}
param_UVB_Rates[0]['key'] = 'A'
param_UVB_Rates[0]['name'] = 'scale_He'
param_UVB_Rates[0]['values'] = [  0.1, 0.3, 0.5,  0.7  ]

param_UVB_Rates[1] = {}
param_UVB_Rates[1]['key'] = 'B'
param_UVB_Rates[1]['name'] = 'scale_H'
param_UVB_Rates[1]['values'] = [ 0.7, 0.8, 0.9, 1.0 ]

param_UVB_Rates[2] = {}
param_UVB_Rates[2]['key'] = 'C'
param_UVB_Rates[2]['name'] = 'deltaZ_He'
param_UVB_Rates[2]['values'] = [ -0.45, -0.15, 0.15, 0.45   ]

param_UVB_Rates[3] = {}
param_UVB_Rates[3]['key'] = 'D'
param_UVB_Rates[3]['name'] = 'deltaZ_H'
param_UVB_Rates[3]['values'] = [ -0.45, -0.15, 0.15, 0.45 ]




from simulation_parameters import grid_name, system
if system == 'Shamrock':
  from parameters_UVB_rates_scale_He  import param_UVB_Rates as param_UVB_Rates_scale_He
  from parameters_UVB_rates_scale_H   import param_UVB_Rates as param_UVB_Rates_scale_H
  from parameters_UVB_rates_deltaZ_He import param_UVB_Rates as param_UVB_Rates_deltaZ_He
  from parameters_UVB_rates_deltaZ_H  import param_UVB_Rates as param_UVB_Rates_deltaZ_H
  from parameters_UVB_rates_deltaZ_H_small  import param_UVB_Rates as param_UVB_Rates_deltaZ_H_small
  from parameters_UVB_rates_grid_16   import param_UVB_Rates as param_UVB_Rates_grid_16
  from parameters_UVB_rates_grid_256  import param_UVB_Rates as param_UVB_Rates_grid_256
  from parameters_UVB_rates_grid_81   import param_UVB_Rates as param_UVB_Rates_grid_81
  from parameters_UVB_rates_grid_36   import param_UVB_Rates as param_UVB_Rates_grid_36
  from parameters_UVB_rates_scale_H_photoion   import param_UVB_Rates as param_UVB_Rates_scale_H_photoion
  from parameters_UVB_rates_scale_H_photoheat  import param_UVB_Rates as param_UVB_Rates_scale_H_photoheat
  

  if grid_name == 'scale_He':  param_UVB_Rates = param_UVB_Rates_scale_He
  if grid_name == 'scale_H':   param_UVB_Rates = param_UVB_Rates_scale_H
  if grid_name == 'deltaZ_He': param_UVB_Rates = param_UVB_Rates_deltaZ_He
  if grid_name == 'deltaZ_H':  param_UVB_Rates = param_UVB_Rates_deltaZ_H
  if grid_name == 'deltaZ_H_small':  param_UVB_Rates = param_UVB_Rates_deltaZ_H_small
  if grid_name == 'grid_16':   param_UVB_Rates = param_UVB_Rates_grid_16
  if grid_name == 'grid_256':  param_UVB_Rates = param_UVB_Rates_grid_256
  if grid_name == 'grid_81':   param_UVB_Rates = param_UVB_Rates_grid_81
  if grid_name == 'scale_H_photoion':   param_UVB_Rates = param_UVB_Rates_scale_H_photoion
  if grid_name == 'scale_H_photoheat':   param_UVB_Rates = param_UVB_Rates_scale_H_photoheat
  if grid_name == 'grid_36':   param_UVB_Rates = param_UVB_Rates_grid_36
  
  
