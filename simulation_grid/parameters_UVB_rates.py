
from simulation_parameters import grid_name, system




# Select UVB parameters from file
if grid_name == '1024_P19':     from parameters_UVB_rates_P19 import param_UVB_Rates
if grid_name == '512_P19m_np3': from parameters_UVB_rates_P19m_np3 import param_UVB_Rates






# 
# if system == 'Shamrock':
#   from parameters_UVB_rates_scale_He  import param_UVB_Rates as param_UVB_Rates_scale_He
#   from parameters_UVB_rates_scale_H   import param_UVB_Rates as param_UVB_Rates_scale_H
#   from parameters_UVB_rates_deltaZ_He import param_UVB_Rates as param_UVB_Rates_deltaZ_He
#   from parameters_UVB_rates_deltaZ_H  import param_UVB_Rates as param_UVB_Rates_deltaZ_H
#   from parameters_UVB_rates_deltaZ_H_small  import param_UVB_Rates as param_UVB_Rates_deltaZ_H_small
#   from parameters_UVB_rates_grid_16   import param_UVB_Rates as param_UVB_Rates_grid_16
#   from parameters_UVB_rates_grid_256  import param_UVB_Rates as param_UVB_Rates_grid_256
#   from parameters_UVB_rates_grid_256_large  import param_UVB_Rates as param_UVB_Rates_grid_256_large
#   from parameters_UVB_rates_grid_81   import param_UVB_Rates as param_UVB_Rates_grid_81
#   from parameters_UVB_rates_grid_36   import param_UVB_Rates as param_UVB_Rates_grid_36
#   from parameters_UVB_rates_scale_H_photoion   import param_UVB_Rates as param_UVB_Rates_scale_H_photoion
#   from parameters_UVB_rates_scale_H_photoheat  import param_UVB_Rates as param_UVB_Rates_scale_H_photoheat
# 
# 
#   if grid_name == 'scale_He':  param_UVB_Rates = param_UVB_Rates_scale_He
#   if grid_name == 'scale_H':   param_UVB_Rates = param_UVB_Rates_scale_H
#   if grid_name == 'deltaZ_He': param_UVB_Rates = param_UVB_Rates_deltaZ_He
#   if grid_name == 'deltaZ_H':  param_UVB_Rates = param_UVB_Rates_deltaZ_H
#   if grid_name == 'deltaZ_H_small':  param_UVB_Rates = param_UVB_Rates_deltaZ_H_small
#   if grid_name == 'grid_16':   param_UVB_Rates = param_UVB_Rates_grid_16
#   if grid_name == 'grid_256':  param_UVB_Rates = param_UVB_Rates_grid_256
#   if grid_name == 'grid_256_large':  param_UVB_Rates = param_UVB_Rates_grid_256_large
#   if grid_name == 'grid_81':   param_UVB_Rates = param_UVB_Rates_grid_81
#   if grid_name == 'scale_H_photoion':   param_UVB_Rates = param_UVB_Rates_scale_H_photoion
#   if grid_name == 'scale_H_photoheat':   param_UVB_Rates = param_UVB_Rates_scale_H_photoheat
#   if grid_name == 'grid_36':   param_UVB_Rates = param_UVB_Rates_grid_36
# 
# 
