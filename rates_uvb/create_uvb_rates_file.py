import os, sys
import numpy as np
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from tools import *
from generate_grackle_uvb_file import Load_Grackle_File, Modify_UVB_Rates, Extend_Rates_Redshift, Copy_Grakle_UVB_Rates, Modify_UVB_Rates, Write_Rates_Grackle_File



data_dir = '/gpfs/alpine/csc434/scratch/bvilasen/'
output_dir = data_dir + 'cosmo_sims/rescaled_P19/1024_50Mpc/'
create_directory( output_dir ) 


# Load the Original Rates
grackle_file_name = 'CloudyData_UVB_Puchwein2019_cloudy.h5'
grackle_data = Load_Grackle_File( grackle_file_name )
max_delta_z = 0.1
rates_data = Extend_Rates_Redshift( max_delta_z, grackle_data )
input_rates = Copy_Grakle_UVB_Rates( rates_data )

parameter_values = { 'scale_He':  0.46,
                     'scale_H':   0.78,
                     'deltaZ_He': 0.28,
                     'deltaZ_H':  0.04 }
# parameter_values = { 'scale_He':  0.46,
#                     'scale_H':   0.78,
#                     'deltaZ_He': 0.0,
#                     'deltaZ_H':  0.0 }
                    
info = 'Rates for '
for p_name in parameter_values.keys():
  p_val = parameter_values[p_name]
  info += f' {p_name}:{p_val}' 
                    
rates_modified = Modify_UVB_Rates( parameter_values, input_rates )

output_rates = {
  'UVBRates':{ 'z':rates_modified['z'], 
               'Photoheating':{ 'piHI':rates_modified['photoheating_HI'], 'piHeI':rates_modified['photoheating_HeI'], 'piHeII':rates_modified['photoheating_HeII'] }, 
               'Chemistry':{ 'k24':rates_modified['photoionization_HI'],  'k26':rates_modified['photoionization_HeI'],  'k25':rates_modified['photoionization_HeII'] },
               'info':info } }
                

out_file_name = output_dir + 'UVB_rates.h5'
Write_Rates_Grackle_File( out_file_name, output_rates )

                
# for key in ['piHI', 'piHeI', 'piHeII' ]:
#   div = output_rates['UVBRates']['Photoheating'][key] / rates_data['UVBRates']['Photoheating'][key] 
#   print( key, div )
# 
# 
# for key in ['k24', 'k26', 'k25' ]:
#   div = output_rates['UVBRates']['Chemistry'][key] / rates_data['UVBRates']['Chemistry'][key] 
#   print( key, div )