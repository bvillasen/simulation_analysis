import os, sys
import numpy as np
import h5py as h5



def Load_Grackle_File( grackle_file_name ):
  grackle_file = h5.File( grackle_file_name, 'r' )

  data_out = {}
  root_key = 'UVBRates'
  data_root = grackle_file[root_key]
  data_out[root_key] = {}
  data_out[root_key]['z'] = data_root['z'][...]
  keys_gk = ['Chemistry', 'Photoheating' ]
  for key_gk in keys_gk:
    data_gk = data_root[key_gk]
    data_out[root_key][key_gk] = {}
    for key in data_gk.keys():
      data_out[root_key][key_gk][key] = data_gk[key][...]
  return data_out
  
  
def Shift_UVB_Rates( delta_z, rates, param_name, interp_log=False ):
  keys_He = { 'Chemistry':['k25'], 'Photoheating':['piHeII']}
  keys_H  = { 'Chemistry':['k24', 'k26'], 'Photoheating':['piHI', 'piHeI']}
  if param_name == 'shift_H': keys = keys_H
  elif param_name == 'shift_He': keys = keys_He
  else:
    print( "ERROR: Wrong parameter name for UVB rates shift")
    exit(-1)
  z_0 = rates['z'].copy()
  z_new = z_0 + delta_z
  root_keys = [ 'Chemistry', 'Photoheating' ]
  for root_key in root_keys:
    for key in keys[root_key]:
      rate_0 = rates[root_key][key].copy()
      if interp_log: rate_0 = np.log10( rate_0 )
      rate_new = np.interp( z_0, z_new, rate_0 )
      if interp_log: rate_new = 10**rate_new
      rates[root_key][key] = rate_new


  


def Modify_Rates_From_Grackle_File( grackle_file_name, parameter_values ):
  grackle_data = Load_Grackle_File( grackle_file_name )
  rates_data = grackle_data.copy()
  
  info = 'Rates for '
  for p_name in parameter_values.keys():
    p_val = parameter_values[p_name]
    info += f' {p_name}:{p_val}' 
  
    if p_name == 'scale_H':
      rates_data['UVBRates']['Chemistry']['k24'] *= p_val
      rates_data['UVBRates']['Chemistry']['k26'] *= p_val
      rates_data['UVBRates']['Photoheating']['piHI'] *= p_val
      rates_data['UVBRates']['Photoheating']['piHeI'] *= p_val
  
    if p_name == 'scale_He':
      rates_data['UVBRates']['Chemistry']['k25'] *= p_val
      rates_data['UVBRates']['Photoheating']['piHeII'] *= p_val

    if p_name == 'scale_H_photoion':
      rates_data['UVBRates']['Chemistry']['k24'] *= p_val
      rates_data['UVBRates']['Chemistry']['k26'] *= p_val
    
    if p_name == 'scale_H_photoheat':
      rates_data['UVBRates']['Photoheating']['piHI'] *= p_val
      rates_data['UVBRates']['Photoheating']['piHeI'] *= p_val
  
    if p_name == 'scale_He_photoion':  rates_data['UVBRates']['Chemistry']['k25'] *= p_val
    if p_name == 'scale_He_photoheat': rates_data['UVBRates']['Photoheating']['piHeII'] *= p_val

      
 
    if p_name  == 'deltaZ_H': Shift_UVB_Rates( p_val, rates_data['UVBRates'], 'shift_H' )
    if p_name  == 'deltaZ_He': Shift_UVB_Rates( p_val, rates_data['UVBRates'], 'shift_He' )

  rates_data['UVBRates']['info'] = info
  return rates_data


def Write_Rates_Grackle_File( out_file_name, rates ):
  root_key = 'UVBRates'
  info = rates[root_key]['info']
  print( f'  Writing {info}' )
  
  len_info = len(info)
  type_info = f'|S{len_info}'
  info = np.array(info, dtype=type_info )
  out_file = h5.File( out_file_name, 'w' )
  root_group = out_file.create_group( root_key )
  root_group.create_dataset( 'Info', data=info )
  root_group.create_dataset( 'z', data=rates[root_key]['z'] )
  data_keys = [ 'Chemistry', 'Photoheating' ]
  for data_key in data_keys:
    group_data = rates[root_key][data_key] 
    data_group = root_group.create_group( data_key )
    for key in group_data.keys():
      data_group.create_dataset( key, data=group_data[key] )
  out_file.close()
  print( f' Saved File: {out_file_name}')
      
      
      
def Generate_Modified_Rates_File( grackle_file_name, out_file_name, parameter_values ):      
  rates = Modify_Rates_From_Grackle_File( grackle_file_name, parameter_values )
  Write_Rates_Grackle_File( out_file_name, rates )
      

      







