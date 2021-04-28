import os, sys
import numpy as np
import h5py as h5
from scipy.interpolate import interp1d


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
  
  
def Shift_UVB_Rates( delta_z, rates, param_name, interp_log=False, kind='linear' ):
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
      z_min = z_0.min()
      z_max = z_0.max()
      if interp_log: rate_0 = np.log10( rate_0 )
      rate_new = np.interp( z_0, z_new, rate_0 )
      # indx_l = np.where( z_0 == z_min )[0]
      # indx_r = np.where( z_0 == z_max )[0]
      # rate_l = rate_0[indx_l]
      # rate_r = rate_0[indx_r]
      # if kind == 'linear':
      # interp_func = interp1d( z_new, rate_0, fill_value='extrapolate' )
      # interp_func = interp1d( z_0, rate_0, bounds_error=False, fill_value=(rate_l, rate_r), kind=kind )
      # z_hr = np.linspace( z_min, z_max, 1000 )
      # rate_hr = interp_func( z_hr )
      # z_new_hr = z_hr + delta_z
      # rate_new = np.interp( z_0, z_new_hr, rate_hr )  
      # rate_new = interp_func( z_0 )
      if interp_log: rate_new = 10**rate_new
      rates[root_key][key] = rate_new


  


def Extend_Redshift( max_delta_z, z ):
  z_new = [z[0]]
  n_z = len( z )
  for i in range( n_z-1 ):
    z_l = z_new[-1]
    z_r = z[i+1]
    while( z_r - z_l ) > max_delta_z:
      z_new.append( 0.5*(z_l + z_r) )
      z_l = z_new[-1]
    z_new.append( z_r )
  z_new = np.array( z_new )
  return z_new
  
def Interpoate_Rate( z_new, z_0, rate, interp_log=False ):
  if interp_log: rate = np.log10(rate)
  rate_new = np.interp( z_new, z_0, rate )
  if interp_log: rate_new = 10**rate_new
  return rate_new 
  
def Copy_Grakle_UVB_Rates( rates_data ):
  grackle_keys = { 'Photoheating':['piHI', 'piHeI', 'piHeII'], 'Chemistry':['k24', 'k25', 'k26'] }
  output_rates = { 'UVBRates':{} }
  output_rates['UVBRates']['z'] = rates_data['UVBRates']['z'].copy()
  for grackle_key in grackle_keys:
    output_rates['UVBRates'][grackle_key] = {}
    field_keys = grackle_keys[grackle_key]
    for field_key in field_keys:
      output_rates['UVBRates'][grackle_key][field_key] = rates_data['UVBRates'][grackle_key][field_key].copy()
  return output_rates  
  
def Extend_Rates_Redshift( max_delta_z, grackle_data ):
  data_out = {}
  root_key = 'UVBRates'
  data_root = grackle_data[root_key].copy()
  z_0 =  data_root['z'][...]
  z_new = Extend_Redshift( max_delta_z, z_0)
  data_out[root_key] = {}
  data_out[root_key]['z'] = z_new
  keys_gk = ['Chemistry', 'Photoheating' ]
  for key_gk in keys_gk:
    data_gk = data_root[key_gk]
    data_out[root_key][key_gk] = {}
    for key in data_gk.keys():
      rate = data_gk[key][...]
      rate_new = Interpoate_Rate( z_new, z_0, rate )
      data_out[root_key][key_gk][key] = rate_new
  return data_out


def Modify_UVB_Rates( parameter_values, rates ):
  input_rates = rates.copy()
  rates_modified = Modify_Rates_From_Grackle_File( None, parameter_values, rates_data=input_rates )
  uvb_rates = rates_modified['UVBRates']
  z = uvb_rates['z']
  heat_HI   = uvb_rates['Photoheating']['piHI']
  heat_HeI  = uvb_rates['Photoheating']['piHeI']
  heat_HeII = uvb_rates['Photoheating']['piHeII']
  ion_HI   = uvb_rates['Chemistry']['k24']
  ion_HeI  = uvb_rates['Chemistry']['k26']
  ion_HeII = uvb_rates['Chemistry']['k25']
  rates_modified =  {}
  rates_modified['z'] = z
  rates_modified['photoheating_HI'] = heat_HI 
  rates_modified['photoheating_HeI'] = heat_HeI 
  rates_modified['photoheating_HeII'] = heat_HeII 
  rates_modified['photoionization_HI'] = ion_HI 
  rates_modified['photoionization_HeI'] = ion_HeI 
  rates_modified['photoionization_HeII'] = ion_HeII 
  return rates_modified
  

def Modify_Rates_From_Grackle_File( grackle_file_name, parameter_values, max_delta_z = 0.1, rates_data=None ):
  if not rates_data:
    grackle_data = Load_Grackle_File( grackle_file_name )
    rates = grackle_data.copy()  
    rates_data = Extend_Rates_Redshift( max_delta_z, rates )
  
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
      
      
      
def Generate_Modified_Rates_File( grackle_file_name, out_file_name, parameter_values, max_delta_z=0.1 ):      
  rates = Modify_Rates_From_Grackle_File( grackle_file_name, parameter_values, max_delta_z=max_delta_z )
  Write_Rates_Grackle_File( out_file_name, rates )
      

      







