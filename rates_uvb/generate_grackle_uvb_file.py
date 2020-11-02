import os, sys
import numpy as np
import h5py as h5



def Load_Grackle_File( grackle_file_name ):
  grackle_file = h5.File( grackle_file_name, 'r' )

  data_out = {}

  # root_key = 'CoolingRates'
  # data_root = grackle_file[root_key]
  # data_out[root_key] = {}
  # 
  # key_gk = 'Primordial'
  # data_gk = data_root[key_gk]
  # data_out[root_key][key_gk] = {}

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
  
  


def Modify_Rates_From_Grackle_File( grackle_file_name, scale_HI, scale_HeII, deltaZ_HI, deltaZ_HeII ):
  grackle_data = Load_Grackle_File( grackle_file_name )
  rates_data = grackle_data.copy()
  
  # rates_data['UVBRates']['z'] = grackle_data['UVBRates']['z']
  
  info = f'Rates for scale_HI:{scale_HI} scale_HeII:{scale_HeII} deltaZ_HI:{deltaZ_HI} deltaZ_HeII:{deltaZ_HeII}'
  rates_data['UVBRates']['info'] = info
  
  
  rates_data['UVBRates']['Chemistry']['k24'] *= scale_HI
  rates_data['UVBRates']['Chemistry']['k25'] *= scale_HI
  rates_data['UVBRates']['Chemistry']['k26'] *= scale_HeII

  rates_data['UVBRates']['Photoheating']['piHI']   *= scale_HI
  rates_data['UVBRates']['Photoheating']['piHeI']  *= scale_HI
  rates_data['UVBRates']['Photoheating']['piHeII'] *= scale_HeII

  return rates_data


def Write_Rates_Grackle_File( out_file_name, rates ):
  root_key = 'UVBRates'
  info = rates[root_key]['info']
  print( f'  Writing {info}' )
  
  out_file = h5.File( out_file_name, 'w' )
  root_group = out_file.create_group( root_key )
  root_group.create_dataset( 'Info', data=rates[root_key]['info'] )
  root_group.create_dataset( 'z', data=rates[root_key]['z'] )
  data_keys = [ 'Chemistry', 'Photoheating' ]
  for data_key in data_keys:
    group_data = rates[root_key][data_key] 
    data_group = root_group.create_group( data_key )
    for key in group_data.keys():
      data_group.create_dataset( key, data=group_data[key] )
  out_file.close()
  print( f' Saved File: {out_file_name}')
      
      
      
def Generate_Modified_Rates_File( grackle_file_name, out_file_name, scale_HI, scale_HeII, deltaZ_HI, deltaZ_HeII  ):      
  rates = Modify_Rates_From_Grackle_File( grackle_file_name, scale_HI, scale_HeII, deltaZ_HI, deltaZ_HeII )
  Write_Rates_Grackle_File( out_file_name, rates )
      

      







