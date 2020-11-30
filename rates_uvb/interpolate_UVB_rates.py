import os, sys
import numpy as np
import h5py as h5
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


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
  
def Interpoate_Rate( z_new, z_0, rate, interp_log=True ):
  if interp_log: rate = np.log10(rate)
  rate_new = np.interp( z_new, z_0, rate )
  if interp_log: rate_new = 10**rate_new
  return rate_new 
  
  
  
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


  

grackle_file_name =  'CloudyData_UVB_Puchwein2019_cloudy.h5' 
grackle_data = Load_Grackle_File( grackle_file_name )  


max_delta_z = 0.1
data_extend = Extend_Rates_Redshift( max_delta_z, grackle_data )

# 
# z_0 = z.copy()
# rate_0 = rate.copy()
# z_min = z_0.min()
# z_max = z_0.max()
# if interp_log: rate_0 = np.log10( rate_0 )
# indx_l = np.where( z_0 == z_min )[0]
# indx_r = np.where( z_0 == z_max )[0]
# rate_l = rate_0[indx_l]
# rate_r = rate_0[indx_r]
# # if kind == 'linear':
# # rate_new = np.interp( z_0, z_new, rate_0 )
# # interp_func = interp1d( z_new, rate_0, fill_value='extrapolate' )
# interp_func = interp1d( z_0, rate_0, bounds_error=False, fill_value=(rate_l, rate_r), kind=kind )
# z_hr = np.linspace( z_min, z_max, 10000 )
# rate_hr = interp_func( z_hr )
# z_new_hr = z_hr + delta_z
# rate_new = np.interp( z_0, z_new_hr, rate_hr )  
# 
# if interp_log: 
#   rate_hr = 10**rate_hr
#   rate_new = 10**rate_new
# 
# 
# nrows = 1
# ncols = 1
# fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))
# 
# 
# ax.plot( z, rate )
# ax.plot( z_hr, rate_hr, '--' )
# ax.plot( z_new_hr, rate_hr, '--' )
# ax.plot( z_0, rate_new, '--' )
# 
# ax.set_yscale('log')
# 
# 
# output_dir = '/home/bruno/Desktop/'
# figure_name = output_dir + 'grid_UVB_rates.png'
# fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
# print( f'Saved Figure: {figure_name}' )
