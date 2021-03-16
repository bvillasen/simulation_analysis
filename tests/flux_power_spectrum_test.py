import os, sys
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from tools import *
from flux_power_spectrum import get_skewer_flux_power_spectrum



data_dir = '/home/bruno/Desktop/ssd_0/data/'
data_dir = '/raid/bruno/data/'
input_dir = data_dir + 'cosmo_sims/sim_grid/1024_P19m_np4_nsim256/S000_A0_B0_C0_D0/analysis_files/'
# output_dir = '/home/bruno/Desktop/'
output_dir = data_dir


n_file = 55
file_name = input_dir + f'{n_file}_analysis.h5'
infile = h5.File( file_name, 'r' )
current_z = infile.attrs['current_z'][0]



data_lya = infile['lya_statistics']
n_skewers = data_lya.attrs['n_skewers'][0]
F_mean = data_lya.attrs['Flux_mean_HI'][0] 

keys = [ 'skewers_x', 'skewers_y', 'skewers_z' ]

los_F = []
for key in keys:
  skewers_data = data_lya[key] 
  vel_Hubble = skewers_data['vel_Hubble'][...]
  los_F_axis = skewers_data['los_transmitted_flux_HI'][...]
  los_F.append( los_F_axis )

los_F = np.concatenate( los_F )

d_log_k = 0.1


n_skewers = 100
skewer_id = 0

ps_all = []
# for skewer_id in range(n_skewers): 
print_line_flush( f'Skewer {skewer_id}/{n_skewers}' )
F = los_F[skewer_id]
delta_F = F / F_mean
k_vals, flux_ps = get_skewer_flux_power_spectrum( vel_Hubble, delta_F, d_log_k=d_log_k ) 
ps_all.append( flux_ps*k_vals / np.pi )



n = len( vel_Hubble )
dv = ( vel_Hubble[-1] - vel_Hubble[0] ) / n
k_vals = 2 *np.pi * np.fft.fftfreq( n, d=dv )
ft = 1./n * np.fft.fft( delta_F )
ft_amp2 = ft.real * ft.real + ft.imag * ft.imag


signal = delta_F




v_max = vel_Hubble.max()
n_points = len( signal )
kvals =  2*np.pi*np.arange(n_points) / v_max


ft_vals = []
for k in range(n):
  A_k = 0
  for m in range(n):
    A_k += signal[m] * np.exp( -1j * kvals[k] * vel_Hubble[m]  )
  ft_vals.append( A_k )
ft_vals = np.array( ft_vals ) / n
ft_amp2 =  ft_vals.real * ft_vals.real + ft_vals.imag * ft_vals.imag   

indices = kvals > 0
kvals = kvals[indices]
ft_amp2 = ft_amp2[indices]




v_max = n*dv
n_points = len( signal )
k_vals_1 =  2*np.pi*np.arange(n_points) / v_max

vel_array = np.arange(n_points) * dv

ft_vals_1 = []
ft_vals_2 = []
for k in range(n):
  A_k_1, A_k_2 = 0, 0
  for m in range(n):
    exponent = 2*np.pi * m * k / n 
    A_k_1 += signal[m] * np.exp( -2*np.pi *1j * m * k / n  )
    k_val = k_vals_1[k] 
    vel = vel_array[m]
    exponent_1 =   k_val * vel
    # print( exponent, exponent_1 )
    A_k_2 += signal[m] * np.exp( -1j * k_val * vel  )
  ft_vals_1.append( A_k_1 )
  ft_vals_2.append( A_k_2 )
ft_vals_1 = np.array( ft_vals_1 ) / n
ft_vals_2 = np.array( ft_vals_2 ) / n
ft_amp2_vals_1 =  ft_vals_1.real * ft_vals_1.real + ft_vals_1.imag * ft_vals_1.imag   
ft_amp2_vals_2 =  ft_vals_2.real * ft_vals_2.real + ft_vals_2.imag * ft_vals_2.imag   

  
  
# ps_all = np.array( ps_all )
# ps_mean = np.mean( ps_all, axis=0 )

# diff = np.abs( ft_amp2_vals_1 - ft_amp2 ) / ft_amp2

# print( f'Diff max: {diff.max()}' )


ncols, nrows = 1, 1
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8*ncols,10*nrows))


ax = ax_l
# ax.plot( k_vals, ft_amp2 )
ax.plot( k_vals, ft_amp2_vals_1, )
ax.plot( k_vals, ft_amp2_vals_2, '--' )
ax.plot( kvals, ft_amp2, '--' )

ax.set_yscale('log')
# ax.set_xscale('log')
# 










file_name = output_dir + 'ps_test.png'
fig.savefig( file_name, bbox_inches='tight', dpi=300 )





