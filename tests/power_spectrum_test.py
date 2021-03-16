import os, sys
import numpy as np
import matplotlib.pyplot as plt
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)

output_dir = '/home/bruno/Desktop/'

L = 1
n_points = 1000
dx = L / n_points
x_range = np.linspace( 0, L, n_points )


lambda_val = 0.1
k = 2 * np.pi / lambda_val 
signal = np.sin( k*x_range )


k_vals = 2*np.pi*np.fft.fftfreq( n_points, d=dx )
ft = 1./n_points * np.fft.fft( signal )
ft_amp2 = ft.real * ft.real + ft.imag * ft.imag

indices = k_vals > 0
k_vals = k_vals[indices]
ft_amp2 = ft_amp2[indices]

k_min = k_vals.min()
k_max = k_vals.max()
# 
# d_log_k = 0.1
# 
# k_min = np.log10( k_min )
# k_max = np.log10( k_max )
# k_start = np.log10( 0.99 * k_vals.min() )
# n_hist_edges = 1
# k_val = k_start
# while k_val < k_max:
#   n_hist_edges += 1
#   k_val += d_log_k
# hist_edges = []
# k_val = k_start
# for i in range( n_hist_edges ):
#   hist_edges.append( 10**k_val )
#   k_val += d_log_k
# intervals = np.array( hist_edges )


ncols, nrows = 1, 2
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,4*nrows))


ax = ax_l[0]
ax.plot( k_vals, ft_amp2 )











file_name = output_dir + 'ps_test.png'
fig.savefig( file_name, bbox_inches='tight', dpi=300 )





