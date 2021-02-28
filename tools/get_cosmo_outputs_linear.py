import sys, os, time
import numpy as np

# output_dir = '/home/bruno/Desktop/'
output_dir = '/Users/bruno/Desktop/'


z_range_all = [ [ 2, 2.4 ], [ 2.5, 4.4 ], [ 4.5, 6.4 ], [6.5, 12] ]
n_zvals_all = [     5,          39,             39,         11    ]
n_ranges = len( z_range_all )
z_vals_all = np.array([])

for i in range( n_ranges ):
  z_range = z_range_all[i]
  n_zvals = n_zvals_all[i]

  z_start, z_end = z_range
  z_vals = np.linspace( z_start, z_end, n_zvals )
  z_vals_all = np.concatenate([ z_vals_all, z_vals ])

n_z_vals = len( z_vals_all )
print(f'n z vals: {n_z_vals} ' )
print(z_vals_all )


z_vals = z_vals_all[::-1] 
a_vals = 1 / ( z_vals + 1 )

file_name = output_dir + 'outputs_cosmo_94.txt' 
np.savetxt( file_name, a_vals )