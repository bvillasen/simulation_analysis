import sys, os, time
import numpy as np

output_dir = '/home/bruno/Desktop/'
# output_dir = '/Users/bruno/Desktop/'


a_vals_2048 = np.loadtxt( output_dir + 'outputs_cosmo_2048.txt' ) 
delta_a_vals_0 = a_vals_2048[1:] - a_vals_2048[:-1]
delta_a_0 = delta_a_vals_0[-1]


# 
a_start, a_end = a_vals_2048[-1], 0.5
a_vals_0 = []
a = a_start
while a < a_end:
  a += delta_a_0
  a_vals_0.append( a )

a_vals_0 = np.array( a_vals_0 )[:-2]



z_range_all = [ [ 0,  1 ],  ]
n_zvals_all = [     93,     ]
n_ranges = len( z_range_all )
z_vals_all = np.array([])

for i in range( n_ranges ):
  z_range = z_range_all[i]
  n_zvals = n_zvals_all[i]

  z_start, z_end = z_range
  z_vals = np.linspace( z_start, z_end, n_zvals )
  if i > 0: z_vals = z_vals[1:]
  z_vals_all = np.concatenate([ z_vals_all, z_vals ])

n_z_vals = len( z_vals_all )
print(f'n z vals: {n_z_vals} ' )
print(z_vals_all )


z_vals = z_vals_all[::-1] 
a_vals = 1 / ( z_vals + 1 )

a_vals = np.concatenate( [ a_vals_0, a_vals ])

delta_a = a_vals[1:] - a_vals[:-1]

file_name = output_dir + f'outputs_cosmo_2048_z0.txt' 
np.savetxt( file_name, a_vals )


# 
# z_range_all = [ [ 2, 2.4 ], [ 2.5, 4.4 ], [ 4.5, 6.4 ], [6.5, 12] ]
# n_zvals_all = [     5,          39,             39,         11    ]
# n_zvals_all = [     5,          20,             20,         11    ]
# n_ranges = len( z_range_all )
# z_vals_all = np.array([])
# 
# for i in range( n_ranges ):
#   z_range = z_range_all[i]
#   n_zvals = n_zvals_all[i]
# 
#   z_start, z_end = z_range
#   z_vals = np.linspace( z_start, z_end, n_zvals )
#   z_vals_all = np.concatenate([ z_vals_all, z_vals ])
# 
# n_z_vals = len( z_vals_all )
# print(f'n z vals: {n_z_vals} ' )
# print(z_vals_all )
# 
# 
# z_vals = z_vals_all[::-1] 
# a_vals = 1 / ( z_vals + 1 )
# 
# file_name = output_dir + f'outputs_cosmo_analysis_{n_z_vals}.txt' 
# np.savetxt( file_name, a_vals )