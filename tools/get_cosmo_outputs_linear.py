import sys, os, time
import numpy as np
import matplotlib.pyplot as plt

output_dir = '/home/bruno/Desktop/'

n_total = 299
z_start = 40.
z_end = 5.


# a_start = 1./(z_start+1)
# a_end = 1./(z_end+1)
# a_vals = np.linspace( a_start, a_end, n_total)
# np.savetxt( 'outputs_cosmo_aLin_400.txt', a_vals)


z_vals = np.linspace( z_start, z_end, n_total)
z_vals = np.concatenate( (np.array([100]), z_vals ) )
a_vals = 1./( z_vals + 1)
np.savetxt( output_dir + 'outputs_cosmo_halos.txt', a_vals)

# a_str = ''
# for a in a_vals:
#   a_str += '{0:.4f}, '.format(a)