import os, sys
from os import listdir
from os.path import isfile, join
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import ImageGrid
import pickle
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from tools import *
from turbo_cmap import *

data_dir = '/raid/bruno/data/'
input_dir = data_dir + 'cosmo_sims/rescaled_P19/2048_50Mpc/skewers/transmitted_flux/'


snapshots = range( 74, 170)

z_vals = []
for n_snap in snapshots:
  file_name = input_dir + f'los_transmitted_flux_{n_snap}_HI.h5'
  file = h5.File( file_name, 'r' )
  z = file.attrs['current_z']
  file.close()
  z_vals.append(z)
z_vals = np.array( z_vals )


