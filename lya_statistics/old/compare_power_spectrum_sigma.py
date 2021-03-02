import os, sys
import numpy as np
import h5py as h5
import pickle
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.transforms as tfrms
import palettable
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from tools import *
from turbo_cmap import *

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'


uvb = 'pchw18'
# uvb = 'hm12'
dataDir = '/home/bruno/Desktop/data/'
# dataDir = '/raid/bruno/data/'
# dataDir = '/data/groups/comp-astro/bruno/'
simulation_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/'
input_dir = simulation_dir + 'transmited_flux_{0}_review/bootstraped_power_spectrum/statistics/'.format(uvb)
output_dir = simulation_dir + '/figures/flux_power_spectrum/'
create_directory( output_dir )

snaps = [ 83, 90,  96, 102,  119, 124, 130, 136, 143, 151, 159, 169, ]
snaps_boss = [  96,  102, 106, 110,  114, 119, 124, 130, 136, 143, 151, 159 ]
snapshots = list( set( snaps_boss ).union(set(snaps)))
snapshots.sort()
print(snapshots)

n_in_sample_list = [ 100, 500, 1000, 2500, 5000, 10000, 25000, 50000, 60000 ]
# n_in_sample_list = [  500,  2500, 10000, 50000,  ]
# for n_in_sample in n_in_sample_list:
snaps_to_plot = [ 169 ]
data_to_plot = {}



for index,n_snap in enumerate(snaps_to_plot):
  data_to_plot[index] = {}
  
  print( f'Plotting Snapshot: {n_snap}')

  stats_file = input_dir + f'stats_{n_snap}.pkl'
  print( f'Loading File: {stats_file}')
  data = pickle.load( open( stats_file, 'rb' ) ) 
  current_z = data['current_z']
  n_iterations = data['n_iterations']
  ps_mean = data['mean']
  k_vals = data['k_vals']
  data_to_plot[index]['mean'] = ps_mean
  data_to_plot[index]['current_z'] = current_z
  data_to_plot[index]['k_vals'] = k_vals

  bootstrap = data['bootstrap']
  
  for n_in_sample in n_in_sample_list:
    data_to_plot[index][n_in_sample] = {}
    stats      = bootstrap[n_in_sample]['statistics']
    data_to_plot[index][n_in_sample]['sigma']  = stats['sigma'] 
    

index = 0
sigma_vals = {}    
for n_in_sample in n_in_sample_list:
  sigma_vals[n_in_sample] = data_to_plot[index][n_in_sample]['sigma']    
