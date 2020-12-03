import os, sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
import h5py as h5
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from load_data import load_snapshot_data_distributed
from power_spectrum_functions import get_power_spectrum


data_name = 'SIMPLE_PPMP_eta035'

# data_dir = '/raid/bruno/data/'
# data_dir = '/home/bruno/Desktop/ssd_0/data/'
data_dir = '/data/groups/comp-astro/bruno/'
enzo_dir = data_dir + 'cosmo_sims/enzo/256_hydro_50Mpc/h5_files/'
cholla_dir = data_dir + f'cosmo_sims/256_hydro_50Mpc/{data_name}/'



precision = np.float64
data_type = 'hydro'
fields = [ 'density' ]

Lbox = 5000.    #kpc/h
proc_grid = [ 2, 2, 2]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 256, 256, 256 ] #Size of the simulation grid
nx, ny, nz = grid_size
dx, dy, dz = Lbox/nx, Lbox/ny, Lbox/nz
subgrid = [ [0, 256], [0, 256], [0, 256] ] #Size of the volume to load



snapshots = [ 0, 10, 29 ]

ps_data = { }

for n_snapshot in snapshots:
  data = load_snapshot_data_distributed( n_snapshot, cholla_dir, data_type, fields, subgrid,  precision, proc_grid,  box_size, grid_size, show_progess=True )
  current_z = data['Current_z']
  dens_ch = data[data_type]['density']

  file_name = enzo_dir + 'snapshot_{0:03}.h5'.format(n_snapshot)
  data_enzo = h5.File( file_name, 'r' )
  dens_en = data_enzo['gas']['density'][...]  


  ps_ch, k_vals, count = get_power_spectrum( dens_ch, Lbox, nx, ny, nz, dx, dy, dz,  n_kSamples=20)
  ps_en, k_vals, count = get_power_spectrum( dens_en, Lbox, nx, ny, nz, dx, dy, dz,  n_kSamples=20)
  diff = ( ps_ch - ps_en ) / ps_en
  print( f'n_snap: {n_snapshot}  diff: {diff}')
  ps_data[n_snapshot] = {}
  ps_data[n_snapshot]['cholla'] = ps_ch
  ps_data[n_snapshot]['enzo'] = ps_en
  ps_data[n_snapshot]['diff'] = diff
  
  

