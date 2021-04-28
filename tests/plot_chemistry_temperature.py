import os, sys
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from tools import *

data_dir = '/raid/bruno/data/'
input_dir = data_dir + 'cosmo_sims/chemistry_test/output_files/'
output_dir = data_dir + 'cosmo_sims/chemistry_test/figures/'
create_directory( output_dir )



n_snapshots = 170

z_all, temp_all = [], []

for n_snap in range(n_snapshots):
  
  in_file_name = input_dir + f'{n_snap}.h5'
  in_file = h5.File( in_file_name, 'r' )
  current_z = in_file.attrs['Current_z'][0]
  temp = in_file['temperature'][...]
  in_file.close()
  
  temp_max = temp.max()
  temp_min = temp.min()
  diff = ( temp_max - temp_min ) / temp_min
  
  z_all.append( current_z )
  temp_all.append( temp_max )
  
ncols, nrows = 1, 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))

ax.plot( z_all, temp_all )


ax.set_xlim( 2, 16 )


file_name = output_dir + 'temp_chemistry_test.png'
fig.savefig( file_name, bbox_inches='tight', dpi=300 )
print( f'Saved Figure: {file_name}'  )


