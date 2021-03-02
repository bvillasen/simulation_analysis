import os, sys
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from load_data import load_analysis_data
from spectra_functions import compute_optical_depth

#Parse Command Parameters
args = sys.argv[1:]
n_args = len(args)

use_mpi = False
if use_mpi:
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  n_procs = comm.Get_size()
else:
  rank = 0
  n_procs = 1

data_dir = '/raid/bruno/data/'
input_dir = data_dir + 'cosmo_sims/sim_grid/512_p19/analysis_files/'
output_dir = input_dir + 'figures/'
if rank == 0: create_directory( output_dir )
if use_mpi: comm.Barrier()



n_snapshot = 3
data = load_analysis_data( n_snapshot, input_dir, phase_diagram=False, lya_statistics=True, load_skewer=True )
cosmology = data['cosmology']

box = data['box'] 
skewer = data['lya_statistics']['skewer']
# skewer['HI_density'] *= 1e-5
HI_density = skewer['HI_density'] 
tau_cholla = skewer['optical_depth']
vel_Hubble_cholla = skewer['vel_Hubble'] 
temperature = skewer['temperature'] 
velocity = skewer['velocity']
F_cholla = skewer['transmitted_flux']

tau_data = compute_optical_depth( cosmology, box, skewer, space='redshift', method='error_function' )
tau = tau_data['tau']
vel_Hubble = tau_data['vel_Hubble']
F = np.exp(-tau)

diff = ( tau_cholla - tau ) / tau
print( f" Diff:  min:{diff.min()}    max:{diff.max()}")


nrows = 5
ncols = 1
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,2.5*nrows))


ax =  ax_l[0]
ax.plot( tau_data['vel_Hubble'], HI_density )
ax.set_yscale('log')

ax =  ax_l[1]
ax.plot( tau_data['vel_Hubble'],velocity )

ax =  ax_l[2]
ax.plot( tau_data['vel_Hubble'], temperature )
# ax.set_yscale('log')

ax =  ax_l[3]
ax.plot( tau_data['vel_Hubble'], tau_cholla )
ax.plot( tau_data['vel_Hubble'], tau, '--' )
ax.set_yscale('log')

ax =  ax_l[4]
ax.plot( tau_data['vel_Hubble'], F_cholla )
ax.plot( tau_data['vel_Hubble'], F, '--' )



figure_name = output_dir + 'skewer_comparison.png'
fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
print( f'Saved Figure: {figure_name}' )
