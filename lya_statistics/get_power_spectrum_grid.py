import sys, os, time
import numpy as np
import h5py as h5
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from tools import *
from load_data import load_snapshot_data_distributed
from spectra_functions import compute_optical_depth

use_mpi = False

if use_mpi :
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  nprocs = comm.Get_size()
else:
  rank = 0
  nprocs = 1

dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/gpfs/alpine/proj-shared/ast149/'

n_points = 2048
nx = n_points
ny = n_points
nz = n_points
ncells = nx * ny * nz

uvb = 'pchw18'
# uvb = 'hm12'

# cosmo_name = 'cosmo_3'
cosmo_name = ''

show_progess = False
if rank == 0: show_progess = True


input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/skewers_grid_{1}/'.format(n_points, uvb, cosmo_name )


snaps = [ 83, 90,  96, 102,  119, 124, 130, 136, 143, 151, 159, 169, ]
snaps_boss = [  96,  102, 106, 110,  114, 119, 124, 130, 136, 143, 151, 159 ]
snapshots = list( set( snaps_boss ).union(set(snaps)))
snapshots.sort()
# print(snapshots)

n_snapshot = snapshots[-2]
snapshot_dir = input_dir + 'snapshot_{0:03}/'.format( n_snapshot ) 
if rank == 0: create_directory( snapshot_dir )

axis = 'x'
    
file_name = snapshot_dir + f'skewers_subgrid_{axis}_{rank}.h5'
file = h5.File( file_name, 'r' )
current_z     = file.attrs['current_z']
subgrid_shape = file.attrs['subgrid_shape']

vel_hubble = file['vel_hubble'][...]
F_subgrid = file[F_subgrid][...]
file.close()

    
    







if rank == 0: print( 'Finished Succesfully' )


