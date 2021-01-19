import sys, os, time
import numpy as np
import h5py as h5
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from tools import *
from load_data import load_snapshot_data_distributed


use_mpi = True  

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


inDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/output_files_{1}/'.format(n_points, uvb, cosmo_name )
output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/skewers_grid_{1}/'.format(n_points, uvb, cosmo_name )
if rank == 0: create_directory( output_dir )

snaps = [ 83, 90,  96, 102,  119, 124, 130, 136, 143, 151, 159, 169, ]
snaps_boss = [  96,  102, 106, 110,  114, 119, 124, 130, 136, 143, 151, 159 ]
snapshots = list( set( snaps_boss ).union(set(snaps)))
snapshots.sort()
print(snapshots)

n_snapshot = snapshots[0]

Lbox = 50000.0     #kpc/h
n_points = 2048
proc_grid = [ 8, 8, 8]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ n_points, n_points, n_points ]


indices = range( n_points )
split_indices = split_indices( indices, rank, nprocs, adjacent=True )
index_start, index_end = split_indices[0], split_indices[-1]
print( f'Rank:{rank}   start:{index_start}   end:{index_end}' )

axis = 'x'

if axis == 'x':
  subgrid_x = [ 0, n_points ]
  subgrid_y = [ 0, n_points ]
  subgrid_z = [ index_start, index_end ]

if axis == 'y':
  subgrid_x = [ 0, npo ]
  subgrid_y = [ 0, n_points ]
  subgrid_z = [ index_z*box_size, (index_z+1)*box_size ]

if axis == 'z':
  subgrid_x = [ index_y*box_size, (index_y+1)*box_size ]
  subgrid_y = [ index_z*box_size, (index_z+1)*box_size ]
  subgrid_z = [ 0, n_points ]


subgrid = [ subgrid_x, subgrid_y, subgrid_z ]
if axis == 'x': vel_field = 'momentum_x'
if axis == 'y': vel_field = 'momentum_y'
if axis == 'z': vel_field = 'momentum_z'


data_type = 'hydro'
precision = np.float64
fields = [ 'HI_density', 'temperature', vel_field  ]
fields = [ 'HI_density' ]
data = load_snapshot_data_distributed( n_snapshot, inDir, data_type, fields, subgrid,  precision, proc_grid,  box_size, grid_size, show_progess=show_progess )
if rank == 0: print( data.keys() )


comm.Barrier()


if rank == 0: print( 'Finished Succesfully' )


