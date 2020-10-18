import os, sys
import numpy as np
import h5py as h5

#Add Modules from other directories
currentDirectory = os.getcwd()
dataDirectory = currentDirectory + "/data_src/"
sys.path.extend([  dataDirectory ] )
from tools import *
from load_data import load_snapshot_data_distributed

use_mpi = True
if use_mpi:
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  n_procs = comm.Get_size()
else:
  rank = 0
  n_procs = 1



nPoints = 1024
dataDir = '/data/groups/comp-astro/bruno/'
inDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/output_files_pchw18/'.format(nPoints)
outDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/output_files_pchw18/statistics/'.format(nPoints)
if rank == 0: create_directory( outDir )




data_type = 'hydro'
fields = ['density' ]



nFiles = 170
indices = range(nFiles)
indices_to_generate = split_indices( indices, rank,  n_procs )
if len(indices_to_generate) == 0: exit()

print( f'pID:{rank}  indxs:{indices_to_generate}' )

stats = {}
for field in fields:
  stats[field] = {}
  stats[field]['min_vals'] = []
  stats[field]['max_vals'] = []
snapshots = []
for nSnap in indices_to_generate:

  snapshots.append(nSnap)
  precision = np.float32
  Lbox = 5000    #kpc/h
  if nPoints == 1024: proc_grid = [ 4, 2, 2]
  if nPoints == 2048: proc_grid = [ 8, 8, 8]
  box_size = [ Lbox, Lbox, Lbox ]
  grid_size = [ nPoints, nPoints, nPoints ] #Size of the simulation grid
  subgrid = [ [0, nPoints], [0, nPoints], [0, nPoints] ] #Size of the volume to load
  data = load_snapshot_data_distributed( nSnap, inDir, data_type, fields, subgrid,  precision, proc_grid,  box_size, grid_size, show_progess=True, get_statistics=True )
  for field in fields:
    stats[field]['min_vals'].append( data[data_type]['statistics'][field]['min'] )
    stats[field]['max_vals'].append( data[data_type]['statistics'][field]['max']  )

snapshots = np.array( snapshots )
for field in fields:
  stats[field]['min_vals'] = np.array(stats[field]['min_vals'])
  stats[field]['max_vals'] = np.array(stats[field]['max_vals']) 

if use_mpi:
  snapshots_all = comm.gather( snapshots, root=0 )
  stats_all = {}
  for field in fields:
    stats_all[field] = {}
    stats_all[field]['min_vals'] = comm.gather( stats[field]['min_vals'], root=0 )
    stats_all[field]['max_vals'] = comm.gather( stats[field]['max_vals'], root=0 )
if rank == 0: 
  snapshots_all = np.concatenate( snapshots_all)
  sort_indxs = np.argsort( snapshots_all )
  for field in fields:
    stats_all[field]['min_vals'] = np.concatenate( stats_all[field]['min_vals'] )[sort_indxs]
    stats_all[field]['max_vals'] = np.concatenate( stats_all[field]['max_vals'] )[sort_indxs]
    print(stats_all[field]['min_vals'] )
    print(stats_all[field]['max_vals'] )

# print( "nSnapshot {0}:  {1} {2}".format( nSnapshot, stats[field]['min_vals'], stats[field]['max_vals']  )
# for field in fields:
#   stats[field]['min_vals'] = np.array( stats[field]['min_vals'] )
#   stats[field]['max_vals'] = np.array( stats[field]['max_vals'] )
#   stats[field]['min_global'] = stats[field]['min_vals'].min()
#   stats[field]['max_global'] = stats[field]['max_vals'].max()
#   print( '{0}: min:{1}'.format( fields, stats[field]['min_vals']  ))
#   print( '{0}: max:{1}'.format( fields, stats[field]['max_vals']  ))
# 
# outFileName = 'stats_{0}.h5'.format(data_type)
# outFile = h5.File( outDir + outFileName, 'w' )
# for field in fields:
#   group = outFile.create_group( field )
#   group.attrs['min_global'] = stats[field]['min_global']
#   group.attrs['max_global'] = stats[field]['max_global']
#   group.create_dataset( 'min_vals', data =stats[field]['min_vals'] )
#   group.create_dataset( 'max_vals', data =stats[field]['max_vals'] )
# 
# 
# outFile.close()
