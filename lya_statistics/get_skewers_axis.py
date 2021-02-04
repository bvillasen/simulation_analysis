import sys, os, time
import numpy as np
import h5py as h5
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from tools import *
from load_data import load_snapshot_data_distributed
from spectra_functions import compute_optical_depth


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

nPoints = 2048
nx = nPoints
ny = nPoints
nz = nPoints
ncells = nx * ny * nz

uvb = 'pchw18'
# uvb = 'hm12'

# cosmo_name = 'cosmo_3'
cosmo_name = ''


inDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/output_files_{1}/'.format(nPoints, uvb, cosmo_name )
output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/skewers_{1}_HeII/'.format(nPoints, uvb, cosmo_name )
create_directory( output_dir )

# snapshots_indices_0 = [83, 86, 90, 93, 96, 99, 102, 106, 110, 114, 119, 124, 130, 136, 143, 151, 159, 169 ]
# snapshots_indices = [74, 77, 80, 83, 86, 90, 93, 96, 99, 102, 106, 110, 114, 119, 124, 130, 136, 143, 151, 159, 169]
# snapshots_indices = [74, 76, 77, 79, 80, 82, 83, 85, 86, 88, 90, 91, 93, 94, 96, 97, 99, 101, 102, 104, 106, 108, 110, 112, 114, 117, 119, 122, 124, 127, 130, 133, 136, 139, 143, 147, 151, 155, 159, 164, 169]
# snapshots_indices = list(range( 74, 170, 1))
snapshots_indices = list(range( 100, 170, 1))


data_type = 'hydro'
show_progess = True





block_size = 256
skewer_stride = 8
n_per_box = block_size // skewer_stride
n_per_dimension = nPoints // skewer_stride
n_boxes = nPoints // block_size

ids_y_local = np.arange(0, block_size, skewer_stride)
ids_z_local = np.arange(0, block_size, skewer_stride)


Lbox = 50000
proc_grid = [ 8, 8, 8]
box_size  = [ Lbox, Lbox, Lbox ]
grid_size = [ 2048, 2048, 2048 ]


for axis in ['x', 'y', 'z']:
  
  if rank >= len(snapshots_indices): continue

  nSnap = snapshots_indices[rank]

  out_file_name = output_dir + 'skewers_{0}_{1}.h5'.format(axis, nSnap)
  file_path = Path(out_file_name)
  if file_path.is_file():
    print( f' Skiping File: {out_file_name} ') 
    continue
  
  
  outFile = h5.File( out_file_name, 'w')


  skewer_ids = []
  for index_y in range(n_boxes):
    for index_z in range(n_boxes):

      if axis == 'x':
        subgrid_x = [ 0, nPoints ]
        subgrid_y = [ index_y*block_size, (index_y+1)*block_size ]
        subgrid_z = [ index_z*block_size, (index_z+1)*block_size ]

      if axis == 'y':
        subgrid_x = [ index_y*block_size, (index_y+1)*block_size ]
        subgrid_y = [ 0, nPoints ]
        subgrid_z = [ index_z*block_size, (index_z+1)*block_size ]

      if axis == 'z':
        subgrid_x = [ index_y*block_size, (index_y+1)*block_size ]
        subgrid_y = [ index_z*block_size, (index_z+1)*block_size ]
        subgrid_z = [ 0, nPoints ]


      subgrid = [ subgrid_x, subgrid_y, subgrid_z ]

      precision = np.float64

      if axis == 'x': vel_field = 'momentum_x'
      if axis == 'y': vel_field = 'momentum_y'
      if axis == 'z': vel_field = 'momentum_z'

      fields = ['density', 'temperature', vel_field, 'HI_density', 'HeII_density' ]
      # data_snapshot = load_snapshot_data_distributed( nSnap, inDir, data_type, fields, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
      data_snapshot = load_snapshot_data_distributed( data_type, fields,  nSnap, inDir,  box_size, grid_size, precision, subgrid=subgrid, proc_grid=proc_grid, show_progess=show_progess  )
      current_z    = data_snapshot['Current_z']
      density      = data_snapshot['density']
      temperature  = data_snapshot['temperature']
      HI_density   = data_snapshot['HI_density']
      HeII_density = data_snapshot['HeII_density']
      velocity     = data_snapshot[vel_field] / density 


      for id_y_local in ids_y_local:
        for id_z_local in ids_z_local:

          id_y_global = index_y * block_size + id_y_local
          id_z_global = index_z * block_size + id_z_local

          skewer_index_y = index_y * n_per_box + id_y_local/skewer_stride
          skewer_index_z = index_z * n_per_box + id_z_local/skewer_stride

          skewer_id = skewer_index_y + skewer_index_z*n_per_dimension
          skewer_ids.append(skewer_id)

          if axis == 'x':
            skewer_density      = density[:,id_y_local, id_z_local]
            skewer_temperature  = temperature[:,id_y_local, id_z_local]
            skewer_HI_density   = HI_density[:,id_y_local, id_z_local]
            skewer_HeII_density = HeII_density[:,id_y_local, id_z_local]
            skewer_velocity     = velocity[:,id_y_local, id_z_local]

          if axis == 'y':
            skewer_density      = density[id_y_local, :, id_z_local]
            skewer_temperature  = temperature[id_y_local, :, id_z_local]
            skewer_HI_density   = HI_density[id_y_local, :, id_z_local]
            skewer_HeII_density = HeII_density[id_y_local, :, id_z_local]
            skewer_velocity     = velocity[id_y_local, :, id_z_local]

          if axis == 'z':
            skewer_density      = density[id_y_local, id_z_local, :]
            skewer_temperature  = temperature[id_y_local, id_z_local, :]
            skewer_HI_density   = HI_density[id_y_local, id_z_local, :]
            skewer_HeII_density = HeII_density[id_y_local, id_z_local, :]
            skewer_velocity     = velocity[id_y_local, id_z_local, :]

          if skewer_density.shape[0]      != nPoints: print("ERROR: Skewer has the wrong length")
          if skewer_temperature.shape[0]  != nPoints: print("ERROR: Skewer has the wrong length")
          if skewer_HI_density.shape[0]   != nPoints: print("ERROR: Skewer has the wrong length")
          if skewer_HeII_density.shape[0] != nPoints: print("ERROR: Skewer has the wrong length")
          if skewer_velocity.shape[0]     != nPoints: print("ERROR: Skewer has the wrong length")


          skewer_group = outFile.create_group( str(skewer_id) )
          skewer_group.attrs['index_y'] = id_y_global
          skewer_group.attrs['index_z'] = id_z_global 
          skewer_group.create_dataset( 'density',      data=skewer_density )
          skewer_group.create_dataset( 'temperature',  data=skewer_temperature )
          skewer_group.create_dataset( 'HI_density',   data=skewer_HI_density )
          skewer_group.create_dataset( 'HeII_density', data=skewer_HeII_density )
          skewer_group.create_dataset( 'velocity',     data=skewer_velocity )

  skewer_ids = np.sort( skewer_ids )
  n_skewers = len(skewer_ids)
  test = list(range( n_skewers))
  diff = skewer_ids - test
  print('Computed {0} skewers.'.format(n_skewers))
  print(' Ids test  min:{0}  max:{1}'.format(min(diff), max(diff) ))



  outFile.attrs['current_z'] = current_z
  outFile.attrs['n'] = n_skewers
  outFile.close()
  print("Saved File: ", out_file_name) 
