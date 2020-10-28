import os, sys
import h5py as h5
import numpy as np
from domain_decomposition import *


def generate_ics_particles( data_in, outDir, outputBaseName, proc_grid, box_size, grid_size ):
  
  print( '\nGenerating ICs: Particles' )
  domain = get_domain_block( proc_grid, box_size, grid_size )

  current_a = data_in['current_a']
  current_z = data_in['current_z']
  # box_size = data_in['box_size']

  data = data_in['dm']
  pos_x = data['pos_x'][...]
  pos_y = data['pos_y'][...]
  pos_z = data['pos_z'][...]
  vel_x = data['vel_x'][...]
  vel_y = data['vel_y'][...]
  vel_z = data['vel_z'][...]
  mass = data['mass'][...]
  particle_mass = mass[0]
  nPart = pos_x.shape[0]
  ids = np.arange(nPart).astype(np.int64)
  print(' Nparticles: ', nPart)

  dx = domain[0]['box']['dx']
  dy = domain[0]['box']['dy']
  dz = domain[0]['box']['dz']
  
  # print(( dx, dy, dz))

  index_x = ( pos_x / dx ).astype(np.int)
  index_y = ( pos_y / dy ).astype(np.int)
  index_z = ( pos_z / dz ).astype(np.int)
  indexs = index_x + index_y * proc_grid[0] + index_z*proc_grid[0]*proc_grid[1]


  n_local_all = []
  nprocs = proc_grid[0] * proc_grid[1] * proc_grid[2]
  for pId in range(nprocs):

    outputFileName = outDir + outputBaseName + ".{0}".format(pId)
    if nprocs == 1: outputFileName = outDir + outputBaseName 
    print(' Writing h5 file: ', outputFileName)
    outFile = h5.File( outputFileName, 'w')
    outFile.attrs['current_a'] = current_a
    outFile.attrs['current_z'] = current_z
    outFile.attrs['particle_mass'] = particle_mass

    indx = np.where(indexs == pId)[0]
    n_local = len(indx)
    n_local_all.append(n_local)
    pos_x_l = pos_x[indx]
    pos_y_l = pos_y[indx]
    pos_z_l = pos_z[indx]
    vel_x_l = vel_x[indx]
    vel_y_l = vel_y[indx]
    vel_z_l = vel_z[indx]
    mass_l = mass[indx]
    ids_l = ids[indx]
    print('  n_local: ', n_local)
    print('  Current_a: ', current_a)
    outFile.attrs['n_particles_local'] = n_local
    # outFile.attrs['N_DM_file'] = np.float(nPart)
    outFile.create_dataset( 'mass', data=mass_l )
    outFile.create_dataset( 'pos_x', data=pos_x_l.astype(np.float64) )
    outFile.create_dataset( 'pos_y', data=pos_y_l.astype(np.float64) )
    outFile.create_dataset( 'pos_z', data=pos_z_l.astype(np.float64) )
    outFile.create_dataset( 'vel_x', data=vel_x_l.astype(np.float64)  )
    outFile.create_dataset( 'vel_y', data=vel_y_l.astype(np.float64)  )
    outFile.create_dataset( 'vel_z', data=vel_z_l.astype(np.float64)  )
    outFile.create_dataset( 'particle_IDs', data=ids_l.astype(np.int64)  )

    outFile.close()
    print('')
  print("Total Particles Saved: ", sum(n_local_all))
