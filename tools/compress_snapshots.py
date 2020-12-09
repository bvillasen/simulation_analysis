import os, sys
from os import listdir
from os.path import isfile, join
from tools import create_directory
import numpy as np
import time

# dataDir = '/home/bruno/Desktop/ssd_0/data/'
# dataDir = '/gpfs/alpine/proj-shared/ast149/'
dataDir = '/data/groups/comp-astro/bruno/'
input_dir = dataDir + 'cosmo_sims/256_hydro_50Mpc/output_files/'
output_dir = dataDir + 'cosmo_sims/256_hydro_50Mpc/snapshots_comparison_nyx/'

hydro = True
particles = False

def split_name( file_name, part=False):
  nSnapshot, name, nBox = file_name.split('.')
  if part:
    indx = nSnapshot.find("_particles")
    nSnapshot = nSnapshot[:indx]
  return [int(nSnapshot), int(nBox)]
  
print(( 'Input Dir: ' + input_dir ))
print(( 'Output Dir: ' + output_dir ))
create_directory( output_dir )
print("")

name_base = 'h5'

if hydro:
  dataFiles = [f for f in listdir(input_dir) if (isfile(join(input_dir, f)) and (f.find('.h5.') > 0 ) and ( f.find('_particles') < 0) ) ]
else:
  dataFiles = [f for f in listdir(input_dir) if (isfile(join(input_dir, f)) and ( f.find('_particles') > 0) ) ]
  
dataFiles = np.sort( dataFiles )
nFiles = len( dataFiles )

type = False if hydro else True
files_names = np.array([ split_name( file_name, part=type ) for file_name in dataFiles ])
snaps, boxes = files_names.T
snapshots_all = np.unique( snaps )
boxes = np.unique( boxes )
snapshots_all.sort()
nSnapshots = len( snapshots_all )
nBoxes = len( boxes )

print(( "Number of snapshots: {0}".format(nSnapshots) ))
print(( "Number of files per snapshot: {0}".format(nBoxes) ))


#Set wich snapshots to compress
snapshots_to_compress = snapshots_all
snapshots_to_compress = [0]
n_to_compress = len(snapshots_to_compress)
print(( "\nNumber of snapshots to compres: {0}".format(n_to_compress) ))
# print(( ' {0}: {1}'.format( index, snapshots_to_compress ) ))


#available Hydro Fields:
#[ density, momentum_x, momentum_y, momentum_z, Enegy, GasEnergy ]
#[ HI_density, HI_density, HeI_density, HeII_density, HeIII_density, e_density, metal_density, temperature, potential ]
# hydro_fields = 'all'
# hydro_fields = ['density' , 'momentum_x', 'HI_density', 'temperature']
# hydro_fields = ['density' ,  'HI_density', 'temperature']
# hydro_fields = ['density' ]
hydro_fields = ['density' , 'GasEnergy' ]
print(( "\nHydro fields: {0}".format(hydro_fields)))

#available Particles Fields:
#[ density, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, mass, particle_IDs ]
# particles_fields = 'all'
# particles_fields = ['density', 'vel_x', 'vel_y', 'vel_z' ]
particles_fields = []
print(( "\nParticles fields: {0}".format(particles_fields)))



Lbox = 5000    #kpc/h
n_points = 256
proc_grid = [ 2, 2, 2]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ n_points, n_points, n_points ] #Size of the simulation grid
subgrid = [ [0, n_points], [0, n_points], [0, n_points] ] #Size of the volume to load
# density = data[data_type]['density']  
# 

precision = np.float64
# # precision = np.float32
# # precision = np.float16
print(( "\nPrecision: {0}".format( precision )))

print( "\nCompressing Snapshots..." )
for nSnap in snapshots_to_compress:
  start = time.time()
  if hydro:
    data = load_snapshot_data_distributed( n_snapshot, input_dir, 'hydro', hydro_fields, subgrid,  precision, proc_grid,  box_size, grid_size, show_progess=True )


  if particles:
    pass

  end = time.time()
  print(( ' Elapsed Time: {0:.2f} min'.format((end - start)/60.) ))
