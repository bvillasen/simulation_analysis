import os, sys, time
from os import listdir
from os.path import isfile, join
import h5py as h5
import numpy as np
import subprocess
import yt
#Extend path to inclide local modules
root_dir = os.path.dirname(os.getcwd())
sub_directories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(sub_directories)
from tools import *
from load_data import get_domain_block
from ics_particles import  Compute_Particles_Domain_Indices, generate_ics_particles_distributed_single_field, Merge_Particles_Fileds
from ics_grid import generate_ics_hydro_distributed_single_field, Merge_Hydro_Fileds

type  = None
field = None
L_MPC = None
for option in sys.argv:
  if option.find("type=") != -1:  type=option[option.find('=')+1:]
  if option.find("field=") != -1: field=option[option.find('=')+1:]
  if option.find("L_MPC=") != -1: L_MPC=int(option[option.find('=')+1:])
  
if not type: 
  print( 'Set type=hydro or type=particles')
  exit(-1)  

print( f'Type: {type} ' )
print( f'Field: {field} ' )
print( f'L_MPC: {L_MPC} ' )

# Box Size
Lbox = L_MPC * 1e3

# Grid Parameters
n_points = 2048
n_boxes  = 512

# n_points = 1024
# n_boxes  = 128

# n_points = 256
# n_boxes  = 8

# L_MPC = int(Lbox/1000)
# data_dir = '/raid/bruno/data/'
# data_dir = '/data/groups/comp-astro/bruno/'
# data_dir = '/home/bruno/Desktop/ssd_0/data/'
data_dir = '/gpfs/alpine/csc434/scratch/bvilasen/'
input_dir  = data_dir + f'cosmo_sims/ics/enzo/{n_points}_{L_MPC}Mpc/'
output_dir = data_dir + f'cosmo_sims/ics/enzo/{n_points}_{L_MPC}Mpc/ics_{n_boxes}_z100/'
print(f'Input Dir: {input_dir}' )
print(f'Output Dir: {output_dir}' )
create_directory( output_dir )


# Load Enzo Dataset
nSnap = 0
snapKey = '{0:03}'.format(nSnap)
inFileName = 'DD0{0}/data0{0}'.format( snapKey)
ds = yt.load( input_dir + inFileName )
data = ds.all_data()
h = ds.hubble_constant
current_z = np.float(ds.current_redshift)
current_a = 1./(current_z + 1)

data_enzo = { 'dm':{}, 'gas':{} }
data_enzo['current_a'] = current_a
data_enzo['current_z'] = current_z

# Domain decomposition
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ n_points, n_points, n_points ]
if n_boxes == 1:   proc_grid = [ 1, 1, 1 ]
if n_boxes == 2:   proc_grid = [ 2, 1, 1 ]
if n_boxes == 8:   proc_grid = [ 2, 2, 2 ]
if n_boxes == 16:  proc_grid = [ 4, 2, 2 ]
if n_boxes == 64:  proc_grid = [ 4, 4, 4 ]
if n_boxes == 128: proc_grid = [ 8, 4, 4 ]
if n_boxes == 512: proc_grid = [ 8, 8, 8 ]

print( 'Output Domain:')
print( f' Box  Size: {box_size}')
print( f' Grid Size: {grid_size}')
print( f' Proc Grid: {proc_grid}')
time.sleep(2)

if type == 'particles':

  Get_Free_Memory( print_out=True )
  particles_domain_indices = Compute_Particles_Domain_Indices( box_size, grid_size, proc_grid, data, ds, output_dir, type_int=np.int16 )
  Get_Free_Memory( print_out=True )


  field_list = [ 'mass', 'pos_x', 'pos_y', 'pos_z', 'vel_x', 'vel_y', 'vel_z' ]
  if field in field_list:
    generate_ics_particles_distributed_single_field( field, particles_domain_indices, proc_grid, grid_size, output_dir, ds, data  )
    Get_Free_Memory( print_out=True )
    
  if field == 'all':
    for field in field_list:
      generate_ics_particles_distributed_single_field( field, particles_domain_indices, proc_grid, grid_size, output_dir, ds, data  )
      Get_Free_Memory( print_out=True )
      time.sleep(2)

  if field == 'merge': Merge_Particles_Fileds( field_list, proc_grid, grid_size, output_dir, output_base_name = 'particles.h5', n_snapshot=0 )


if type == 'hydro':
  
  data_grid = ds.covering_grid( level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions )
  

  field_list = [ 'density', 'velocity_x', 'velocity_y', 'velocity_z', 'thermal_energy' ]

  if field in field_list:
    Get_Free_Memory( print_out=True )
    generate_ics_hydro_distributed_single_field( field, proc_grid, output_dir, ds, data_grid )

    
  if field == 'all':
    for field in field_list:
      Get_Free_Memory( print_out=True )
      generate_ics_hydro_distributed_single_field( field, proc_grid, output_dir, ds, data_grid )
      print('')
      time.sleep(2)
   
  if field == 'merge': Merge_Hydro_Fileds( field_list, proc_grid, output_dir, output_base_name='h5', n_snapshot=0 ) 

    
    
    
    
    
    
    