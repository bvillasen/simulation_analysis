import os, sys
import numpy as np
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from parameters_UVB_rates import param_UVB_Rates
from simulation_grid import Simulation_Grid
from simulation_parameters import *
from plot_UVB_Rates import Plot_Grid_UVB_Rates

use_mpi = True
if use_mpi:
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  n_procs = comm.Get_size()
else:
  rank = 0
  n_procs = 1
  
  

SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
sim_ids = SG.sim_ids

precision = np.float32

type = 'hydro'


fields_hydro = [ 'density', 'temperature', 'HI_density', 'HeII_density' ]
fields_particles = [ 'density' ]

if type == 'hydro': fields_list = fields_hydro
if type == 'particles': fields_list = fields_particles

snaps_dir = SG.root_dir + f'snapshot_files_{type}/'
reduced_dir = SG.root_dir + f'reduced_snapshot_files/'
if rank == 0: create_directory( reduced_dir )
if use_mpi: comm.Barrier() 

sim_ids = np.array( list(sim_ids) )
n_sims = len( sim_ids )
indices_local = split_indices( range(n_sims), rank, n_procs  )
sims_local =  sim_ids[indices_local]
print( sims_local )


n_sims_local = len( sims_local )
file_counter = 0
time_start = time.time()

# sim_id = sims_local[0]
for sim_id in sims_local:
  simulation = SG.Grid[sim_id]
  sim_key = simulation['key']

  input_dir  = snaps_dir + f'{sim_key}/'
  output_dir = reduced_dir + f'{sim_key}/'
  create_directory( output_dir )


  files = listdir( input_dir )
  n_files_per_sim = len( files )
  n_files_local = n_sims_local * n_files_per_sim

  if file_counter == 0:
    if rank == 0: print(f'N files per snapshot: {n_files_per_sim}')
    if rank == 0: print( f'Splitting over {n_procs} processes ' )
    
    
  for file_name in files:
    in_file  = h5.File( input_dir + file_name, 'r' )
    out_file = h5.File( output_dir + file_name, 'w' )
    
    # Copy the header
    for key in in_file.attrs.keys():
      out_file.attrs[key] = in_file.attrs[key]
      
    # Copy the fields
    for field in fields_list:
      data = in_file[field][...].astype( precision )
      out_file.create_dataset( field, data=data )

    # Close
    in_file.close()
    out_file.close() 
    file_counter += 1
    if rank == 0: print_progress( file_counter, n_files_local, time_start )


if rank == 0: 
  print( '\nFinised Successfully')
