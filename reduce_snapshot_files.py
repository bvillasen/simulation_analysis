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

use_mpi = False
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


fields_hydro = [ 'density', 'temperature', 'HI_density', 'HeI_density', 'HeII_density' ]
fields_particles = [ 'density' ]

if type == 'hydro': fields = fields_hydro
if type == 'particles': fields = fields_particles

snaps_dir = SG.root_dir + f'snapshot_files_{type}/'
reduced_dir = SG.root_dir + f'snapshot_files_{type}_reduced/'
if rank == 0: create_directory( reduced_dir )
if use_mpi: coMM.Barrier() 

sim_ids = np.array( sim_ids )
n_sims = len( sim_ids )
indices_local = split_indices( range(n_sims), rank, n_procs  )
sims_local =  sim_ids[indices_local]
print( sims_local )

# simulation = SG.Grid[sim_id]

 



