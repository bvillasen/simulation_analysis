import os, sys
import numpy as np
import h5py as h5
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from parameters_UVB_rates import param_UVB_Rates
from simulation_grid import Simulation_Grid
from simulation_parameters import *
from normalize_power_spectrum  import Normaliz_Flux_Power_Spectrum



use_mpi = True
if use_mpi:
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  n_procs = comm.Get_size()
else:
  rank = 0
  n_procs = 1
  

sim_ids = range(400)

sim_ids_proc = split_indices( sim_ids, rank, n_procs, adjacent=False )

print( f'Id: {rank} sim_ids:{sim_ids_proc}' )  
  

for sim_id in sim_ids_proc:  

  SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
  SG.Load_Grid_Analysis_Data( sim_ids = [sim_id], load_fit=False )
  ps_out_dir = SG.root_dir + 'flux_power_spectrum_files/'
  # if rank == 0: create_directory( ps_out_dir )

  Normaliz_Flux_Power_Spectrum( sim_id, ps_out_dir, SG )

