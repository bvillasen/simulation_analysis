import os, sys
from pathlib import Path
import numpy as np
import h5py as h5
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from lya_statistics_functions import Compute_Lya_Statistics_Distribution_From_Simulation
from tools import *


input_dir = '/raid/bruno/data/cosmo_sims/256_hydro_50Mpc/analysis_files/'

use_mpi = False
if use_mpi:
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  n_procs = comm.Get_size()
else:
  rank = 0
  n_procs = 1

output_dir = input_dir + 'lya_statistics/'
if rank == 0: create_directory( output_dir )
if use_mpi: comm.Barrier()


files = [f for f in listdir(input_dir) if (isfile(join(input_dir, f)) and ( f.find('_analysis') > 0) ) ]
indices = [ '{0:03}'.format( int(file.split('_')[0]) ) for file in files ]

indices.sort()
n_files = len( files )
if rank == 0: print( f' N_Analysis_Files: {n_files}' )

indices_to_generate = split_indices( indices, rank,  n_procs )
if len(indices_to_generate) == 0: exit()




for n_file in  indices_to_generate: 
  Compute_Lya_Statistics_Distribution_From_Simulation( n_file, input_dir, output_dir )





