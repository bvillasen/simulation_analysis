import sys, os
import numpy
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from transfer_grid_functions import *

data_dir = '/data/groups/comp-astro/bruno/cosmo_sims/sim_grid/'
src_dir = data_dir + '1024_P19m_np4_nsim256/' 
dst_dir = data_dir + '1024_P19m_np4_nsim320/'





src_params = Get_Grid_Parameter_Values( src_dir )
dst_params = Get_Grid_Parameter_Values( dst_dir )
Link_Simulation_dirctories( src_params, dst_params )



files_to_copy = ['run_output.log']

n_copied, n_skipped = 0, 0
n_dst_sims = len( dst_params )

for sim_id in range( n_dst_sims ):
  dst_sim = dst_params[sim_id]
  dst_dir = dst_sim['dir']
  if dst_sim['src_dir']:
    src_id = dst_sim['src_id']
    src_sim  = src_params[src_id]
    failed = False 
    for param_name in dst_sim['parameters']:
      if dst_sim['parameters'][param_name] != src_sim['parameters'][param_name]:
        print( "ERROR: Parameters denot match {dst_sim['parameters']}, {src_sim['parameters']}") 
        failed = True
    if not failed:
      print( f"Copying: {src_sim['parameters']} ->  {dst_sim['parameters']}  ")
      n_copied += 1
  else:
    n_skipped += 1




