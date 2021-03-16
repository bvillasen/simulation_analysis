import sys, os
import numpy
from shutil import copyfile, copytree
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from transfer_grid_functions import *


# data_dir = '/data/groups/comp-astro/bruno/cosmo_sims/sim_grid/'
data_dir = '/gpfs/alpine/csc434/scratch/bvilasen/cosmo_sims/sim_grid/'
src_dir = data_dir + '1024_P19m_np4_nsim256/' 
dst_dir = data_dir + '1024_P19m_np4_nsim320/'

copy_reduced_files = False
if copy_reduced_files:
  src_reduced = src_dir + 'reduced_files/'
  dst_reduced = dst_dir + 'reduced_files/'
  create_directory( dst_reduced )



src_params = Get_Grid_Parameter_Values( src_dir )
dst_params = Get_Grid_Parameter_Values( dst_dir )
Link_Simulation_dirctories( src_params, dst_params )



files_to_copy = ['run_output.log', 'param.txt', 'uvb_params.txt']



n_copied, n_skipped = 0, 0
n_dst_sims = len( dst_params )

for sim_id in range( n_dst_sims ):
  dst_sim = dst_params[sim_id]
  dst_dir = dst_sim['dir']
  src_dir = dst_sim['src_dir']
  if src_dir:
    src_id = dst_sim['src_id']
    src_sim  = src_params[src_id]
    failed = False 
    for param_name in dst_sim['parameters']:
      if dst_sim['parameters'][param_name] != src_sim['parameters'][param_name]:
        print( "ERROR: Parameters dont match {dst_sim['parameters']}, {src_sim['parameters']}") 
        failed = True
    if failed: break
    
    print( f"\nCopying: {src_sim['parameters']} ->  {dst_sim['parameters']}  ")
    src_dir_short = src_dir[src_dir.find('sim_grid')+9:]+'/'
    dst_dir_short = dst_dir[dst_dir.find('sim_grid')+9:]+'/' 
    for file in files_to_copy:
      copyfile(src_dir + '/' + file, dst_dir + '/' + file )
      print( f' Copied  {src_dir_short+file} -> {dst_dir_short+file} ' )
    
    if copy_reduced_files:
      src_red_dir = src_reduced + src_sim['name']
      dst_red_dir = dst_reduced + dst_sim['name']
      src_red_short = src_red_dir[src_red_dir.find('sim_grid')+9:]+'/'
      dst_red_short = dst_red_dir[dst_red_dir.find('sim_grid')+9:]+'/' 
      copytree(src_red_dir, dst_red_dir )
      print( f' Copied  {src_red_short} -> {dst_red_short} ' )
    
      
    
    
    n_copied += 1
  else:
    n_skipped += 1



print( f'\nSuccessfully Transfered Grid:   n_copied:{n_copied}    n_skipped:{n_skipped}')


