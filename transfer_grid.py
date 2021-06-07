import sys, os, time
import numpy
from shutil import copyfile, copytree, move
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from transfer_grid_functions import *


data_dir = '/data/groups/comp-astro/bruno/cosmo_sims/sim_grid/'
# data_dir = '/gpfs/alpine/csc434/scratch/bvilasen/cosmo_sims/sim_grid/'
# data_dir = '/raid/bruno/data/cosmo_sims/sim_grid/'


# src_grid_dir = data_dir + '1024_P19m_np4_nsim256/' 
# dst_grid_dir = data_dir + '1024_P19m_np4_nsim320/'
 
src_grid_dir = data_dir + '1024_P19m_np4_nsim320/' 
dst_grid_dir = data_dir + '1024_P19m_np4_nsim400/'


copy_reduced_files = False

if copy_reduced_files:
  src_reduced = src_grid_dir + 'reduced_files/'
  dst_reduced = dst_grid_dir + 'reduced_files/'
  create_directory( dst_reduced )

copy_power_spectrum_files = False
if copy_power_spectrum_files:
  src_ps = src_grid_dir + 'flux_power_spectrum_files/'
  dst_ps = dst_grid_dir + 'flux_power_spectrum_files/'
  create_directory( dst_ps )



move_reduced_snapshots = True
if move_reduced_snapshots:
  src_reduced_snapshots = src_grid_dir + 'reduced_snapshot_files/'
  dst_reduced_snapshots = dst_grid_dir + 'reduced_snapshot_files/'
  create_directory( dst_reduced_snapshots )



  # params_default = {'deltaZ_H':0.0, 'deltaZ_He':0.2 }
src_params = Get_Grid_Parameter_Values( src_grid_dir )
dst_params = Get_Grid_Parameter_Values( dst_grid_dir )
Link_Simulation_dirctories( src_params, dst_params )




files_to_copy = []

# files_to_copy = ['run_output.log', 'param.txt', 'uvb_params.txt']
# directories_to_copy = [ 'analysis_files' ]

files_to_copy = []
directories_to_copy = [ ]

copy_simulation_directory = False



n_copied, n_skipped = 0, 0
n_dst_sims = len( dst_params )
dst_ids_to_transfer = range( n_dst_sims )
# dst_ids_to_transfer = [ 64 ]

for sim_id in dst_ids_to_transfer:
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
  
  
    for dir in directories_to_copy:
      dst_indir = dst_dir + '/' + dir 
      dst_dir_content = os.listdir(dst_indir)
      if len(dst_dir_content) == 0:
        print( f' Deleting Empty: {dst_indir}')
        os.rmdir( dst_dir + '/' + dir )
        print( f' Copying  {src_dir_short+dir} -> {dst_dir_short+dir} ' )  
        copytree(src_dir + '/' + dir, dst_indir )
        print( f' Copied   {src_dir_short+dir} -> {dst_dir_short+dir} ' )  
  

    if move_reduced_snapshots:
      src_red_dir = src_reduced_snapshots + src_sim['name']
      dst_red_dir = dst_reduced_snapshots + dst_sim['name']
      src_red_short = src_red_dir[src_red_dir.find('sim_grid')+9:]+'/'
      dst_red_short = dst_red_dir[dst_red_dir.find('sim_grid')+9:]+'/'
      print( f' Moving  {src_red_short} -> {dst_red_short} ' ) 
      src_exits = os.path.isdir( src_red_dir)
      if not src_exits: 
        print( f'ERROR: Directory Doesnt Exists: { src_red_dir} ')
        break
      dst_exits = os.path.isdir( dst_red_dir)
      if dst_exits: 
        print( f'ERROR: Directory Exists: {dst_red_dir} ')
        continue
      files = os.listdir( src_red_dir )
      n_files = len(files)
      print( f' N Files in src dir: {n_files} ' )
      if n_files != 1920:
        print( 'ERROR: n_files != 1920 ')
        break
      # move( src_red_dir, dst_red_dir )
      # # print ( f'{src_red_dir} -> {dst_red_dir} ')
      # files = os.listdir( dst_red_dir )
      # n_files = len(files)
      # print( f' N Files in dst dir: {n_files} ' )
      # if n_files != 1920:
      #   print( 'ERROR: n_files != 1920 ')
      #   break
  
    if copy_reduced_files:
      src_red_dir = src_reduced + src_sim['name']
      dst_red_dir = dst_reduced + dst_sim['name']
      src_red_short = src_red_dir[src_red_dir.find('sim_grid')+9:]+'/'
      dst_red_short = dst_red_dir[dst_red_dir.find('sim_grid')+9:]+'/' 
      dst_dir_content = os.listdir(dst_red_dir)
      # print(dst_dir_content )
      dst_analysis = dst_red_dir + '/analysis_files'
      dst_mcmc = dst_red_dir + "/analysis_files/fit_mcmc"
      copy_directory = False
      empty_mcmc = False
      if os.path.isdir( dst_analysis ):
        if os.path.isdir( dst_mcmc ):
          if len( os.listdir(dst_mcmc) ) == 0:
            empty_mcmc = True
            print( f' Deleting Empty: {dst_red_dir + "/analysis_files/fit_mcmc"}')
            os.rmdir( dst_red_dir + '/analysis_files/fit_mcmc' )
            copy_directory = True
        if len( os.listdir(dst_analysis) ) == 0:
          print( f' Deleting Empty: {dst_red_dir + "/analysis_files"}')
          os.rmdir( dst_red_dir + '/analysis_files' )
          copy_directory = True
  # 
      if copy_directory:
        print( f' Deleting Empty: { dst_red_dir }')
        os.rmdir( dst_red_dir )
        copytree(src_red_dir, dst_red_dir )
        print( f' Copied  {src_red_short} -> {dst_red_short} ' )
      else: print( f'Skipped: {dst_red_dir}' )
      # time.sleep(0.1)
  
    if copy_power_spectrum_files:
      src_ps_dir = src_ps + src_sim['name']
      dst_ps_dir = dst_ps + dst_sim['name']
      src_ps_short = src_ps_dir[src_ps_dir.find('sim_grid')+9:]+'/'
      dst_ps_short = dst_ps_dir[dst_ps_dir.find('sim_grid')+9:]+'/' 
      dst_dir_content = os.listdir(dst_ps_dir)
      if len(dst_dir_content) == 0:
        os.rmdir( dst_ps_dir )
        copytree(src_ps_dir, dst_ps_dir )
        print( f' Copied  {src_ps_short} -> {dst_ps_short} ' )
        time.sleep( 0.1 )
  
  
  
    n_copied += 1
  else:
    n_skipped += 1



print( f'\nSuccessfully Transfered Grid:   n_copied:{n_copied}    n_skipped:{n_skipped}')


