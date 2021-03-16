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



param_names = dst_params[0]['parameters'].keys()
n_params = len( param_names )
n_src_sims = len( src_params )
n_dst_sims = len( sdt_params )
dst_array = np.zeros([ n_dst_sims, n_params ])
src_array = np.zeros([ n_src_sims, n_params ]) 


for sim_id in range( n_src_sims ):
  sim_params = src_params[sim_id]['parameters']
  for param_id, param_name in enumerate(param_names):
    src_array[sim_id, param_id] = sim_params[param_name]
    

for sim_id in range( n_dst_sims ):
  sim_params = dst_params[sim_id]['parameters']
  for param_id, param_name in enumerate(param_names):
    dst_array[sim_id, param_id] = sim_params[param_name]
    


