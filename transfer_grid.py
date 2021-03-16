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

