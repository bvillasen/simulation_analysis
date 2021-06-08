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


SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )


params = { 'scale_He':None, 'deltaZ_He':None, 'scale_H':0.86, 'deltaZ_H':0.0 }
print( f'Selecting: {params} ' )
selected_sims = SG.Select_Simulations( params, tolerance=5e-3 )
n_sims = len( selected_sims )
print( f'N Selected Sims: {n_sims} ' )



sim_id = selected_sims[0]
data_sim = SG.Grid[sim_id]
sim_key = data_sim['key']

snapshot_files = None

reduced_snaps_dir = SG.root_dir + 'reduced_snapshot_files/'
input_dir = reduced_snaps_dir + f'{sim_key}/'

if not snapshot_files:
  in_dir_files = os.listdir( input_dir)
  snaps = [ file.split('.')[0] for file in in_dir_files  ]
  boxes = [ file.split('.')[-1] for file in in_dir_files  ]
  snaps = list( set( snaps ) )
  boxes = list( set( boxes ) )
  snaps = [ int(snap) for snap in snaps ]
  boxes = [ int(snap) for snap in boxes ]
  print( f' N Snapshots: {len(snaps)} ')
  print( f' N Boxes: {len(boxes)} ' )

  box = bo

  









