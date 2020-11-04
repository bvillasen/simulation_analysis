import os, sys
import numpy as np
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from parameters_UVB_rates import param_UVB_Rates
from simulation_grid import Simulation_Grid
from simulation_parameters import *

create_directory( root_dir )


SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Get_Grid_Status()
SG.Load_Grid_UVB_Rates()

import matplotlib.pyplot as plt



output_dir = root_dir + 'figures/'
create_directory( output_dir )

nrows = 2
ncols = 3
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))


font_size = 15


sim_ids = SG.Grid.keys()
# sim_ids = [ 0 ]


root_keys = [ 'Chemistry', 'Photoheating' ]

for j, root_key in enumerate(root_keys):

  for sim_id in sim_ids:
    sim_rates = SG.Grid[sim_id]['UVB_rates']
    rates = sim_rates[root_key]
    
    for i,key in enumerate(rates.keys()):
      ax = ax_l[j][i]
      z = sim_rates['z'][...]
      rates_data = sim_rates[root_key][key][...]
      ax.plot( z, rates_data )
      
      
  for i,key in enumerate(rates.keys()):
    ax = ax_l[j][i]
    ax.set_ylabel( key, fontsize=font_size  )
    ax.set_xlabel( r'$z$', fontsize=font_size )  
    ax.set_yscale('log')






figure_name = output_dir + 'grid_UVB_rates.png'
fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
print( f'Saved Figure: {figure_name}' )



