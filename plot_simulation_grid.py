import os, sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from parameters_UVB_rates import param_UVB_Rates
from simulation_grid import Simulation_Grid
from simulation_parameters import *

import matplotlib
# matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'


output_dir = root_dir + 'figures/'
create_directory( output_dir )

SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )


sim_ids = [ 0, 1 ]
SG.Load_Grid_Analysis_Data( sim_ids=sim_ids, load_fit=True )



nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,2.5*nrows))

for sim_id in sim_ids:
  z = SG.Grid[sim_id]['analysis']['z']
  T0 = SG.Grid[sim_id]['analysis']['T0']
  ax.plot( z, T0 )



figure_name = output_dir + 'grid_phase_diagram.png'
fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
print( f'Saved Figure: {figure_name}' )





nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,2.5*nrows))


for sim_id in sim_ids:
  z = SG.Grid[sim_id]['analysis']['z']
  F = SG.Grid[sim_id]['analysis']['F_mean']
  tau = - np.log( F )
  ax.plot( z, tau )



figure_name = output_dir + 'grid_optical_depth.png'
fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
print( f'Saved Figure: {figure_name}' )