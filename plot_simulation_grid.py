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
matplotlib.use('Agg') 
# matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'


output_dir = root_dir + 'figures/'
create_directory( output_dir )

SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )


sim_ids = [ 0, 1, 2, 3, 4, 5 ]
scales_He = [ 1., .9, .8, .7, .6, .5 ]
SG.Load_Grid_Analysis_Data( sim_ids=sim_ids, load_fit=True )


font_size = 15

nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))

for sim_id in sim_ids:
  z = SG.Grid[sim_id]['analysis']['z']
  T0 = SG.Grid[sim_id]['analysis']['T0']
  label = r'$\beta_{HeII}$' + ' $= {0}$'.format(scales_He[sim_id])
  ax.plot( z, T0 , label=label )

ax.set_ylabel( r'$T_0$', fontsize=font_size  )
ax.set_xlabel( r'$z$', fontsize=font_size )
leg = ax.legend(loc=1, frameon=False, fontsize=font_size)
ax.set_xlim( 2, 12 )
ax.set_ylim( 3000, 18000)
figure_name = output_dir + 'grid_phase_diagram.png'
fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
print( f'Saved Figure: {figure_name}' )





nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))


for sim_id in sim_ids:
  z = SG.Grid[sim_id]['analysis']['z']
  F = SG.Grid[sim_id]['analysis']['F_mean']
  tau = - np.log( F )
  label = r'$\beta_{HeII}$' + ' $= {0}$'.format(scales_He[sim_id])
  ax.plot( z, tau, label=label )


ax.set_ylabel( r'$\tau_{eff}$', fontsize=font_size )
ax.set_xlabel( r'$z$', fontsize=font_size )
leg = ax.legend(loc=2, frameon=False, fontsize=font_size)
ax.set_xlim( 2, 6 )
ax.set_yscale('log')
ax.set_ylim( 0.1, 10 )

figure_name = output_dir + 'grid_optical_depth.png'
fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
print( f'Saved Figure: {figure_name}' )