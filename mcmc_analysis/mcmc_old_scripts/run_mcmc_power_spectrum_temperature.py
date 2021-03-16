import os, sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from parameters_UVB_rates import param_UVB_Rates
from simulation_grid import Simulation_Grid
from simulation_parameters import *
from mcmc_functions import *
from mcmc_data_functions import *
from data_thermal_history import *
from mcmc_plotting_functions import *

# data_sets = [ 'Boss', 'Walther', 'Boera', 'Viel' ]
data_sets = [ 'Boss' ]
# data_sets = [ 'Walther' ]
# data_sets = [ 'Boera' ]
# data_sets = [ 'Boss', 'Walther' ]
# data_sets = [ 'Walther', 'Boera' ]
# data_sets = [ 'Walther', 'Viel' ]


name = ''
for data_set in data_sets:
  name += data_set + '_'
name = name[:-1] 

field = 'P(k)+T0'
ps_data_dir = 'lya_statistics/data/'
output_dir = root_dir + f'fit_results_{field}_{name}/'
create_directory( output_dir )

# load_mcmc_results = False
load_mcmc_results = True


SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Load_Grid_Analysis_Data()
ps_range = SG.Get_Power_Spectrum_Range( kmax=0.01 )
sim_ids = SG.sim_ids

z_min = 2.0
z_max = 5.0 

comparable_data = Get_Comparable_Power_Spectrum_T0( ps_data_dir, z_min, z_max, data_sets, ps_range  )
comparable_grid =  Get_Comparable_Power_Spectrum_T0_from_Grid( comparable_data, SG )



stats_file = output_dir + 'fit_mcmc.pkl'
samples_file = output_dir + 'samples_mcmc.pkl'

fields = [ 'P(k)' ]
z_vals = [ 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.6, 5.0,  ]
data_grid_ps = Get_Data_Grid_Power_spectrum( z_vals, SG )
data_grid_T0 = Get_Data_Grid( ['T0'], SG )


params = SG.parameters

if load_mcmc_results:
  print( f'Loading File: {stats_file}')
  stats = pickle.load( open( stats_file, 'rb' ) )
  param_stats = {}
  for p_id in params.keys():
    p_name = params[p_id]['name']
    p_stats = stats[p_name]
    params[p_id]['mean'] = p_stats['mean']
    params[p_id]['sigma'] = p_stats['standard deviation']
  print( f'Loading File: {samples_file}')
  param_samples = pickle.load( open( samples_file, 'rb' ) )

else:
  nIter = 200000 
  nBurn = nIter / 5
  nThin = 1
  # model, params_mcmc = mcmc_model_3D( comparable_data, comparable_grid, field, 'mean', SG )
  model, params_mcmc = mcmc_model_4D( comparable_data, comparable_grid, field, 'mean', SG )
  MDL = pymc.MCMC( model )  
  MDL.sample( iter=nIter, burn=nBurn, thin=nThin )
  stats = MDL.stats()
  param_stats = {}
  for p_id in params.keys():
    p_name = params[p_id]['name']
    p_stats = stats[p_name]
    params[p_id]['mean'] = p_stats['mean']
    params[p_id]['sigma'] = p_stats['standard deviation']
  Plot_MCMC_Stats( stats, MDL, params_mcmc,  stats_file, output_dir, plot_corner=False )
  param_samples = Write_MCMC_Results( stats, MDL, params_mcmc,  stats_file, samples_file,  output_dir  )




# labels = { 'scale_He':r'$\beta_{\mathrm{He}}$', 'scale_H':r'$\beta_{\mathrm{H}}$', 'deltaZ_He':r'$\Delta z_{\mathrm{He}}$'    }
labels = { 'scale_He':r'$\beta_{\mathrm{He}}$', 'scale_H':r'$\beta_{\mathrm{H}}$', 'deltaZ_He':r'$\Delta z_{\mathrm{He}}$', 'deltaZ_H':r'$\Delta z_{\mathrm{H}}$'    }
Plot_Corner( param_samples, labels, output_dir  )





hpi_sum = 0.95
# Obtain distribution of the power spectrum
# samples_ps = Sample_Power_Spectrum_from_Trace( param_samples, data_grid_ps, SG, hpi_sum=hpi_sum, n_samples=1000 )
# Plot_Power_Spectrum_Sampling( samples_ps, ps_data_dir, output_dir, scales='large', system=system )

# Obtain distribution of the temperature
samples_T0 = Sample_T0_from_Trace( param_samples, data_grid_T0, SG, hpi_sum=hpi_sum, n_samples=1000 )

samples = samples_T0

# Load T0 and tau data
comparable_data = Get_Comparable_Composite_T0_tau()

import pylab
import matplotlib
import matplotlib.font_manager
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'

if system == 'Lux':      prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/brvillas/fonts', "Helvetica.ttf"), size=12)
if system == 'Shamrock': prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=12)

nrows = 1
ncols = 1

font_size = 18
label_size = 16
alpha = 0.5

c_pchw18 = pylab.cm.viridis(.7)
c_hm12 = pylab.cm.cool(.3)

c_boss = pylab.cm.viridis(.3)
c_walther = pylab.cm.viridis(.3)
c_viel = 'C1'
c_boera = pylab.cm.Purples(.7)

text_color  = 'black'
color_line = c_pchw18
color_data = c_boss

fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))

obs_name = 'T0'
z = samples['z']
mean = samples['mean']
high = samples['higher']
low = samples['lower']
ax.plot( z, mean, zorder=1 )
ax.fill_between( z, high, low, alpha=alpha, zorder=1 )  

data_set = comparable_data[obs_name]
data_z = data_set['z']
data_mean = data_set['mean'] 
data_error = data_set['sigma'] 
ax.errorbar( data_z, data_mean, yerr=data_error, fmt='none',  alpha=0.8, ecolor= color_data, zorder=2)
ax.scatter( data_z, data_mean, label='Gaikwad et al.', alpha=0.8, color= color_data, zorder=2) 

ax.tick_params(axis='both', which='major', direction='in', labelsize=label_size )
ax.tick_params(axis='both', which='minor', direction='in' )
ax.set_ylabel( r'$T_0   \,\,\, [\,\mathrm{K}\,]$', fontsize=font_size  )
ax.set_xlabel( r'$z$', fontsize=font_size )
leg = ax.legend(loc=1, frameon=False, fontsize=font_size, prop=prop)
ax.set_xlim( 1.8, 12 )
ax.set_ylim( 3000, 18000)

figure_name = output_dir + f'fig_T0_sampling.png'
fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
print( f'Saved Figure: {figure_name}' )



