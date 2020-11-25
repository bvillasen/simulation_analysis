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
import palettable

output_dir = root_dir + 'fit_results_composite/'
create_directory( output_dir )

load_mcmc_stats = False


SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Load_Grid_Analysis_Data()
sim_ids = SG.sim_ids

comparable_data = Get_Comparable_Composite_T0_tau()
comparable_grid = Get_Comparable_Composite_T0_tau_from_Grid( comparable_data, SG )


field = 'T0+tau'
stats_file = output_dir + 'fit_mcmc.pkl'

fields = [ 'T0', 'tau' ]
data_grid = Get_Data_Grid( fields, SG )
  
params = SG.parameters

if load_mcmc_stats:
  print( f'Loading File: {stats_file}')
  stats = pickle.load( open( stats_file, 'rb' ) )
  param_stats = {}
  for p_id in params.keys():
    p_name = params[p_id]['name']
    p_stats = stats[p_name]
    params[p_id]['mean'] = p_stats['mean']
    params[p_id]['sigma'] = p_stats['standard deviation']

else:
  nIter = 200000
  nBurn = nIter / 5
  nThin = 1
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


  # Plot_MCMC_Stats( stats, MDL, params_mcmc,  stats_file, output_dir )

  samples = {} 
  for p_id in params_mcmc.keys():
    param = params_mcmc[p_id]
    samples[p_id] = {}
    samples[p_id]['name'] = param['name']
    samples[p_id]['trace'] = param['sampler'].trace() 


labels = { 'scale_He':r'$\beta_{\mathrm{He}}$', 'scale_H':r'$\beta_{\mathrm{H}}$', 'deltaZ_He':r'$\Delta z_{\mathrm{He}}$', 'deltaZ_H':r'$\Delta z_{\mathrm{H}}$'    }
    
param_ids = samples.keys()

n_param = len( param_ids )
color = 'C0'
data_color = 'C9'
font_size = 15
label_size = 18
alpha = 0.6
fig_size = 5
space = 0.05
tick_label_size = 12
tick_length = 7
tick_width = 2
border_width = 2.0
n_bins_1D = 20
hist_1D_line_width = 2
hist_1D_line_color = 'C0'

n_bins_2D = 30
hist_2D_colormap = palettable.cmocean.sequential.Ice_20_r.mpl_colormap


fig, ax_l = plt.subplots(nrows=n_param, ncols=n_param, figsize=(fig_size*n_param,fig_size*n_param))
fig.subplots_adjust( wspace=space, hspace=space )

for j in range( n_param ):
  for i in range( n_param ):
    
    ax = ax_l[j][i]
    plot_y_lables, plot_x_lables = False, False
    plot_y_ticks  = True
    if i == 0: plot_y_lables = True
    if j == n_param-1: plot_x_lables = True 
    if i == j: 
      plot_y_lables = False
      plot_y_ticks = False
    
    if i > j:
      ax.axis("off")
      continue
    
    if plot_x_lables: ax.tick_params(axis='x', which='major', direction='in', labelsize=tick_label_size, length=tick_length, width=tick_width )
    else:             ax.tick_params(axis='x', which='major', direction='in', labelsize=0, length=tick_length, width=tick_width )
    if plot_y_lables: ax.tick_params(axis='y', which='major', direction='in', labelsize=tick_label_size, length=tick_length, width=tick_width )
    else:             ax.tick_params(axis='y', which='major', direction='in', labelsize=0, length=tick_length, width=tick_width )
    if not plot_y_ticks: ax.tick_params(axis='y', which='major', length=0 )
    
    
    
    if plot_y_lables:
      if j == 0: y_label = ''
      else: y_label = labels[samples[j]['name']]  
      ax.set_ylabel( y_label, fontsize=label_size )
    
    if plot_x_lables:
      x_label = labels[samples[i]['name']]  
      ax.set_xlabel( x_label, fontsize=label_size )
      
    if i == j: plot_type = '1D'
    if i < j:  plot_type = '2D'
    
    
    if plot_type == '1D':
      name  = samples[j]['name']
      trace = samples[j]['trace']
      hist, bin_edges = np.histogram( trace, bins=n_bins_1D ) 
      bin_centers = ( bin_edges[:-1] + bin_edges[1:] ) / 2.
      bin_width = bin_centers[0] - bin_centers[1]  
      ax.step( bin_centers, hist, where='mid', color=hist_1D_line_color, linewidth=hist_1D_line_width  )
      
    if plot_type == '2D':
      trace_y = samples[j]['trace']
      trace_x = samples[i]['trace']
      hist, x_edges, y_edges = np.histogram2d( trace_x, trace_y, bins=[n_bins_2D, n_bins_2D] )
      hist = hist.astype( np.float ) 
      hist = hist.T / hist.sum()
      hist_mean = hist.mean()
      hist_sigma = np.sqrt( ((hist-hist_mean)**2).mean() )
      extent = [ trace_x.min(), trace_x.max(), trace_y.min(), trace_y.max() ]
      ax.imshow( hist[::-1], cmap=hist_2D_colormap, extent=extent, aspect='auto' )
      ax.contour( hist, [hist_sigma, 2* hist_sigma], extent=extent )
          
    [sp.set_linewidth(border_width) for sp in ax.spines.values()]


figure_name = output_dir + 'corner.png'
fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
print( f'Saved Figure: {figure_name}' )





# # Obtain distribution of observables
# n_samples = 10000
# observables = [ 'T0', 'tau' ]
# observables_samples = Sample_Observables( n_samples, observables, params, data_grid, SG  )
# Plot_Observables( observables_samples, comparable_data, params, SG, 'sampling', output_dir)











