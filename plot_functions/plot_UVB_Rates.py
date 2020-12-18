import os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

def Plot_Grid_UVB_Rates( SG, output_dir, ids_to_plot=None, plot_label=None ):


  nrows = 2
  ncols = 3
  fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))


  font_size = 15


  sim_ids = SG.Grid.keys()
  if ids_to_plot == None: ids_to_plot = sim_ids
  root_keys = [ 'Chemistry', 'Photoheating' ]
  
  second_keys = {'Chemistry':['k24', 'k26', 'k25'], 'Photoheating':['piHI', 'piHeI', 'piHeII', ] }

  for j, root_key in enumerate(root_keys):

    for sim_id in ids_to_plot:
      sim_rates = SG.Grid[sim_id]['UVB_rates']
      rates = sim_rates[root_key]
      
      for i,key in enumerate(second_keys[root_key]):
        ax = ax_l[j][i]
        z = sim_rates['z'][...]
        rates_data = sim_rates[root_key][key][...]
        label = ''
        if plot_label != None:
          val = SG.Grid[sim_id]['parameters'][plot_label]
          label = plot_label + ' = {0:.2f}'.format(val)
          
        ax.plot( z, rates_data, label = label )
        
        if plot_label: ax.legend( frameon = False, fontsize=font_size )
        
    for i,key in enumerate(rates.keys()):
      ax = ax_l[j][i]
      ax.set_ylabel( key, fontsize=font_size  )
      ax.set_xlabel( r'$z$', fontsize=font_size )  
      ax.set_yscale('log')

  figure_name = output_dir + 'grid_UVB_rates.png'
  fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
  print( f'Saved Figure: {figure_name}' )



