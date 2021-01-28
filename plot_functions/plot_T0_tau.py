import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import palettable
import pylab
import pickle
from matplotlib.legend_handler import HandlerTuple
import os, sys
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from tools import *
from mcmc_data_functions import Get_Comparable_Composite_T0_tau


def plot_T0_and_tau( output_dir, sim_data_sets=None, system=None ):

  # Load T0 and tau data
  comparable_data = Get_Comparable_Composite_T0_tau()

  import matplotlib
  import matplotlib.font_manager
  matplotlib.rcParams['mathtext.fontset'] = 'cm'
  matplotlib.rcParams['mathtext.rm'] = 'serif'

  if system == 'Lux':      prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/brvillas/fonts', "Helvetica.ttf"), size=12)
  if system == 'Shamrock': prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=12)

  nrows = 1
  ncols = 2

  font_size = 18
  label_size = 16
  alpha = 0.8

  c_pchw18 = pylab.cm.viridis(.7)
  c_hm12 = pylab.cm.cool(.3)

  c_boss = pylab.cm.viridis(.3)
  c_walther = pylab.cm.viridis(.3)
  c_viel = 'C1'
  c_boera = pylab.cm.Purples(.7)

  text_color  = 'black'
  color_line = c_pchw18
  color_data = c_boss

  fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))

  ax = ax_l[0]
  obs_name = 'T0'
  
  for sim_data in sim_data_sets:
    z = sim_data['z']
    vals = sim_data[obs_name]
    ax.plot( z, vals , label=sim_data['plot_label'], zorder=1 )
    # ax.plot( z, vals , c=color_line, label=sim_data['plot_label'], zorder=1 )


  data_set = comparable_data[obs_name]
  data_z = data_set['z']
  data_mean = data_set['mean'] 
  data_error = data_set['sigma'] 
  ax.errorbar( data_z, data_mean, yerr=data_error, fmt='none',  alpha=0.8, ecolor= color_data, zorder=2)
  ax.scatter( data_z, data_mean, label='Data for MCMC Fit', alpha=0.8, color= color_data, zorder=2) 

  ax.tick_params(axis='both', which='major', direction='in', labelsize=label_size )
  ax.tick_params(axis='both', which='minor', direction='in' )
  ax.set_ylabel( r'$T_0   \,\,\, [\,\mathrm{K}\,]$', fontsize=font_size  )
  ax.set_xlabel( r'$z$', fontsize=font_size )
  leg = ax.legend(loc=1, frameon=False, fontsize=font_size, prop=prop)
  ax.set_xlim( 2, 12 )
  ax.set_ylim( 3000, 18000)


  ax = ax_l[1]
  obs_name = 'tau'


  for sim_data in sim_data_sets:
    z = sim_data['z']
    vals = sim_data[obs_name]
    ax.plot( z, vals,  label=sim_data['plot_label'], zorder=1 )
    # ax.plot( z, vals, c=color_line, label=sim_data['plot_label'], zorder=1 )

  data_set = comparable_data[obs_name]
  data_z = data_set['z']
  data_mean = data_set['mean'] 
  data_error = data_set['sigma'] 
  ax.errorbar( data_z, data_mean, yerr=data_error, fmt='none',  alpha=0.8, ecolor= color_data, zorder=2)
  ax.scatter( data_z, data_mean, label='Data for MCMC Fit', alpha=0.8, color= color_data, zorder=2) 

  ax.tick_params(axis='both', which='major', direction='in', labelsize=label_size )
  ax.tick_params(axis='both', which='minor', direction='in' )
  ax.set_ylabel( r'$\tau_{eff}$', fontsize=font_size  )
  ax.set_xlabel( r'$z$', fontsize=font_size )
  leg = ax.legend(loc=2, frameon=False, fontsize=font_size, prop=prop)
  ax.set_xlim( 2, 6 )
  ax.set_ylim( 0.1, 8)
  ax.set_yscale('log')

  figure_name = output_dir + f'fig_T0_tau.png'
  fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
  print( f'Saved Figure: {figure_name}' )




def plot_tau_HeII( output_dir, sim_data_sets=None, system=None ):

  # Load T0 and tau data
  # comparable_data = Get_Comparable_Composite_T0_tau()

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
  alpha = 0.8

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

  
  obs_name = 'tau_HeII'
  for sim_data in sim_data_sets:
    z = sim_data['z']
    vals = sim_data[obs_name]
    ax.plot( z, vals , label=sim_data['plot_label'], zorder=1 )
    # ax.plot( z, vals , c=color_line, label=sim_data['plot_label'], zorder=1 )

  # 
  # data_set = comparable_data[obs_name]
  # data_z = data_set['z']
  # data_mean = data_set['mean'] 
  # data_error = data_set['sigma'] 
  # ax.errorbar( data_z, data_mean, yerr=data_error, fmt='none',  alpha=0.8, ecolor= color_data, zorder=2)
  # ax.scatter( data_z, data_mean, label='Data for MCMC Fit', alpha=0.8, color= color_data, zorder=2) 

  ax.tick_params(axis='both', which='major', direction='in', labelsize=label_size )
  ax.tick_params(axis='both', which='minor', direction='in' )
  ax.set_ylabel( r'$\tau_{eff} \,\, \mathrm{HeII}$', fontsize=font_size  )
  ax.set_xlabel( r'$z$', fontsize=font_size )
  leg = ax.legend(loc=1, frameon=False, fontsize=font_size, prop=prop)
  ax.set_xlim( 2, 6 )
  ax.set_ylim( 1, 10 )



  figure_name = output_dir + f'fig_tau_HeII.png'
  fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
  print( f'Saved Figure: {figure_name}' )



