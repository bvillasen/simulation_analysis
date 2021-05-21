import os, sys
from pathlib import Path
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
cosmo_dir = os.path.dirname(os.getcwd()) + '/'
sys.path.append(cosmo_dir + 'tools')
from tools import *
#Append analysis directories to path
extend_path()
from load_data import load_analysis_data
from phase_diagram_functions import fit_thermal_parameters_mcmc, get_density_temperature_values_to_fit



# data_dir = '/home/bruno/Desktop/data/'
# data_dir = '/home/brunoq/Desktop/ssd_0/data/'
data_dir = '/raid/bruno/data/'

sim_name = '1024_50Mpc'
sim_name = '2048_100Mpc'
input_dir = data_dir + f'cosmo_sims/rescaled_P19/{sim_name}/analysis_files/'
output_dir = data_dir + f'cosmo_sims/rescaled_P19/figures/phase_diagram_{sim_name}/fit/'
create_directory( output_dir )

n_file = 55
for n_file in range(3,56):
  data = load_analysis_data( n_file, input_dir )
  z = data['cosmology']['current_z'] 
  data_pd = data['phase_diagram']
  fit_mcmc = Load_Pickle_Directory( input_dir + f'fit_mcmc/fit_{n_file}.pkl' )

  values_to_fit = get_density_temperature_values_to_fit( data['phase_diagram'], delta_min=-1.5, delta_max=3, n_samples_line=50, fraction_enclosed=0.67 )




  dens = values_to_fit['density']
  temp = values_to_fit['temperature']
  temp_sigma_p = values_to_fit['temperature_sigma_p']
  temp_sigma_m = values_to_fit['temperature_sigma_l']
  temp_error = np.array([ temp_sigma_m, temp_sigma_p ])


  fit_all = {}
  fit_all[0] = { 'xmin':-1, 'xmax':0, 'order':1 }
  fit_all[1] = { 'xmin':0, 'xmax':1,  'order':1 }
  fit_all[2] = { 'xmin':-1, 'xmax':1,  'order':2 }


  for index in fit_all.keys():
    fit = fit_all[index]
    xmin, xmax = fit['xmin'], fit['xmax']
    order = fit['order']
    indices_0 = ( dens >= xmin ) * ( dens <= xmax )
    fit_res =  np.polyfit( dens[indices_0], temp[indices_0], order   )[::-1]
    x_vals = np.linspace( xmin, xmax, 100 )
    y_vals = np.zeros_like( x_vals )
    # fit_res[0]*x_vals + fit_res[1]
    fit['fit_res'] = fit_res
    for i,coeff in enumerate( fit_res ):
      y_vals +=  coeff * x_vals**i
    fit['xvals'], fit['yvals'] = x_vals, y_vals
    if order == 1: label = r'$T_0 \,= \, {0:.2f} \times 10^4 \,\,\, \gamma \, = \,{1:.2f} $'.format( 10**fit_res[0]/ 10**4 , fit_res[1]+1) 
    if order == 2:label = r'Parabolic: $T_0 \,= \, {0:.2f} \times 10^4 \,\,\, C_1\,=\, {1:.2f} \,\,\, C_2\,=\, {2:.2f}  $'.format( 10**fit_res[0]/ 10**4, fit_res[1], fit_res[2] ) 
    fit['label'] = label
    


  import matplotlib
  matplotlib.rcParams['mathtext.fontset'] = 'cm'
  matplotlib.rcParams['mathtext.rm'] = 'serif'


  system = 'Shamrock'
  if system == 'Lux':      prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/brvillas/fonts', "Helvetica.ttf"), size=16)
  if system == 'Shamrock': prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=16)


  fig_width = 8
  fig_dpi = 300
  label_size = 18
  figure_text_size = 18
  legend_font_size = 16
  tick_label_size_major = 15
  tick_label_size_minor = 13
  tick_size_major = 5
  tick_size_minor = 3
  tick_width_major = 1.5
  tick_width_minor = 1
  border_width = 1

  text_color = 'black'

  color_lines = [ 'C0', 'C2',  'C4' ]

  color_line = 'C3'

  font_size = 18
  label_size = 16
  alpha = 0.5

  nrows, ncols = 2, 2
  fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))
  plt.subplots_adjust( hspace = 0.1, wspace=0.2)


  for i in range( ncols ):
    for j in range( nrows):
    
      ax = ax_l[j][i]

      if i == 0 and j == 0:
        
        
        T0    = fit_mcmc['T0']['mean']
        gamma = fit_mcmc['gamma']['mean']
        mcmc_x = np.linspace(-1, 1, 10)
        mcmc_y = gamma * mcmc_x  + T0 
        label = r'$T_0 \,= \, {0:.2f} \times 10^4 \,\,\, \gamma \, = \,{1:.2f} $'.format( 10**T0/ 10**4 , gamma+1) 
        ax.plot( mcmc_x, mcmc_y, c=color_line, lw=2.5, label=label, zorder=1  )
        for index in [0,1]:
          fit = fit_all[index]
          ax.plot( fit['xvals'], fit['yvals'], '--', c=color_lines[index], lw=2.5, label=fit['label'], zorder=2 )
              
        ax.errorbar( dens, temp, yerr=temp_error, fmt='o',  ms=4, c='C1', zorder=1)
        # ax.scatter( dens, temp,  s=2, c='C1', zorder=3 )
        
      if i == 1 and j == 0:
        T0    = fit_mcmc['T0']['mean']
        gamma = fit_mcmc['gamma']['mean']
        mcmc_points = gamma * dens + T0 
        exp_mcmc = 10**mcmc_points
        exp_data = 10**temp 
        y =  ( exp_data - exp_mcmc ) / exp_mcmc
        exp_temp_m =  10**(temp - temp_error[0])
        exp_temp_p =  10**(temp + temp_error[1])
        exp_error = [np.abs( exp_temp_m - exp_mcmc ) / exp_mcmc,  np.abs( exp_temp_p - exp_mcmc ) / exp_mcmc ]
        label = r'$T_0 \,= \, {0:.2f} \times 10^4 \,\,\, \gamma \, = \,{1:.2f} $'.format( 10**T0/ 10**4 , gamma+1)
        ax.plot( [-1, 1], [0, 0], lw=2.5, label=label, color='C3',   )
        ax.errorbar( dens, y , yerr=exp_error, fmt='o',  ms=4, c='C1', zorder=1 )
        
        for index in [0,1]:
          fit = fit_all[index]
          fit_res = fit['fit_res']
          xmin, xmax = fit['xmin'], fit['xmax']
          x = np.linspace( xmin, xmax, 100 )
          y_mcmc = gamma * x + T0
          y_fit = fit_res[0] + fit_res[1]*x
          temp_mcmc = 10**y_mcmc
          temp_fit  = 10**y_fit 
          diff = ( temp_fit - temp_mcmc ) / temp_mcmc
          ax.plot( x, diff, '--', c=color_lines[index], lw=2.5, label=fit['label'], zorder=2 )
      
          x = np.linspace( -1, 1, 100 )
          y_mcmc = gamma * x + T0
          y_fit = fit_res[0] + fit_res[1]*x
          temp_mcmc = 10**y_mcmc
          temp_fit  = 10**y_fit 
          diff = ( temp_fit - temp_mcmc ) / temp_mcmc
          ax.plot( x, diff, '--', c=color_lines[index], lw=2.5, zorder=2, alpha=0.3 )
          
      if i == 0 and j == 1:
        
        
        T0    = fit_mcmc['T0']['mean']
        gamma = fit_mcmc['gamma']['mean']
        mcmc_x = np.linspace(-1, 1, 10)
        mcmc_y = gamma * mcmc_x  + T0 
        label = r'$T_0 \,= \, {0:.2f} \times 10^4 \,\,\, \gamma \, = \,{1:.2f} $'.format( 10**T0/ 10**4 , gamma+1) 
        ax.plot( mcmc_x, mcmc_y, c=color_line, lw=2.5, label=label, zorder=1  )
        for index in [2]:
          fit = fit_all[index]
          ax.plot( fit['xvals'], fit['yvals'], '--', c=color_lines[index], lw=2.5, label=fit['label'], zorder=2 )
              
        ax.errorbar( dens, temp, yerr=temp_error, fmt='o',  ms=4, c='C1', zorder=1)
        # ax.scatter( dens, temp,  s=2, c='C1', zorder=3 )
      
      if i == 1 and j == 1:
        T0    = fit_mcmc['T0']['mean']
        gamma = fit_mcmc['gamma']['mean']
        mcmc_points = gamma * dens + T0 
        exp_mcmc = 10**mcmc_points
        exp_data = 10**temp 
        y =  ( exp_data - exp_mcmc ) / exp_mcmc
        exp_temp_m =  10**(temp - temp_error[0])
        exp_temp_p =  10**(temp + temp_error[1])
        exp_error = [np.abs( exp_temp_m - exp_mcmc ) / exp_mcmc,  np.abs( exp_temp_p - exp_mcmc ) / exp_mcmc ]
        label = r'$T_0 \,= \, {0:.2f} \times 10^4 \,\,\, \gamma \, = \,{1:.2f} $'.format( 10**T0/ 10**4 , gamma+1)
        ax.plot( [-1, 1], [0, 0], lw=2.5, label=label, color='C3',   )
        ax.errorbar( dens, y , yerr=exp_error, fmt='o',  ms=4, c='C1', zorder=1 )
        
        for index in [2]:
          fit = fit_all[index]
          fit_res = fit['fit_res']
          xmin, xmax = fit['xmin'], fit['xmax']
          x = np.linspace( xmin, xmax, 100 )
          y_mcmc = gamma * x + T0
          y_fit = fit_res[0] + fit_res[1]*x + fit_res[2]*x**2 
          temp_mcmc = 10**y_mcmc
          temp_fit  = 10**y_fit 
          diff = ( temp_fit - temp_mcmc ) / temp_mcmc
          ax.plot( x, diff, '--', c=color_lines[index], lw=2.5, label=fit['label'], zorder=2 )
      
      

        
        
      if i == 0: ax.set_ylabel(r'$\log_{10} \, T \,\,[\,\mathrm{K}\,]$', fontsize=label_size , color=text_color)
      if i == 1: ax.set_ylabel(r'$\Delta T / T $', fontsize=label_size , color=text_color)
      ax.set_xlabel(r'$\log_{10} \, \Delta$ ', fontsize=label_size , color=text_color )

      text  = r'$z = {0:.1f}$'.format( z ) 
      ax.text(0.03, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size, color=text_color )

    


      ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in')
      ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')
      [sp.set_linewidth(border_width) for sp in ax.spines.values()]
      x_min, x_max = -1.6, 2.2
      if i == 0: y_min, y_max =  3.4, 5.0
      if i == 1: y_min, y_max =  -1, 1
      ax.set_xlim( x_min, x_max )
      ax.set_ylim( y_min, y_max )

      ax.legend( loc=4, frameon=False, prop=prop)

  figure_name = output_dir + f'pd_fit_{n_file}.png'
  fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
  print( f'Saved Figure: {figure_name}' )

