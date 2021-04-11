import os, sys
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from parameters_UVB_rates import param_UVB_Rates
from simulation_grid import Simulation_Grid
from simulation_parameters import *
from normalize_power_spectrum  import Normaliz_Flux_Power_Spectrum
from data_optical_depth import Compute_analytical_TauEff_Becker
from flux_power_spectrum import get_skewer_flux_power_spectrum






def ReScale_Optical_Depth_To_F_Mean_Diff( alpha, F_mean, tau_los  ):
  print(alpha)
  tau_los_rescaled = tau_los * alpha
  F_los_rescaled = np.exp( - tau_los_rescaled )
  F_mean_rescaled = F_los_rescaled.mean()
  diff = F_mean_rescaled - F_mean
  return diff


def ReScale_Optical_Depth_To_F_Mean( tau_los, F_mean ):
  from scipy import optimize
  guess = tau_eff / tau_los.mean()
  alpha = optimize.newton(ReScale_Optical_Depth_To_F_Mean_Diff, guess, args=(F_mean, tau_los ) ) 
  tau_los_rescaled = alpha * tau_los
  F_los_rescaled = np.exp( - tau_los_rescaled )
  F_mean_rescaled = F_los_rescaled.mean()
  diff = np.abs( F_mean_rescaled - F_mean ) / F_mean
  if diff > 1e-6: print( 'WARNING: Rescaled F_mean mismatch: {F_mean_rescaled}   {f_mean}')
  return  tau_los_rescaled




sim_id = 0


SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Load_Grid_Analysis_Data( sim_ids = [sim_id], load_fit=False )
ps_out_dir = SG.root_dir + 'flux_power_spectrum_files/'

# sim_name = SG.Grid[sim_id]['key']
# output_dir = ps_out_dir + sim_name + '/'
output_dir = SG.root_dir + 'figures/rescaled_ps/'
create_directory( output_dir )

data_sim  = SG.Grid[sim_id]['analysis']
data_ps  = SG.Grid[sim_id]['analysis']['power_spectrum']
available_indices = data_sim['ps_available_indices'] 
# available_indices = [30]





sim_dir = SG.Get_Simulation_Directory( sim_id )
analysis_dir = sim_dir + 'analysis_files/'


n_file = 55
for n_file in available_indices:


  file_name = analysis_dir + f'{n_file}_analysis.h5'
  print( f'Loading File: {file_name}' )
  file = h5.File( file_name, 'r' )
  current_z = file.attrs['current_z'][0]
  print( f' current_z: {current_z}' )


  lya_data = file['lya_statistics']
  ps_data = lya_data['power_spectrum']


  skewers_key = 'skewers_x'

  skewers_data = lya_data['skewers_x']
  vel_Hubble_x = skewers_data['vel_Hubble'][...]
  flux_HI_x    = skewers_data['los_transmitted_flux_HI'][...]

  skewers_data = lya_data['skewers_y']
  vel_Hubble_y = skewers_data['vel_Hubble'][...]
  flux_HI_y    = skewers_data['los_transmitted_flux_HI'][...]

  skewers_data = lya_data['skewers_z']
  vel_Hubble_z = skewers_data['vel_Hubble'][...]
  flux_HI_z    = skewers_data['los_transmitted_flux_HI'][...]

  n_skewers = lya_data.attrs['n_skewers']
  F_mean_HI = lya_data.attrs['Flux_mean_HI']

  k_vals_0 = ps_data['k_vals'][...]
  ps_mean_0 = ps_data['p(k)'][...]
  indices = ps_mean_0 > 1e-10
  k_vals_sim  = k_vals_0[indices]
  ps_mean_sim = ps_mean_0[indices]


  vel_Hubble = vel_Hubble_x
  F_los_all = np.concatenate([ flux_HI_x, flux_HI_y, flux_HI_z ])


  tau_eff_simulation = -np.log(F_mean_HI[0])
  tau_eff_Becker =  Compute_analytical_TauEff_Becker( current_z )

  print(f' Tau Effective Simulation: {tau_eff_simulation}' )
  print(f' Tau Effective Becker:     {tau_eff_Becker}' )



  d_log_k = 0.1

  tau_eff = tau_eff_simulation
  F_mean = np.exp( -tau_eff )

  # Substitute the F=0 pixels 
  F_min = 1e-200
  F_los_all[ F_los_all < F_min ] = F_min


  alpha_vals = [ 1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7 ]
  ps_mean_rescaled = {}

  for alpha in alpha_vals:
    print( f'ReScaling P(k) by: {alpha}' )
    flux_ps_all_rescaled = []
    # Compute the rescaled Mean Transmitted Flux
    tau_los_all = -np.log( F_los_all )
    tau_los_all_rescaled = alpha * tau_los_all
    F_los_all_rescaled = np.exp( -tau_los_all_rescaled )
    F_mean_rescaled = F_los_all_rescaled.mean()

    # los_id = 0 
    # F_los_rescaled = F_los_all_rescaled[los_id]
    for F_los_rescaled in F_los_all_rescaled:
      delta_F = ( F_los_rescaled - F_mean_rescaled ) / F_mean_rescaled
      k_vals, flux_power_spectrum = get_skewer_flux_power_spectrum( vel_Hubble, delta_F, d_log_k=d_log_k )
      flux_ps_all_rescaled.append( flux_power_spectrum )

    flux_ps_all_rescaled = np.array( flux_ps_all_rescaled )  
    ps_mean_rescaled_alpha = flux_ps_all_rescaled.mean( axis=0 )  
    ps_mean_rescaled[alpha] = ps_mean_rescaled_alpha


  diff_ps = {}
  for alpha in alpha_vals:
    if alpha == 1.0: continue
    print( f'Computing Diff for: {alpha}' )
    ps_resc = ps_mean_rescaled[alpha] 
    ps_0 = ps_mean_rescaled[1.0]
    diff = ( ps_resc - ps_0 ) / ps_0
    diff_ps[alpha] = diff
    
    

  # Plot the fractional differences 
  system = 'Shamrock'

  import matplotlib
  matplotlib.rcParams['mathtext.fontset'] = 'cm'
  matplotlib.rcParams['mathtext.rm'] = 'serif'
  if system == 'Lux':      prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/brvillas/fonts', "Helvetica.ttf"), size=12)
  if system == 'Shamrock': prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=12)



  nrows, ncols = 1, 1
  fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,5*nrows) )

  for alpha in alpha_vals:
    if alpha == 1.0: continue
    diff = diff_ps[alpha]
    label  = r'$\alpha\, = \, {0:.1f} $'.format(alpha)
    ax.plot( k_vals, diff, label=label )


  ax.axhline( y=0, color='r', ls='--' )

  ax.text(0.95, 0.95, r'$z={0:.1f}$'.format(current_z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=13 ) 

  legend_loc = 2
  ax.legend( loc=legend_loc, frameon=False, prop=prop, ncol=3 )
  ax.set_ylabel( r'$ \Delta P\,(k) / P\,(k) $',  fontsize=12 )
  ax.set_xlabel( r'$ k   \,\,\,  [\mathrm{s}\,\mathrm{km}^{-1}] $',  fontsize=12 )
  ax.set_ylim( -0.5, 0.5 )
  ax.set_xscale('log')

  figure_name = output_dir + f'fig_rescaled_power_spectrum_{n_file}.png'
  fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
  print( f'Saved Figure: {figure_name}' )


