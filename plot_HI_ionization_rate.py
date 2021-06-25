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
from mcmc_plotting_functions import *
from mcmc_sampling_functions import *
from data_photoionization_HI import data_photoionization_HI_gaikwad_2017, data_photoionization_HI_becker_bolton_2013, data_photoionization_HI_dalosio_2018, data_photoionization_HI_gallego_2021, data_photoionization_HI_calverley_2011, data_photoionization_HI_wyithe_2011
from colors import *

data_name = 'fit_results_P(k)+tau_HeII_Boss'



ps_data_dir = 'lya_statistics/data/'
mcmc_dir = root_dir + 'fit_mcmc/'
input_dir = mcmc_dir + f'{data_name}/observable_samples/' 
output_dir = mcmc_dir + f'{data_name}/observable_figures/'
create_directory( output_dir )
rescale_walter_file = ps_data_dir + 'rescale_walther_to_boss.pkl' 


stats_file = input_dir + 'fit_mcmc.pkl'
samples_file = input_dir + 'samples_mcmc.pkl'

params = Load_Pickle_Directory( input_dir + 'params.pkl' )

print( f'Loading File: {stats_file}')
stats = pickle.load( open( stats_file, 'rb' ) )
for p_id in params.keys():
  p_name = params[p_id]['name']
  p_stats = stats[p_name]
  params[p_id]['mean'] = p_stats['mean']
  params[p_id]['sigma'] = p_stats['standard deviation']
print( f'Loading File: {samples_file}')
param_samples = pickle.load( open( samples_file, 'rb' ) )


# Get the Highest_Likelihood parameter values 
params_HL = Get_Highest_Likelihood_Params( param_samples, n_bins=100 )

# # Obtain distribution of the power spectrum
# file_name = input_dir + 'samples_power_spectrum.pkl'
# samples_ps = Load_Pickle_Directory( file_name )
# 
# # Obtain distribution of the other fields
# file_name = input_dir + 'samples_fields.pkl' 
# field_list = ['T0', 'tau', 'tau_HeII']
# samples_fields = Load_Pickle_Directory( file_name )

# Obtain distribution of the UVBRates
file_name = input_dir + 'samples_uvb_rates_new.pkl' 
samples_uvb_rates = Load_Pickle_Directory( file_name )

z = samples_uvb_rates['photoheating_HI']['z']
heat_HI   = samples_uvb_rates['photoheating_HI']
heat_HeI  = samples_uvb_rates['photoheating_HeI']
heat_HeII = samples_uvb_rates['photoheating_HeII']
ion_HI    = samples_uvb_rates['photoionization_HI']
ion_HeI   = samples_uvb_rates['photoionization_HeI']
ion_HeII  = samples_uvb_rates['photoionization_HeII']

heat_HI_m = heat_HI['Highest_Likelihood']
heat_HI_l = heat_HI['lower']
heat_HI_h = heat_HI['higher']

heat_HeI_m = heat_HeI['Highest_Likelihood']
heat_HeI_l = heat_HeI['lower']
heat_HeI_h = heat_HeI['higher']

heat_HeII_m = heat_HeII['Highest_Likelihood']
heat_HeII_l = heat_HeII['lower']
heat_HeII_h = heat_HeII['higher']

ion_HI_m = ion_HI['Highest_Likelihood']
ion_HI_l = ion_HI['lower']
ion_HI_h = ion_HI['higher']

ion_HeI_m = ion_HeI['Highest_Likelihood']
ion_HeI_l = ion_HeI['lower']
ion_HeI_h = ion_HeI['higher']

ion_HeII_m = ion_HeII['Highest_Likelihood']
ion_HeII_l = ion_HeII['lower']
ion_HeII_h = ion_HeII['higher']


data_out_m = [ z, heat_HI_m, heat_HeI_m, heat_HeII_m, ion_HI_m, ion_HeI_m, ion_HeII_m  ]
data_out_h = [ z, heat_HI_h, heat_HeI_h, heat_HeII_h, ion_HI_h, ion_HeI_h, ion_HeII_h  ]
data_out_l = [ z, heat_HI_l, heat_HeI_l, heat_HeII_l, ion_HI_l, ion_HeI_l, ion_HeII_l  ]


file_name = input_dir + 'uvb_rates.txt'
np.savetxt( file_name, data_out_m )

file_name = input_dir + 'uvb_rates_percentile2.5.txt'
np.savetxt( file_name, data_out_l )

file_name = input_dir + 'uvb_rates_percentile97.5.txt'
np.savetxt( file_name, data_out_h )




system = 'Shamrock'

field_list = ['photoheating_HI', 'photoheating_HeI', 'photoheating_HeII', 'photoionization_HI', 'photoionization_HeI', 'photoionization_HeII' ]

titles = [ 'Photoheating HI', 'Photoheating HeI', 'Photoheating HeII', 'Photoionization HI', 'Photoionization HeI', 'Photoionization HeII' ]

y_labels = [ r'HI Photoheating Rate [$\mathrm{eV\, s^{-1}}$]',  r'HeI Photoheating Rate [$\mathrm{eV\, s^{-1}}$]', r'HeII Photoheating Rate [$\mathrm{eV\, s^{-1}}$] ',
             r'HI Photoionization Rate [$\mathrm{s^{-1}}$]',  r'HeI Photoionization Rate [$\mathrm{s^{-1}}$]', r'HeII Photoionization Rate [$\mathrm{s^{-1}}$] ']

from scipy import interpolate as interp 
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'
if system == 'Lux':      prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/brvillas/fonts', "Helvetica.ttf"), size=14)
if system == 'Shamrock': prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=14)


font_size = 18

tick_size_major, tick_size_minor = 6, 4
tick_label_size_major, tick_label_size_minor = 14, 12
tick_width_major, tick_width_minor = 1.5, 1




data_sets = [ data_photoionization_HI_becker_bolton_2013, data_photoionization_HI_dalosio_2018, data_photoionization_HI_calverley_2011, data_photoionization_HI_wyithe_2011]
colors_data = [ 'C1', 'C2', 'C5', 'C4' ]


black_background = True
if black_background:
  text_color = 'white'
  green = greens[4]
  yellow = yellows[2]
  purple = purples[3] 
  blue = yellow_green_blues[-2]
  red = reds[5]
  colors_data = [ 'C1', green, red, yellow,  purple ]
  
  colormap =  palettable.colorbrewer.sequential.Blues_9_r
  colors = colormap.mpl_colors
  color_line = colors[4]
  color_bar  = colors[4]

nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))


rates_data = samples_uvb_rates['photoionization_HI']
rates_z  = rates_data['z']
rates_HL = rates_data['Highest_Likelihood'] / 1e-12
high     = rates_data['higher'] / 1e-12 
lower    = rates_data['lower'] / 1e-12

ax.plot( rates_z, rates_HL,  zorder=1, color=color_line, label='This Work (Modified P19)'  )
ax.fill_between( rates_z, high, lower, alpha=0.6, zorder=1, color=color_bar  )

for i,data_set in enumerate(data_sets):
  data_z = data_set['z']
  data_mean = data_set['mean']
  data_high = data_set['high']
  data_low  = data_set['low']
  data_name = data_set['name']
  yerr = [  data_mean-data_low, data_high-data_mean, ]
  color_data = colors_data[i]
  ax.errorbar( data_z, data_mean, yerr=yerr, fmt='o',  c= color_data, zorder=2, label=data_name)


color_gallego = colors_data[-1]
data_set = data_photoionization_HI_gallego_2021
data_z = data_set['z'][1]
data_mean = data_set['mean'][1]
data_high = data_set['high'][1]
data_low  = data_set['low'][1]
data_name = data_set['name']
yerr = [  [data_mean-data_low], [data_high-data_mean] ]
color_data = colors_data[i]
ax.errorbar( [data_z], [data_mean], yerr=yerr, fmt='o',  c= color_gallego, zorder=2, label=data_name)
indices = [0, 2]
data_z = data_set['z'][indices]
data_mean = data_set['mean'][indices]
data_high  = data_set['mean'][indices]
data_low  = data_set['low'][indices] * 1.2
yerr = [  data_mean-data_low, data_high-data_mean ]
ax.errorbar( data_z, data_mean, yerr=yerr, fmt='o', uplims=True, c=color_gallego, zorder=2, )


color_gaikwad = 'C4'
data_set = data_photoionization_HI_gaikwad_2017
data_z = data_set['z']
data_mean = data_set['mean']
data_high = data_set['high']
data_low  = data_set['low']
data_name = data_set['name']
yerr = [  data_mean-data_low, data_high-data_mean ]
color_data = colors_data[i]
ax.errorbar( data_z, data_mean, yerr=yerr, fmt='o',  c= color_gaikwad, zorder=2, label=data_name)


# ax.text(0.8, 0.95, text, horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=font_size )
legend_loc = 0
leg = ax.legend(  loc=legend_loc, frameon=False, prop=prop    )
for text in leg.get_texts():
  plt.setp(text, color = text_color)

ax.set_yscale('log')
ax.set_xlabel( r'$z$', fontsize=22, color=text_color)
ax.set_ylabel( r'$\Gamma_{\mathrm{HI}} \,\,\,\,\, [\,\mathrm{10^{-12} \,s^{-1} } ]$', fontsize=22 , color=text_color)

ax.tick_params(axis='both', which='major', direction='in', color=text_color, labelcolor=text_color, labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major  )
ax.tick_params(axis='both', which='minor', direction='in', color=text_color, labelcolor=text_color, labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor  )


ax.set_xlim(0.05, 6.5)
ax.set_ylim(0.02, 2.0)

if black_background: 
  fig.patch.set_facecolor('black') 
  ax.set_facecolor('k')
  [ spine.set_edgecolor(text_color) for spine in list(ax.spines.values()) ]


figure_name = output_dir + 'fig_phothoionization_HI_black_new.png'
fig.savefig( figure_name, bbox_inches='tight', dpi=300, facecolor=fig.get_facecolor() )
print( f'Saved Figure: {figure_name}' )


# 
# params_HL = params_HL.flatten()
# beta_He, beta_H, deltaZ_He, deltaZ_H = params_HL
# label =  r'$\beta_{\mathrm{He}}:$' + f'{beta_He:.2f}' + '\n' 
# label += r'$\beta_{\mathrm{H}}:$' + f' {beta_H:.2f}' + '\n' 
# label += r'$\Delta z_{\mathrm{He}}:$' + f'{deltaZ_He:.2f}' + '\n' 
# label += r'$\Delta z_{\mathrm{H}}:$' + f'{deltaZ_H:.2f}' 
# 
# 
# 
# 
# Plot_Power_Spectrum_Sampling( samples_ps, ps_data_dir, output_dir, scales='large',  system=system, label=label )
# Plot_Power_Spectrum_Sampling( samples_ps, ps_data_dir, output_dir, scales='middle', system=system, label=label, rescaled_walther=True, rescale_walter_file=rescale_walter_file )
# Plot_Power_Spectrum_Sampling( samples_ps, ps_data_dir, output_dir, scales='all',    system=system, label=label, rescaled_walther=True, rescale_walter_file=rescale_walter_file )
# 
# 
# 
# Plot_T0_Sampling( samples_fields['T0'], output_dir, system=system, label=label, plot_splines=True )
# 
# Plot_tau_HeII_Sampling( samples_fields, output_dir, system=system, label=label )
