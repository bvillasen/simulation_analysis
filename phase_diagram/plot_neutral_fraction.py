import os, sys, time
from pathlib import Path
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pylab
import palettable
cwd = os.getcwd()
cosmo_dir = cwd[: cwd.find('simulation_analysis')] + 'simulation_analysis/'
tools_dir = cosmo_dir + 'tools'
sys.path.append( tools_dir )
from tools import *
#Append analysis directories to path
extend_path()

data_dir = '/raid/bruno/data/'
input_dir_0 = data_dir + 'cosmo_sims/2048_hydro_50Mpc/neutral_fraction_pchw18/'
input_dir_1 = data_dir + 'cosmo_sims/rescaled_P19/2048_50Mpc/neutral_fraction_P19m/'
output_dir  = data_dir + 'cosmo_sims/rescaled_P19/2048_50Mpc/figures/'
create_directory( output_dir )


dirs = { 'P19': input_dir_0, 'P19m':input_dir_1 }


data_all = {}
for data_name in dirs:
  input_dir = dirs[data_name]
  z_vals, frac_HI, frac_HeI, frac_HeII = [], [], [], []
  for n_file in range( 169 ):
    data = Load_Pickle_Directory( input_dir + f'fraction_{n_file}.pkl')
    z_vals.append( data['current_z'] ) 
    frac_HI.append( data['HI']['dens_fraction'] )
    frac_HeI.append( data['HeI']['dens_fraction'] )
    frac_HeII.append( data['HeII']['dens_fraction'] )
  z_vals = np.array( z_vals )
  frac_HI = np.array( frac_HI )
  frac_HeI = np.array( frac_HeI )
  frac_HeII = np.array( frac_HeII )
  data_all[data_name] = { 'z':z_vals, 'HI':frac_HI, 'HeI':frac_HeI, 'HeII':frac_HeII }
  
  
  
  
black_background = True
if black_background:
  text_color = 'white'

chem_type = 'HeII'
# chem_type = 'HI'

tick_size_major, tick_size_minor = 6, 4
tick_label_size_major, tick_label_size_minor = 14, 12
tick_width_major, tick_width_minor = 1.5, 1

blues = palettable.colorbrewer.sequential.Blues_9_r
blue = blues.mpl_colors[4]
oranges = palettable.colorbrewer.sequential.YlOrBr_9.mpl_colors
orange = oranges[3]

green = pylab.cm.viridis(.7)

colors = {'P19':orange, 'P19m':blue }
labels = { 'P19': 'Puchwein et al. 2019', 'P19m':'This Work (Modified P19)'}

font_size = 16
label_size = 16
legend_size = 12
prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=legend_size )
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'

line_width = 2

ncols, nrows = 1, 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12*ncols,3.*nrows))


for data_name in [ 'P19', 'P19m' ]: 
  z = data_all[data_name]['z']
  if chem_type == 'HI': ionized_frac = 1 - data_all[data_name]['HI'] 
  if chem_type == 'HeII': ionized_frac = 1 - data_all[data_name]['HeI'] - data_all[data_name]['HeII']   

  ionized_frac[ ionized_frac < 0 ] = 0
  ionized_frac[ ionized_frac > 1 ] = 1
  z = z[::-1]
  ionized_frac = ionized_frac[::-1]
  z_interp = np.linspace( z[0], z[-1], 1000000)
  ionized_frac_interp = np.interp( z_interp, z, ionized_frac)
  ionized = np.where( ionized_frac_interp > 0.97)[0]
  z_ionized = z_interp[ionized.max()]
  
  color = colors[data_name]
  label = labels[data_name]
  ax.plot( z, ionized_frac, c=color, linewidth=line_width, label=label )
  ax.axvline( x=z_ionized, ls='--', c=color )
  text = r'$z_{\mathrm{98\%}}$' + '$ \, = \, {0:.2f}$'.format( z_ionized)
  if chem_type == 'HeII':
    if data_name == 'P19': text_pos_x = 0.155
    if data_name == 'P19m': text_pos_x = 0.225
  if chem_type == 'HI':
    if data_name == 'P19': text_pos_x = 0.165
    if data_name == 'P19m': text_pos_x = 0.205
  text_pos_y = 0.5
  ax.text(text_pos_x, text_pos_y, text,rotation=90, horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=font_size, color=color )

  
if chem_type == 'HeII': ax.set_ylabel( r'$\mathrm{ HeII \,\, Ionization \,\, Fraction }$', fontsize=font_size, color=text_color )
if chem_type == 'HI': ax.set_ylabel( r'$\mathrm{ HI \,\, Ionization \,\, Fraction }$', fontsize=font_size, color=text_color )
ax.set_xlabel( r'$z$', fontsize=font_size, color=text_color )
if chem_type == 'HI': ax.set_xlim( 5, 12 )
if chem_type == 'HeII': ax.set_xlim( 2.5, 5 )
ax.set_ylim( -0.0, 1.02 )


ax.tick_params(axis='both', which='major', direction='in', color=text_color, labelcolor=text_color, labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major  )
ax.tick_params(axis='both', which='minor', direction='in', color=text_color, labelcolor=text_color, labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor  )


leg = ax.legend(loc=1, frameon=False, fontsize=22, prop=prop)
for text in leg.get_texts():
  plt.setp(text, color = text_color)

if black_background: 
  fig.patch.set_facecolor('black') 
  ax.set_facecolor('k')
  [ spine.set_edgecolor(text_color) for spine in list(ax.spines.values()) ]
    


figure_name = output_dir + f'fig_ionized_fraction_{chem_type}.png'
fig.savefig( figure_name, bbox_inches='tight', dpi=300, facecolor=fig.get_facecolor() )
print( f'Saved Figure: {figure_name}' )



