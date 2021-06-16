import os, sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
import h5py as h5
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from colors import *

output_dir = '/home/bruno/Desktop/'

file_name = 'CloudyData_UVB_Puchwein2019_cloudy.h5'
file = h5.File( file_name, 'r' )
rates = file['UVBRates']

z = rates['z'][...]
rates_heating = rates['Photoheating']
heating_HI   = rates_heating['piHI'][...]
heating_HeI  = rates_heating['piHeI'][...]
heating_HeII = rates_heating['piHeII'][...]
rates_ionization = rates['Chemistry']
ionization_HI   = rates_ionization['k24'][...]
ionization_HeI  = rates_ionization['k26'][...]
ionization_HeII = rates_ionization['k25'][...]


scale_He_vals  = [ 0.3,   0.5, 0.8, 1.0 ]
scale_H_vals   = [ 0.6,  0.73, 0.86, 1.0 ]
deltaZ_He_vals = [ -0.1,  0.2, 0.5, 0.8 ]
deltaZ_H_vals  = [ -0.4, -0.2, 0.0, 0.2 ]

import matplotlib.font_manager
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'

prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/', "Helvetica.ttf"), size=16)


font_size = 20
line_width = 2

tick_size_major, tick_size_minor = 6, 4
tick_label_size_major, tick_label_size_minor = 14, 12
tick_width_major, tick_width_minor = 1.5, 1

black_background = True
if black_background:
  text_color = 'white'
  c_0 = yellows[5]
  c_1 = yellows[2]
  c_2 = greens[5]
  c_3 = blues[5]
  c_0 = matter[0]
  c_1 = matter[4]
  c_2 = greens[4]
  c_3 = blues[4]
  colors = [ c_0, c_1, c_2, c_3  ] 


nrows, ncols = 1, 2
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))
plt.subplots_adjust(  wspace=0.15 ) 

chem_type = 'H'
chem_type = 'He'

mod_type = 'scale'
# mod_type = 'delta_z'


if chem_type == 'H':
  ionization_rate = ionization_HI
  heating_rate = heating_HI
  scale_vals = scale_H_vals
  deltaz_vals = deltaZ_H_vals
  label_scale = '\\beta_{\mathrm{H}}'
  label_deltaz = '\\Delta z _{\mathrm{H}}'
  
if chem_type == 'He':
  ionization_rate = ionization_HeII
  heating_rate = heating_HeII
  scale_vals = scale_He_vals
  deltaz_vals = deltaZ_He_vals
  label_scale = '\\beta_{\mathrm{He}}'
  label_deltaz = '\\Delta z _{\mathrm{H}}'
  
  
  
ax = ax_l[0]
if mod_type == 'scale':
  for i,scale_val in enumerate(scale_vals):
    str = '{0} \, = \, {1:.1f}'.format(label_scale, scale_val)
    label = r'${0}$'.format(str)
    ax.plot( z, ionization_rate * scale_val , color=colors[i], lw=line_width, label=label )

if mod_type == 'delta_z':
  for i,deltaz_val in enumerate(deltaz_vals):
    str = '{0} \, = \, {1:.1f}'.format(label_deltaz, deltaz_val)
    label = r'${0}$'.format(str)
    ax.plot( z+deltaz_val, ionization_rate  , color=colors[i], lw=line_width, label=label )
    

if chem_type == 'H': text_rate = r'$\mathrm{ Rate \,\, for \,\, HI }$'
if chem_type == 'He': text_rate = r'$\mathrm{ Rate \,\, for \,\, HeII }$'
ax.text(0.85, 0.95, text_rate, horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=18, color=text_color) 



ax.set_yscale('log')

ax.set_xlabel( r'$z$', fontsize=24, color=text_color )
ax.set_ylabel( r'$\mathrm{Photoionization \,\, Rate \,\,\, [s^{-1}]}$', fontsize=font_size, color=text_color )
ax.tick_params(axis='both', which='major', direction='in', color=text_color, labelcolor=text_color, labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major  )
ax.tick_params(axis='both', which='minor', direction='in', color=text_color, labelcolor=text_color, labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor  )

ax.set_xlim( 0.8, 14 )

leg = ax.legend(loc=3, frameon=False, fontsize=26, prop=prop)
for text in leg.get_texts():
  plt.setp(text, color = text_color)

if black_background: 
  fig.patch.set_facecolor('black') 
  ax.set_facecolor('k')
  [ spine.set_edgecolor(text_color) for spine in list(ax.spines.values()) ]


ax = ax_l[1]
if mod_type == 'scale':
  label_root = label_scale
  for i,scale_val in enumerate(scale_vals):
    str = '{0} \, = \, {1:.1f}'.format(label_scale, scale_val)
    label = r'${0}$'.format(str)
    ax.plot( z, heating_rate * scale_val , color=colors[i], lw=line_width, label=label )

if mod_type == 'delta_z':
  label_root = label_deltaz
  for i,deltaz_val in enumerate(deltaz_vals):
    label = r'${0} \, = \, {1:.1f}$'.format(label_deltaz, deltaz_val)
    ax.plot( z+deltaz_val, heating_rate  , color=colors[i], lw=line_width, label=label )


ax.text(0.85, 0.95, text_rate, horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=18, color=text_color) 


ax.set_yscale('log')

ax.set_xlabel( r'$z$', fontsize=24, color=text_color )
ax.set_ylabel( r'$\mathrm{Photoheating \,\, Rate \,\,\, [eV \, s^{-1}]}$', fontsize=font_size, color=text_color )
ax.tick_params(axis='both', which='major', direction='in', color=text_color, labelcolor=text_color, labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major  )
ax.tick_params(axis='both', which='minor', direction='in', color=text_color, labelcolor=text_color, labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor  )

ax.set_xlim( 0.8, 14 )

leg = ax.legend(loc=3, frameon=False, fontsize=26, prop=prop)
for text in leg.get_texts():
  plt.setp(text, color = text_color)

if black_background: 
  fig.patch.set_facecolor('black') 
  ax.set_facecolor('k')
  [ spine.set_edgecolor(text_color) for spine in list(ax.spines.values()) ]



    

figure_name = output_dir + f'UVB_rates_{chem_type}_{mod_type}.png'
fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
print( f'Saved Figure: {figure_name}' )







