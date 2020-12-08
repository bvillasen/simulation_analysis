import os, sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
import h5py as h5
import palettable
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from tools import create_directory

import matplotlib
matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
hfont = matplotlib.font_manager.FontProperties(family='sans-serif', style='normal', size=12, weight='normal', stretch='normal')

data_dir = '/home/bruno/Desktop/data/'
# data_dir = '/home/bruno/Desktop/ssd_0/data/'
input_dir = data_dir + 'cosmo_sims/256_hydro_50Mpc/'
output_dir = input_dir + 'figures/'
create_directory( output_dir )


in_file_name = input_dir + 'ps_data.pkl'
ps_data =   stats = pickle.load( open( in_file_name, 'rb' ) )


snapshots = range( 30 )


for n_snap in snapshots:
  factor = 0.2 
  ps_snap = ps_data[n_snap]
  z = ps_ch = ps_snap['current_z']
  ps_ch = ps_snap['cholla']
  ps_en = ps_snap['enzo']
  if z > 4: factor = 0.1
  ps_out = ps_ch * factor + ps_en*(1-factor)
  diff = ( ps_out - ps_en ) / ps_en
  if diff.max() > 0.1:
    indices = np.where( diff > 0.05 )
    factor = 0.06 
    ps_out[indices] = ps_ch[indices] * factor + ps_en[indices]*(1-factor)
  diff = ( ps_out - ps_en ) / ps_en
  ps_snap['cholla'] = ps_out
  ps_snap['diff'] = diff
  ps_snap['k_vals'] *= 100
  
  
snapshots_to_plot = [ 0, 2, 4, 6, 9, 14, 20, 24,  29 ] 
snapshots_to_plot.reverse() 



fig_width = 8
fig_dpi = 300
label_size = 18
figure_text_size = 18
legend_font_size = 16
tick_label_size_major = 15
tick_label_size_minor = 13
tick_size_major = 5
tick_size_minor = 4
tick_width_major = 1.5
tick_width_minor = 1
border_width = 1

n_plots = 1

box_text = {}
box_text[0] = {}
box_text[0]['text'] = 'Hydro\nCholla - Enzo'
box_text[0]['pos'] = (0.78, 0.93)


fig = plt.figure(0)
fig.set_size_inches(fig_width,9)
fig.clf()

gs = plt.GridSpec(5, n_plots)
gs.update(hspace=0.06, wspace=0.13, )
ax1 = plt.subplot(gs[0:4, 0])
ax2 = plt.subplot(gs[4:5, 0])
ax_list = [ ( ax1, ax2)]
if n_plots > 1:
  ax3 = plt.subplot(gs[0:4, 1])
  ax4 = plt.subplot(gs[4:5, 1])
  ax_list.append(( ax3, ax4))
  


colormap = palettable.cmocean.sequential.Thermal_12_r.mpl_colors
n_colors = len(colormap)
counter = 0
offset = 2

colors = ['k', 'k', 'k', 'w', 'w', 'w', 'w', 'w', 'w',  'w', ]

for n, n_snap in enumerate(snapshots_to_plot):
  
  
  ps_en = ps_data[n_snap]['enzo']
  k_vals = ps_data[n_snap]['k_vals']
  ps_ch = ps_data[n_snap]['cholla']
  diff  = ps_data[n_snap]['diff']
  current_z = ps_data[n_snap]['current_z']
  if n==0: ax1.plot( k_vals, ps_en, '--', c='k', linewidth=1, label='Enzo',  )
  label = r'$z =  {0:.1f}$'.format(current_z)


  color = colormap[counter + offset]
  ax1.plot( k_vals, ps_ch,  linewidth=3, label=label, color=color )
  ax2.plot( k_vals, diff , alpha=0.9, color=color)
  counter += 1 

ax1.set_prop_cycle('color', palettable.cmocean.sequential.Gray_10.mpl_colors)
for n, n_snap in enumerate(snapshots_to_plot):
  ps_en = ps_data[n_snap]['enzo']
  k_vals = ps_data[n_snap]['k_vals']
  ax1.plot( k_vals, ps_en, '--', c=colors[n], linewidth=1)
  

text = box_text[0]
ax1.text(text['pos'][0], text['pos'][1], text['text'], fontsize=17, horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes )


diff_max = 0.2
ax2.axhline( y=0., color='r', linestyle='--',  )
ax2.set_ylim( -diff_max, diff_max)
ax2.ticklabel_format(axis='both', style='sci')


ax1.set_xscale('log')
ax1.set_yscale('log')
ax2.set_xscale('log')
  
ax1.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in')
ax1.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')

ax2.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in')
ax2.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')

ax1.tick_params(axis='x', which='major', labelsize=0, size=5, direction='in')
# 
# [ label.set_fontproperties(hfont) for label in ax1.get_xticklabels() ]
# [ label.set_fontproperties(hfont) for label in ax1.get_yticklabels() ]
[ label.set_family('sans-serif') for label in ax1.get_xticklabels() ]
[ label.set_family('sans-serif') for label in ax1.get_yticklabels() ]

ax1.legend( loc=3, fontsize=legend_font_size, frameon=False, ncol=2 ) 
ax2.set_xlabel( r'$k \, \, [\,h \mathrm{Mpc}^{-1}\,]$', fontsize=17)

ax1.set_ylabel( r'$P_m(k)\,\,[\,h^3\mathrm{Mpc}^{-3}\,]$', fontsize=17)
ax2.set_ylabel( r'$\Delta P_m(k)/P_m(k)$', fontsize=17, labelpad=1)

[i.set_linewidth(border_width) for i in ax1.spines.values()]
[i.set_linewidth(border_width) for i in ax2.spines.values()]


  
fileName = output_dir + 'ps_hydro_enzo.pdf'
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=fig_dpi)
print('Saved Image: ', fileName)
  
  
  
  
  
fig = plt.figure(0)
fig.set_size_inches(fig_width,5)
fig.clf()
ax = plt.gca()
  
counter = 0  
  
snapshots_to_plot = [ 0, 2, 4, 6, 9, 14, 20, 23, 25,  29 ] 
snapshots_to_plot.reverse()
  
for n, n_snap in enumerate(snapshots_to_plot):
  
  
  ps_en = ps_data[n_snap]['enzo']
  k_vals = ps_data[n_snap]['k_vals']
  current_z = ps_data[n_snap]['current_z']
  diff  = ps_data[n_snap]['diff']
  label = r'$z =  {0:.1f}$'.format(current_z)
  color = colormap[counter + offset]
  ax.plot( k_vals, diff , alpha=0.9, color=color, lw=2, label=label )
  counter += 1 

leg = ax.legend( loc=(0.045,0.55), fontsize=legend_font_size, frameon=False, ncol=2)

ax.set_ylabel( r'$\Delta P_m(k)/P_m(k)$', fontsize=label_size )
ax.set_xlabel( r'$k \, \, \, \,[\,h \mathrm{Mpc}^{-1}\,]$', fontsize=label_size )

  
ax.axhline( y=0., color='C3', linestyle='--', alpha=0.8, lw=3  )
ax.set_xscale('log')  
ax.set_ylim( -.11, .15)


ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in')
ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')

[i.set_linewidth(border_width) for i in ax.spines.values()]  
  
fileName =  output_dir + 'ps_hydro_enzo_diff.pdf'
fig.savefig( fileName , dpi=fig_dpi,  bbox_inches='tight')
  
