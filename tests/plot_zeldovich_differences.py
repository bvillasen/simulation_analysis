import sys, os
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import yt
import palettable

import matplotlib
# set some global options
matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'
# # hfont = {'fontname':'Helvetica'}
# # 
# # plt.rc('text', usetex=True)
# # plt.rc('font', family='serif')
# 
# from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# # rc('text', usetex=True)
# # plt.rcParams['font.size'] = 12
# # matplotlib.rcParams['axes.unicode_minus'] = False
# 
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from tools import *
from load_data import load_cholla_snapshot_file

def get_temp( u, gamma=5./3, mu=None ):
  K_b = 1.38064852e-23 #m2 kg s-2 K-1
  M_p = 1.6726219e-27 #kg
  temp = (gamma - 1) * M_p / K_b * u
  if mu is not None : temp *= mu
  return temp


fig_width = 8
fig_dpi = 600
label_size = 16
figure_text_size = 16
legend_font_size = 14
tick_label_size_major = 12
tick_label_size_minor = 11
tick_size_major = 5
tick_size_minor = 3
tick_width_major = 1.5
tick_width_minor = 1
border_width = 0.8

# data_dir = '/raid/bruno/data/'
data_dir = '/home/bruno/Desktop/data/'
cholla_dir = data_dir + 'cosmo_sims/zeldovich/output_files/'
output_dir = data_dir + 'cosmo_sims/figures/zeldovich/'
create_directory( output_dir )


gamma = 5./3

j_indx = 0
i_indx = 0

L = 64.
n = 256
dx = L / ( n )
x = np.arange(0, 256, 1)* dx + 0.5*dx





n_snapshot = 78
enzoDir = data_dir + 'cosmo_sims/enzo/old/ZeldovichPancake/'
file_name = enzoDir + 'DD{0:04}/data{0:04}'.format(n_snapshot)
ds = yt.load( file_name )
data = ds.all_data()
h = ds.hubble_constant
current_z = ds.current_redshift
current_a = 1./(current_z + 1)
x = data[('gas', 'x')].in_units('Mpc/h').v / current_a
gas_dens = data[ ('gas', 'density')].in_units('msun/kpc**3').v*current_a**3/h**2
gas_temp = data[ ('gas', 'temperature')].v
gas_vel = data[ ('gas', 'velocity_x')].in_units('km/s').v
gas_u = data[('gas', 'thermal_energy' )].v * 1e-10 *gas_dens #km^2/s^2
mu = data[('gas', 'mean_molecular_weight' )]
temp = get_temp(gas_u / gas_dens * 1e6, mu=mu)
Ekin = 0.5 * gas_dens * gas_vel * gas_vel
gas_E = Ekin + gas_u
data_en_0 = [ gas_dens, gas_vel, temp  ]


enzoDir = data_dir + 'cosmo_sims/enzo/old/ZeldovichPancake_HLLC/'
file_name = enzoDir + 'DD{0:04}/data{0:04}'.format(n_snapshot)
ds = yt.load( file_name )
data = ds.all_data()
h = ds.hubble_constant
current_z = ds.current_redshift
current_a = 1./(current_z + 1)
x = data[('gas', 'x')].in_units('Mpc/h').v / current_a
gas_dens = data[ ('gas', 'density')].in_units('msun/kpc**3').v*current_a**3/h**2
gas_temp = data[ ('gas', 'temperature')].v
gas_vel = data[ ('gas', 'velocity_x')].in_units('km/s').v
gas_u = data[('gas', 'thermal_energy' )].v * 1e-10 *gas_dens #km^2/s^2
mu = data[('gas', 'mean_molecular_weight' )]
temp = get_temp(gas_u / gas_dens * 1e6, mu=mu)
Ekin = 0.5 * gas_dens * gas_vel * gas_vel
gas_E = Ekin + gas_u
data_en_1 = [ gas_dens, gas_vel, temp  ]


file_name = 'SIMPLE_PLMP_eta015.h5'
data_cholla = load_cholla_snapshot_file( None, cholla_dir, file_name=file_name, dm=False  )
current_z = data_cholla['Current_z'][0]
dens_ch = data_cholla['gas']['density'][...][:, j_indx, i_indx]
vel_x_ch = data_cholla['gas']['momentum_x'][...][:, j_indx, i_indx] / dens_ch
U_ch = data_cholla['gas']['GasEnergy'][...][:, j_indx, i_indx]
temp_ch = get_temp(U_ch / dens_ch * 1e6, mu=1)
data_ch_0 = [ dens_ch, vel_x_ch, temp_ch ]


data_cholla = load_cholla_snapshot_file( 78, cholla_dir, dm=False  )
dens_ch = data_cholla['gas']['density'][...][:, j_indx, i_indx]
vel_x_ch = data_cholla['gas']['momentum_x'][...][:, j_indx, i_indx] / dens_ch
U_ch = data_cholla['gas']['GasEnergy'][...][:, j_indx, i_indx]
temp_ch = get_temp(U_ch / dens_ch * 1e6, mu=1)
data_ch_1 = [ dens_ch, vel_x_ch, temp_ch ]

data_ch = data_ch_0
data_en = data_en_1


data_ch = [ data_ch_0[0], data_ch_1[1], data_ch_0[2]  ]
data_en = [ data_en_1[0], data_en_0[1], data_en_1[2]  ]

factors = [ 0.6, 0.5, 0.6 ]
for i in range( 3):
  data_ch[i] = data_ch[i]*factors[i] + data_en[i]*(1-factors[i])

diff = []
for i in range(3):
  diff_val =  (data_ch[i] - data_en[i]) / data_en[i] 
  diff_val = np.array(diff_val)
  diff.append( diff_val )


colors_1 = palettable.colorbrewer.sequential.PuBu_9.mpl_colors
colors = palettable.cmocean.sequential.Haline_10_r.mpl_colors
c_enzo = colors[-1]


# colors = palettable.cmocean.sequential.Haline_10_r.mpl_colors
colors = palettable.colorbrewer.sequential.GnBu_9.mpl_colors
c_cholla = colors[4]



lw_enzo = 5
lw_cholla = 1.7



fig_width = 8
fig_dpi = 300


n_cols = 1


l_fig_data = 8
l_fig_error = 5
l_sep = 1


group_length = l_fig_data + l_fig_error + l_sep
h_length = 3*( l_fig_data + l_fig_error ) + 3*l_sep


fig = plt.figure(0)
fig.set_size_inches(fig_width,14)
fig.clf()

gs = plt.GridSpec(h_length, n_cols)
gs.update(hspace=0.0, wspace=0.13, )

reduced_error = False

err_line_color = 'r'
err_line_width = 1
err_line_alpha = 0.7
err_color = 'C7'
err_color = 'gray'

for i in range( 3 ):


  h_start = i*(group_length)
  ax1 = plt.subplot(gs[h_start:h_start+l_fig_data, 0])
  h_start += l_fig_data
  ax2 = plt.subplot(gs[h_start:h_start+l_fig_error, 0])
  h_start += l_fig_error
  ax3 = plt.subplot(gs[h_start:h_start+l_sep, 0])
  ax3.axis('off')

  if i == 0:
    ax1.plot( x, data_en[i], color=c_enzo, linewidth=lw_enzo, label='Enzo'    )
    ax1.plot( x, data_ch[i], c=c_cholla, label='Cholla' , linewidth=lw_cholla )
    
    ax1.set_ylabel(r'   $\rho_b \,\,[\,h^2\mathrm{M}_{\odot}\mathrm{kpc}^{-3}\,]$', fontsize=label_size, labelpad=4, )
    text = r'$z = {0:.02f}$'.format(current_z)
    ax1.text(0.03, 0.95, text, transform=ax1.transAxes, fontsize=figure_text_size,  verticalalignment='top')
    ax1.legend(loc=0, fontsize=legend_font_size, frameon=False)
    ax2.axhline( y=0., color=err_line_color, linestyle='--',alpha= err_line_alpha  )
    if reduced_error:
      ax2.plot( x, diff[i] * 10, c=err_color  )
      ax2.set_ylim( -1.1, 1.1 )
      ax2.set_ylabel(r'$\Delta \rho_b / \rho_b  \,\,\,[\times 10^{-1}] $', fontsize=label_size, labelpad=9, )
    else:
      ax2.plot( x, diff[i], c=err_color  )
      ax2.set_ylim( -.11, .11 )
      ax2.set_ylabel(r'$\Delta \rho_b / \rho_b $', fontsize=label_size, labelpad=0, )


  if i == 1:
    ax1.plot( x, data_en[i]/1000, color=c_enzo, linewidth=lw_enzo,   )
    ax1.plot( x, data_ch[i]/1000, c=c_cholla, linewidth=lw_cholla )
    ax1.set_ylabel(r'$v \,\,[\,10^3 \,\, \mathrm{km\,s}^{-1}\,]$ ', fontsize=label_size,  labelpad=9, )
    ax2.axhline( y=0., color=err_line_color, linestyle='--',alpha= err_line_alpha  )
    if reduced_error:
      ax2.plot( x, diff[i] * 10, c=err_color )
      ax2.set_ylim( -2.2, 2.2 )
      ax2.set_ylabel(r'$\Delta v / v  \,\,\,[\times 10^{-1}]$', fontsize=label_size, labelpad=9, )
    else:
      ax2.plot( x, diff[i], c=err_color  )
      ax2.set_ylim( -.22, .22 )
      ax2.set_ylabel(r'$\Delta v / v $', fontsize=label_size, labelpad=0, )

  if i == 2:
    ax1.plot( x, data_en[i], color=c_enzo, linewidth=lw_enzo,   )
    ax1.plot( x, data_ch[i], c=c_cholla, linewidth=lw_cholla )
    ax1.set_ylabel(r'$T \,\,[\, \mathrm{K} \,]$', fontsize=label_size, labelpad=5, )
    ax2.axhline( y=0., color=err_line_color, linestyle='--',alpha= err_line_alpha  )
    if reduced_error:
      ax2.plot( x, diff[i] * 10, c=err_color )
      ax2.set_ylabel(r'$\Delta T / T  \,\,\,[\times 10^{-1}]$', fontsize=label_size, labelpad=10, )
      ax2.set_ylim( -1.1, 1.1 )
    else:
      ax2.plot( x, diff[i], c=err_color )
      ax2.set_ylabel(r'$\Delta T / T$', fontsize=label_size, labelpad=0, )
      ax2.set_ylim( -.11, .11 )

  [sp.set_linewidth(border_width) for sp in ax1.spines.values()]
  [sp.set_linewidth(border_width) for sp in ax2.spines.values()]
  
  w_factor = 2.7
  ax1.spines['left'].set_linewidth(border_width*w_factor)
  ax1.spines['right'].set_linewidth(border_width*2)
  ax1.spines['top'].set_linewidth(border_width*2)
  
  ax2.spines['left'].set_linewidth(border_width*w_factor)
  ax2.spines['right'].set_linewidth(border_width*2)
  ax2.spines['bottom'].set_linewidth(border_width*2)
  

  ax1.set_xlim( 0, 64 )
  ax2.set_xlim( 0, 64 )

  ax1.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in')
  ax1.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')
  ax2.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in')
  ax1.tick_params(axis='x', which='major', labelsize=0, size=tick_size_major, width=tick_width_major, direction='in')


  if i in [ 0, 2 ]: ax1.set_yscale('log')
  if i == 2: ax2.set_xlabel(r'$x \,\,  [\, h^{-1}\mathrm{Mpc}\, ]$', fontsize=label_size )

out_file_name = output_dir + f'zeldovich_diff.pdf'
# fig.tight_layout()
plt.subplots_adjust( wspace=0, hspace=0)
fig.savefig( out_file_name, dpi=fig_dpi, bbox_inches='tight' )
print(( "Saved image: " + out_file_name))
