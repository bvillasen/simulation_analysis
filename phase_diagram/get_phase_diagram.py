import os, sys
from os import listdir
from os.path import isfile, join
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import ImageGrid
import pickle
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from load_data import load_snapshot_data_distributed
from tools import *
from phase_diagram_functions import *
from turbo_cmap import *

use_mpi = True
if use_mpi :
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  nprocs = comm.Get_size()
else:
  rank = 0
  nprocs = 1
  

show_progess = False
if rank == 0: show_progess = True

data_dir = '/raid/bruno/data/'

n_points = 256
nx = n_points
ny = n_points
nz = n_points
ncells = nx * ny * nz

input_grackle = data_dir + 'cosmo_sims/256_hydro_50Mpc/data_grackle/'
input_chem    = data_dir + 'cosmo_sims/256_hydro_50Mpc/output_files/'
output_dir = data_dir + 'cosmo_sims/256_hydro_50Mpc/figures/phase_diagram/'
if rank == 0: create_directory( output_dir )

Lbox = 50000.0 #kpc/h
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ n_points, n_points, n_points ]
precision = np.float64



input_dirs = {'grackle': input_grackle, 'chem':input_chem }

types = [ 'grackle', 'chem' ]

pd_data = {}

# n_snaps = 170
# snaps_all = range( n_snaps )
snaps_all = np.arange( 80, 109, 1 )
snaps_proc = np.array(split_indices( snaps_all, rank, nprocs ))
snaps_proc = snaps_all[snaps_proc]

if len(snaps_proc) == 0: exit()

print( f'pID: {rank}  snaps:{snaps_proc} ' )

for n_snap in snaps_proc:
  

  for type in types:

    input_dir = input_dirs[type]
    fields = [ 'density', 'temperature'  ]
    data_gas = load_snapshot_data_distributed( 'hydro', fields, n_snap, input_dir, box_size, grid_size,  precision, show_progess=show_progess )
    current_z = data_gas['Current_z']
    density = data_gas['density']
    temperature = data_gas['temperature']
    dens_mean = np.mean( density )

    #Get Overdensity
    density = density / dens_mean

    #Get Bin Egdes for the histogram
    dens_start, dens_end = -2, 4
    temp_start, temp_end = -2, 8
    nbins = 800
    bins_dens = np.logspace( dens_start, dens_end, nbins, base=10 )
    bins_temp = np.logspace( temp_start, temp_end, nbins, base=10 )

    print(" Generating Phase Diagram ,   n_bins:{0}".format(nbins))
    centers_dens, centers_temp, phase = get_phase_diagram_bins( density, temperature, bins_dens, bins_temp  )
    dens_points, temp_points = np.meshgrid( centers_dens, centers_temp )
    pd_data['data'] = phase
    pd_data['dens_min'], pd_data['dens_max'], pd_data['n_dens'] = 10**centers_dens.min(), 10**centers_dens.max(), len( centers_dens)
    pd_data['temp_min'], pd_data['temp_max'], pd_data['n_temp'] = 10**centers_temp.min(), 10**centers_temp.max(), len( centers_temp)  
    temp_points = temp_points.flatten()
    dens_points = dens_points.flatten()
    phase_1D = phase.flatten() 
    indices = np.where(phase_1D > 0 )
    phase_1D = phase_1D[indices]
    dens_points = dens_points[indices]
    temp_points = temp_points[indices]
    pd_data[type] = { 'dens_points':dens_points, 'temp_points':temp_points, 'phase_1D':phase_1D, 'phase_2D':phase,  }



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

  n_data = 2
  n_rows = 1
  n_cols = n_data

  title_all = [ 'Grackle', 'Chem' ]

  import matplotlib
  matplotlib.rcParams['mathtext.fontset'] = 'cm'
  matplotlib.rcParams['mathtext.rm'] = 'serif'

  # Set up figure and image grid
  fig = plt.figure(0, figsize=(fig_width*n_cols,10*n_rows),  )
  grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                   nrows_ncols=(n_rows,n_cols),
                   axes_pad=0.2,
                   share_all=True,
                   cbar_location="right",
                   cbar_mode="single",
                   cbar_size="5%",
                   cbar_pad=0.1,
                   )

  colormap =  'turbo'
  alpha = 0.6

  x_min, x_max = dens_start, dens_end
  y_min, y_max = temp_start, temp_end

  for i, type in enumerate( types ):


    ax = grid[i]

    dens_points = pd_data[type]['dens_points']
    temp_points = pd_data[type]['temp_points']
    phase_1D = np.log10(pd_data[type]['phase_1D'])    
    phase_2D = np.log10(pd_data[type]['phase_2D'])    
    v_min, v_max = min( np.log10(pd_data['grackle']['phase_1D']).min(), np.log10(pd_data['chem']['phase_1D']).min() ), max( np.log10(pd_data['grackle']['phase_1D']).max(), np.log10(pd_data['chem']['phase_1D']).max() )    

    im = ax.scatter( dens_points, temp_points, c=phase_1D, s=0.1, vmin=v_min, vmax=v_max, alpha=alpha, cmap=colormap  )
    im = ax.imshow( phase_2D[::-1], extent=[ dens_start, dens_end, temp_start, temp_end ] )



    cb = ax.cax.colorbar(im,   )
    cb.ax.tick_params(labelsize=tick_label_size_major, size=tick_size_major, color=text_color, width=tick_width_major, length=tick_size_major, labelcolor=text_color, direction='in' )
    ax.cax.toggle_label(True)
    [sp.set_linewidth(border_width) for sp in cb.ax.spines.values()]

    ax.set_aspect( 0.6 )
    # ax.set_aspect( 1.0)

    font = {'fontname': 'Helvetica',
        'color':  text_color,
        'weight': 'normal',
        'size': label_size,
        'ha':'center'
        }
    cb.set_label_text( r'$\log_{10}  \,\, P\,(\Delta, T\,) $', fontdict=font )
    ax.set_ylabel(r'$\log_{10} \, T \,\,[\,\mathrm{K}\,]$', fontsize=label_size , color=text_color)
    ax.set_xlabel(r'$\log_{10} \, \Delta$ ', fontsize=label_size , color=text_color )


    text  = r'$z = {0:.2f}$'.format( current_z ) 
    if i == 0: ax.text(0.05, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size, color=text_color)

    text = title_all[i]
    ax.text(0.95, 0.95, text, horizontalalignment='right',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size, color=text_color)


    ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in')
    ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')
    [sp.set_linewidth(border_width) for sp in ax.spines.values()]
    ax.set_xlim( x_min, x_max )
    ax.set_ylim( y_min, y_max )

  out_fileName = output_dir + f'phase_diagram_{n_snap}.png'
  fig.savefig( out_fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
  print(( 'Saved Image: ' + out_fileName ))


