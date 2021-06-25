import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import palettable

import matplotlib
# set some global options
matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'
# hfont = {'fontname':'Helvetica'}
# 


data = np.array([ #n_per_gpu   #n_gpus #T_hydro  
                  [ 128, 8,     48.494, 148.336, 74.485, 75.753 ],
                  [ 128, 64,    48.494, 178.288, 81.616, 77.971 ],
                  [ 128, 512,   49.287, 204.754, 87.956, 81.312 ],
                  [ 128, 1024,  48.972, 217.751, 89.223, 83.201 ],
                  [ 128, 2048,  49.762, 223.772, 91.125, 82.567 ],
                  [ 128, 4096,  50.238, 231.854, 90.333, 85.737 ],
                  [ 128, 8192,  50.396, 239.303, 94.929, 80.032 ],
                  [ 128, 16384, 52.773, 257.211, 94.612, 76.862 ],
                  [ 256, 8,     42.631, 168.463, 50.872, 79.873 ],
                  [ 256, 64,    44.532, 192.235, 54.041, 80.032 ],
                  [ 256, 512,   43.582, 216.957, 59.746, 86.529 ],
                  [ 256, 1024,  43.899, 229.635, 58.637, 88.431 ],
                  [ 256, 2048,  44.057, 237.401, 58.321, 87.797 ],
                  [ 256, 4096,  44.691, 241.363, 59.429, 93.027 ],
                  [ 256, 8192,  45.642, 253.883, 67.512, 85.261 ],
                  [ 256, 16384, 47.544, 264.025, 68.463, 87.322 ]  ]).T


plot_grav_over_logN = False

#size nproc  
n_per_gpu   = data[0]
n_procs     = data[1]
t_hydro     = data[2]
t_grav      = data[3]
t_mpi       = data[4]
t_particles = data[5]
t_total = t_hydro + t_grav + t_mpi + t_particles


indx_128 = n_per_gpu == 128
n_procs_128 = n_procs[indx_128]
t_hydro_128 = t_hydro[indx_128] 
t_mpi_128 = t_mpi[indx_128]
t_grav_128 = t_grav[indx_128]
t_particles_128 = t_particles[indx_128]
t_total_128 = t_total[indx_128]

indx_256 = n_per_gpu == 256
n_procs_256 = n_procs[indx_256]
t_hydro_256 = t_hydro[indx_256] 
t_mpi_256 = t_mpi[indx_256] 
t_grav_256 = t_grav[indx_256]
t_particles_256 = t_particles[indx_256]
t_total_256 = t_total[indx_256]



fig = plt.figure(0)
fig.set_size_inches(7,6)
fig.clf()
ax = plt.gca()

c_hydro = 'C0'
c_mpi = 'C4'
c_grav = 'C3'
c_particles = 'C2'
c_total = 'k'



ms = 4
ax.plot( n_procs_128, t_hydro_128, '--', c=c_hydro,  alpha=0.8, marker='o', markersize=ms)
ax.plot( n_procs_128, t_mpi_128, '--', c=c_mpi, alpha=0.8, marker='o', markersize=ms)
ax.plot( n_procs_128, t_particles_128, '--', c=c_particles,  alpha=0.8, marker='o', markersize=ms)
ax.plot( n_procs_128, t_grav_128, '--', c=c_grav, alpha=0.8, marker='o', markersize=ms)
ax.plot( n_procs_128, t_total_128, '--', c=c_total,  alpha=0.8, marker='o', markersize=ms, label=r'128$^3$ / GPU')

ax.plot( n_procs_256, t_hydro_256, c=c_hydro, alpha=0.6, marker='D', markersize=ms)
ax.plot( n_procs_256, t_mpi_256, c=c_mpi, alpha=0.6, marker='D', markersize=ms)
ax.plot( n_procs_256, t_particles_256, c=c_particles, alpha=0.6, marker='D', markersize=ms)
ax.plot( n_procs_256, t_grav_256, c=c_grav, alpha=0.6, marker='D', markersize=ms)
ax.plot( n_procs_256, t_total_256, c=c_total, alpha=0.6, marker='D', markersize=ms,  label=r'256$^3$ / GPU')


ax.legend( loc=2, frameon=False, fontsize=9)

fs = 8
ax.text(0.05, 0.05, 'Hydro', fontsize=fs, color=c_hydro, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

ax.text(0.05, 0.11, 'MPI comm', fontsize=fs, color=c_mpi, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

ax.text(0.05, 0.17, 'Particles', fontsize=fs, color=c_particles, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)


ax.text(0.05, 0.34, 'Poisson', fontsize=fs, color=c_grav, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
ax.text(0.05, 0.68, 'Total', fontsize=fs, color=c_total, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)



ax.text(0.45, 0.95, 'Cholla Weak Scaling on Summit 2021', fontsize=12, color=c_total, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)


# ax.set_ylim(0, 280)
# ax.set_ylim(0, 300)

ax.set_xlim(6, 20000)
ax.set_ylim(0., 550)
# ax.set_yscale('log')

ax.tick_params(axis='both', which='major', direction='in' )
ax.tick_params(axis='both', which='minor', direction='in' )


fs = 10
ax.set_ylabel( r'Milliseconds / 128$^3$ Cells / GPU', fontsize=fs)
ax.set_xlabel( r'Number of GPUs', fontsize=fs)
ax.set_xscale('log')

output_dir = './'
fileName = output_dir + 'scaling_summit_adiabatic_FOM.png'
fig.savefig(  fileName ,  bbox_inches='tight', dpi=300)

