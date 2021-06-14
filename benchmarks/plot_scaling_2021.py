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


data = np.array([
  [ 128, 8,   256,  256,  256,  7, 25, 0.61058, 13.5332, 4.37364, 20.8606, 0.57814, 0.23009, 2.70891, 0.386009, 0.0103664, 0.0313950, 43.323  ],
  [ 128, 64,  512,  512,  512,  7, 25, 0.61211, 14.2115, 9.45942, 43.2449, 1.31321, 0.23378, 3.68409, 0.639925, 0.0147057, 0.0329208, 73.4466 ], 
  [ 128, 512, 1024, 1024, 1024, 7, 25, 0.57550, 14.3294, 14.0297, 72.5421, 1.91953, 0.24733, 5.06349, 0.753593, 0.0143433, 0.0409985, 109.516 ],
  [ 256, 8,   512,  512,  512,  7, 25, 5.30103, 98.4441, 16.9867, 164.173, 1.85945, 1.88611, 10.7939, 0.882177, 0.0176144, 0.0294209, 300.374 ], 
  [ 256, 64,  1024, 1024, 1024, 7, 25, 5.47482, 98.5077, 40.2253, 343.83,  4.68699, 3.11683, 16.1255, 1.914600, 0.0217533, 0.0384808, 513.942 ], 
  [ 256, 512, 2048, 2048, 2048, 7, 25, 8.62242, 98.7881, 51.2002, 553.771, 6.14125, 22.6115, 22.4021, 2.403140, 0.0201893, 0.0444412, 766.005 ], 
]).T



#size nproc  dt  hydo  bound  grav_pot  pot_bound  part_dens  part_bound  part_dens_boud  part_adv_1  part_adv_2  total  
n_per_gpu = data[0]
n_procs = data[1]
t_dt = data[7]
t_hydro = data[8]
t_bound = data[9]
t_pot = data[10]
t_pot_bound = data[11]
t_part_dens = data[12]
t_part_bound = data[13]
t_part_dens_bound = data[14]
t_part_adv1 = data[15]
t_part_adv2 = data[16]
t_total = data[17]


t_hydro = t_hydro + t_dt
t_mpi = t_bound + t_pot_bound + t_part_bound + t_part_dens_bound
t_grav = t_pot
t_particles = t_part_dens + t_part_adv1 + t_part_adv2
t_total_1 = t_hydro + t_mpi + t_grav + t_particles


indx_128 = n_per_gpu == 128
n_procs_128 = n_procs[indx_128]
t_hydro_128 = t_hydro[indx_128] 
t_mpi_128 = t_mpi[indx_128]
t_grav_128 = t_grav[indx_128]
t_particles_128 = t_particles[indx_128]
t_total_128 = t_total[indx_128]


indx_256 = n_per_gpu == 256
n_procs_256 = n_procs[indx_256]
t_hydro_256 = t_hydro[indx_256] / 8 
t_mpi_256 = t_mpi[indx_256] / 8
t_grav_256 = t_grav[indx_256]/ 8
t_particles_256 = t_particles[indx_256]/ 8
t_total_256 = t_total[indx_256]/ 8


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
ax.plot( n_procs_128, t_grav_128, '--', c=c_grav, alpha=0.8, marker='o', markersize=ms)
ax.plot( n_procs_128, t_particles_128, '--', c=c_particles,  alpha=0.8, marker='o', markersize=ms)
ax.plot( n_procs_128, t_total_128, '--', c=c_total,  alpha=0.8, marker='o', markersize=ms, label=r'128$^3$ / GPU')


ax.plot( n_procs_256, t_hydro_256, c=c_hydro, alpha=0.6, marker='D', markersize=ms)
ax.plot( n_procs_256, t_mpi_256, c=c_mpi, alpha=0.6, marker='D', markersize=ms)
ax.plot( n_procs_256, t_grav_256, c=c_grav, alpha=0.6, marker='D', markersize=ms)
ax.plot( n_procs_256, t_particles_256, c=c_particles, alpha=0.6, marker='D', markersize=ms)
ax.plot( n_procs_256, t_total_256, c=c_total, alpha=0.6, marker='D', markersize=ms,  label=r'256$^3$ / GPU')

ax.legend( loc=2, frameon=False, fontsize=9)

fs = 8
ax.text(0.05, 0.015, 'Hydro', fontsize=fs, color=c_hydro, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

ax.text(0.05, 0.21, 'MPI comm', fontsize=fs, color=c_mpi, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

ax.text(0.05, 0.18, 'Particles', fontsize=fs, color=c_particles, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

ax.text(0.05, 0.15, 'Poisson', fontsize=fs, color=c_grav, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
ax.text(0.05, 0.32, 'Total', fontsize=fs, color=c_total, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

ax.text(0.4, 0.93, 'Cholla Weak Scaling on Summit 2021', fontsize=12, color=c_total, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)


# ax.set_ylim(0, 280)
ax.set_ylim(0, 200)

fs = 10
ax.set_ylabel( r'Milliseconds / 128$^3$ Cells / GPU', fontsize=fs)
ax.set_xlabel( r'Number of GPUs', fontsize=fs)
ax.set_xscale('log')

output_dir = '/home/bruno/Desktop/'

fileName = output_dir + 'scaling_summit_adiabatic_2021.png'
fig.savefig(  fileName ,  bbox_inches='tight', dpi=300)

