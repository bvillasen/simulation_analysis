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


data = np.array([                        # dt    hydo     bound  grav_pot  pot_bou  part_dens   part_bd part_ds_  part_adv_1  part_adv_2  total
  [ 256, 8,    512,  512,  512,  7, 25, 5.18826,  98.5283, 16.9028, 161.221, 1.80492, 1.46012,  9.62755, 0.843534, 0.0111008, 0.0298786, 295.617 ],
  [ 256, 64,   1024, 1024, 1024, 7, 25, 5.18076,  98.5481, 40.8246, 340.91 , 4.9009,  1.47477,  9.60053, 1.98144,  0.0128269, 0.0335979, 503.468 ],
  [ 256, 512,  2048, 2048, 2048, 7, 25, 5.19629,  98.8552, 48.7151, 548.141, 5.71922, 1.48909,  9.62788, 2.19488,  0.0178242, 0.043745,  720  ], 
  [ 256, 1024, 4096, 2048, 2048, 7, 25, 5.17821,  98.7652, 77.3116, 597.763, 9.48107, 1.49278,  9.58447, 3.65285,  0.0169182, 0.045414,  803.292],
  [ 128, 8,    256,  256,  256,  7, 25, 0.605288, 13.5077, 4.35181, 21.0332, 0.57868, 0.221691, 2.0584,  0.383415, 0.0104427, 0.0319099, 42.7826 ], 
  [ 128, 64,   512,  512,  512,  7, 25, 0.600977, 14.1521, 9.63328, 43.6517, 1.32642, 0.233574, 2.15286, 0.628576, 0.012331,  0.0352573, 72.4271 ],
  [ 128, 512,  1024, 1024, 1024, 7, 25, 0.588837, 14.3254, 16.5381, 72.1548, 2.32326, 0.242958, 2.28973, 0.713434, 0.01791,   0.0421333, 109.237 ],
  [ 128, 1024, 2048, 1024, 1024, 7, 25, 0.591526, 14.3312, 19.9013, 77.1555, 2.75222, 0.25157,  2.3151,  0.841799, 0.0170517, 0.043087 , 118.2 ],
  [ 256, 8192, 8192, 4096, 4096, 7, 25, 5.18116,  98.8896, 83.2331, 689.307, 10.0834, 1.49766,  9.61277, 3.59948,  0.0193501, 0.0525856, 901.476 ],
  [ 128, 2048, 2048, 2048, 1024, 7, 25, 0.585804, 14.3372, 18.098,  82.8333, 2.31804, 0.25588,  2.34041, 0.684147, 0.0203896, 0.0452805, 121.518 ], 
  [ 256, 2048, 4096, 4096, 2048, 7, 25, 5.18372, 98.8032, 76.4487, 612.993, 9.04003, 1.52037, 9.60188, 3.21218, 0.0198078, 0.048399, 816.871 ],
  [ 256, 4096, 4096, 4096, 4096, 7, 25, 5.13445, 98.8145, 78.0825, 648.206, 9.51975, 1.50611, 9.57257, 3.3093, 0.0261402, 0.0579548, 854.23 ]  
]).T



#size nproc  
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
ax.text(0.05, 0.615, 'Hydro', fontsize=fs, color=c_hydro, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

ax.text(0.05, 0.49, 'MPI comm', fontsize=fs, color=c_mpi, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

ax.text(0.05, 0.14, 'Particles', fontsize=fs, color=c_particles, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

ax.text(0.05, 0.69, 'Poisson', fontsize=fs, color=c_grav, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
ax.text(0.05, 0.78, 'Total', fontsize=fs, color=c_total, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

ax.text(0.45, 0.95, 'Cholla Weak Scaling on Summit 2021', fontsize=12, color=c_total, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)


# ax.set_ylim(0, 280)
# ax.set_ylim(0, 300)

ax.set_xlim(6, 16200)
ax.set_ylim(0.1, 400)
ax.set_yscale('log')

fs = 10
ax.set_ylabel( r'Milliseconds / 128$^3$ Cells / GPU', fontsize=fs)
ax.set_xlabel( r'Number of GPUs', fontsize=fs)
ax.set_xscale('log')

output_dir = '/home/bruno/Desktop/'

fileName = output_dir + 'scaling_summit_adiabatic_2021_log_new.png'
fig.savefig(  fileName ,  bbox_inches='tight', dpi=300)

