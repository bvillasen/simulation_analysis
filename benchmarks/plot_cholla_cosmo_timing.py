import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import matplotlib
# set some global options
matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'

output_dir = '/home/bruno/Desktop/'

data = np.loadtxt( 'timing_cosmo_data.txt').T

n_mpi, nx, ny, nz, n_omp, n_steps, dt, hydro, bound, grav_pot, pot_bound,  part_dens,  part_bound,  part_dens_boud,  part_adv_1,  part_adv_2,  total = data

time_hydro = hydro
time_mpi = bound + pot_bound + part_dens_boud + part_bound
time_poisson = grav_pot
time_particles = dt + part_dens + part_adv_1 + part_adv_2 
time_total = total

label_size = 14


xlabels = ['Hydro', 'MPI', 'Poisson (Paris)', 'Particles', 'Total']

names = [ 'CPU', 'PARTICLES_GPU', 'PARTICLES_GPU + GRAVITY_GPU', 'PARTICLES_GPU + GRAVITY_GPU + GPU_MPI' ]
n = len(hydro)

times_all, speed_all = [], []
for i in range(n):
  times = np.array([ time_hydro[i], time_mpi[i], time_poisson[i], time_particles[i], time_total[i]])
  times_all.append(times)
  speed = times_all[0] / times  
  speed_all.append( speed ) 
  
speed_all = np.array( speed_all )
  

x = np.arange(len(xlabels))  # the label locations
width = 0.9  # the width of the bars


ncols, nrows = 1, 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12*ncols,8*nrows))
 
for i in range(n):
  times = times_all[i]
  rects = ax.bar(x - width/2 + width*(i+0.5)/4 , times, width/4, label=names[i])
  labels = [ f'{speed:.1f}x' if (speed < 10 and speed > 1) or speed < 0.95 else f'{speed:.0f}x' for speed in speed_all[i] ]
  print(labels)
  ax.bar_label(rects, labels=labels, padding=3)



ax.tick_params(axis='both', which='major', labelsize=14, size=5, width=1.5, direction='in')


ax.set_title('CHOLLA Cosmological Simulation Timing', fontsize=14 )
ax.set_xticks(x)
ax.set_xticklabels(xlabels, fontsize=label_size)
ax.set_ylabel('Time per Iteration  [ms]', fontsize=label_size)
ax.legend( frameon=False, fontsize=12 )
# 
# # ax.bar_label(rects2, padding=3)
# 
# 
# 
figure_name = output_dir + f'cholla_cosmo_timing.png'
fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
print( f'Saved Figure: {figure_name}' )