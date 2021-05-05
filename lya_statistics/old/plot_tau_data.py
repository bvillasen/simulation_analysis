import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import palettable
import pylab

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from data_optical_depth import *


import matplotlib
# set some global options
matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'

fig_width = 8
fig_dpi = 300
label_size = 20
figure_text_size = 18
legend_font_size = 16
tick_label_size_major = 15
tick_label_size_minor = 13
tick_size_major = 5
tick_size_minor = 3
tick_width_major = 1.5
tick_width_minor = 1
border_width = 1





c_0 = pylab.cm.viridis(.7)
c_1 = pylab.cm.cool(.3)
c_10 = pylab.cm.cool(.3)

# c_2 = 'C1'
# c_3 = 'C9'
# c_3 = purples[-1]
# c_4 = yellows[3]
c_2 = pylab.cm.inferno(.75)
c_3 = pylab.cm.viridis(.7)
c_3 = pylab.cm.hsv(.5)


dataDir = '/home/bruno/Desktop/ssd_0/data/'
output_dir = dataDir + 'cosmo_sims/sim_grid/figures/'
create_directory( output_dir )


nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_width,6*nrows))
fs = 22

text_color ='black'



c_1 = pylab.cm.viridis(.3)
c_2 = c_3
c_3 = 'C3'
c_4 = 'C9'
c_5 = 'C1'

#Add data Becker
data_tau = data_optical_depth_Becker_2013
z = data_tau['z']
tau = data_tau['tau']
tau_p = data_tau['tau_sigma_p']
tau_m = data_tau['tau_sigma_m']
tau_error = np.array([ tau_p, tau_m])
ax.errorbar( z, tau, yerr=tau_error, fmt='o', c=c_1, label=data_tau['name'] )

#Add data Boera
data_tau = data_optical_depth_Boera_2019 
z = data_tau['z']
tau = data_tau['tau']
tau_p = data_tau['tau_sigma_p']
tau_m = data_tau['tau_sigma_m']
tau_error = np.array([ tau_p, tau_m])
ax.errorbar( z, tau, yerr=tau_error, fmt='o', c=c_2, label=data_tau['name'] )



#Add data Bosman
data_tau = data_optical_depth_Bosman_2018 
z = data_tau['z']
tau = data_tau['tau']
tau_p = data_tau['tau_sigma_p']
tau_m = data_tau['tau_sigma_m']
tau_error = np.array([ tau_p, tau_m])
ax.errorbar( z, tau, yerr=tau_error, fmt='o', c=c_3, label=data_tau['name'] )

data_tau = data_optical_depth_Keating_2020 
z = data_tau['z']
tau = data_tau['tau']
tau_p = data_tau['tau_sigma_p']
tau_m = data_tau['tau_sigma_m']
tau_error = np.array([ tau_p, tau_m])
ax.errorbar( z, tau, yerr=tau_error, fmt='o', c=c_4, label=data_tau['name'] )

data_tau = data_optical_depth_Jiani
z = data_tau['z']
tau = data_tau['tau']
tau_p = data_tau['tau_sigma_p']
tau_m = data_tau['tau_sigma_m']
tau_error = np.array([ tau_p, tau_m])
ax.errorbar( z, tau, yerr=tau_error, fmt='o', c=c_5, label=data_tau['name'] )




ax.tick_params(color=text_color, labelcolor=text_color, labelsize=15,  length=8)
ax.tick_params(which='minor', color=text_color, labelcolor=text_color, labelsize=15,  length=4)
for spine in list(ax.spines.values()):
    spine.set_edgecolor(text_color)
    


leg = ax.legend( loc=2, fontsize=legend_font_size, frameon=False)
for text in leg.get_texts():
    plt.setp(text, color = text_color)
# ax.set_yscale('log')

ax.set_ylabel( r'$\tau_{eff} $', fontsize=label_size, color= text_color  )
ax.set_xlabel(r'$z$', fontsize=label_size, color= text_color )

ax.set_yscale('log')
ax.set_xlim(2, 6.05)
ax.set_ylim(.1, 10)

[sp.set_linewidth(border_width) for sp in ax.spines.values()]

ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in' )
ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')



# ax.grid(True, which="both",)

fileName = output_dir + 'optical_depth_log.png'


fig.savefig( fileName ,  pad_inches=0.1, facecolor=fig.get_facecolor(), bbox_inches='tight', dpi=300)
print('Saved Image: ', fileName)




