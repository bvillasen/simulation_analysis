import os, sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#Extend path to inclide local modules
root_dir = os.path.dirname(os.getcwd())
sub_directories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(sub_directories)
from tools import *



data_dir = '/raid/bruno/data/'
input_dir_base = data_dir + 'cosmo_sims/ics/enzo/wdm/256_hydro_50Mpc_'
output_dir = data_dir + 'cosmo_sims/ics/enzo/wdm/figures/'
create_directory( output_dir )


sim_names = [ 'cdm', 'wdm_m0.5kev', 'wdm_m1.0kev', 'wdm_m3.0kev' ]
labels = ['CDM', 'WDM m = 0.5 keV', 'WDM m = 1.0 keV', 'WDM m = 3.0 keV' ]

data_all = {}
for i,sim_name in enumerate(sim_names):
  input_dir = input_dir_base + f'{sim_name}/'
  in_file_name = input_dir + 'input_powerspec.txt'
  data = np.loadtxt( in_file_name )
  k, P_cdm, P_vcdm, P_total = data.T
  data_all[i] = { 'k':k, 'P_cdm':P_cdm, 'P_vcdm':P_vcdm, 'P_total':P_total, 'label':labels[i] }
  

nrows, ncols = 1, 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))

matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'
prop = matplotlib.font_manager.FontProperties( fname=os.path.join('/home/bruno/fonts/Helvetica', "Helvetica.ttf"), size=12)
font_size = 18
label_size = 16

type = 'P_vcdm'

for i in data_all:
  data = data_all[i]
  k = data['k']
  p = data[type]
  lambda_ps = 1/(2*np.pi**3) * k**3 * p* 10e6
  label = data['label']
  ax.plot( k, lambda_ps, label=label )

ax.set_xlim( 1e-3, 100 )
ax.set_ylim( 1e-4, 300 )

ax.set_xscale('log')
ax.set_yscale('log')
ax.tick_params(axis='both', which='major', direction='in', labelsize=label_size )
ax.tick_params(axis='both', which='minor', direction='in' )
ax.set_ylabel( r'$ \frac{1}{2\pi^3} k^3 P(k)  $', fontsize=font_size  )
ax.set_xlabel(  r'$ k   \,\,\,  [h\,\mathrm{Mpc}^{-1}] $', fontsize=font_size )
leg = ax.legend(loc=2, frameon=False, fontsize=font_size, prop=prop)


figure_name = output_dir + f'fig_wdm_{type}.png'
fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
print( f'Saved Figure: {figure_name}' )