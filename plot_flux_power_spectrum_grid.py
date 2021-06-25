import os, sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
import h5py as h5
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from plot_flux_power_spectrum import plot_power_spectrum_grid

data_dir = '/raid/bruno/data/'
output_dir = data_dir + 'cosmo_sims/figures/'
create_directory( output_dir )

ps_data_dir = 'lya_statistics/data/'

file_name = ps_data_dir + 'data_power_spectrum_simulations/data_HM12_P19.h5'
file = h5.File( file_name, 'r' )
z_vals = file['z_vals'][...]
data_ps = file['pchw18']

ps_data = {'z_vals':z_vals, 'label':'P19' }
for indx in data_ps:
  ps_data[int(indx)] = { 'k_vals': data_ps[indx]['kvals'][...], 'ps_mean':data_ps[indx]['delta_power'][...]/data_ps[indx]['kvals'][...]*np.pi}


# plot_power_spectrum_grid( ps_data_dir, output_dir, ps_data=ps_data, scales='large_middle', sim_data_sets=None, system='Shamrock', black_background=True  )
# plot_power_spectrum_grid( ps_data_dir, output_dir, ps_data=ps_data, scales='large_reduced', sim_data_sets=None, system='Shamrock', black_background=True  )
# plot_power_spectrum_grid( ps_data_dir, output_dir, ps_data=ps_data, scales='small_reduced', sim_data_sets=None, system='Shamrock', black_background=True  )


root_dir = data_dir + 'cosmo_sims/sim_grid/1024_P19m_np4_nsim400/'
mcmc_dir = root_dir + 'fit_mcmc/'
data_name = 'fit_results_P(k)+tau_HeII_Boss_Irsic_Boera'
input_dir = mcmc_dir + f'{data_name}/observable_samples/' 
file_name = input_dir + 'samples_power_spectrum.pkl'

samples_ps = Load_Pickle_Directory( file_name )

z_vals = np.array([ samples_ps[i]['z'] for i in samples_ps  ])
ps_data = { 'label': 'This Work (Modified P19)', 'z_vals':z_vals }
for indx in samples_ps:
  ps_data[indx] = { 'k_vals':samples_ps[indx]['k_vals'], 'ps_mean':samples_ps[indx]['Highest_Likelihood']/samples_ps[indx]['k_vals']*np.pi   }


ps_data = {0:ps_data}
plot_power_spectrum_grid( ps_data_dir, output_dir, ps_data=ps_data, scales='large_middle', sim_data_sets=None, system='Shamrock', black_background=True  )
plot_power_spectrum_grid( ps_data_dir, output_dir, ps_data=ps_data, scales='large_reduced', sim_data_sets=None, system='Shamrock', black_background=True  )
plot_power_spectrum_grid( ps_data_dir, output_dir, ps_data=ps_data, scales='small_reduced', sim_data_sets=None, system='Shamrock', black_background=True  )




