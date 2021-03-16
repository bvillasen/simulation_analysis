import os, sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from simulation_parameters import root_dir, system
from plot_sampling import Plot_Sampling_Corner, Plot_Sampling_T0, Plot_Sampling_Power_Spectrum


ps_data_dir = 'lya_statistics/data/'
mcmc_dir = root_dir + f'fit_mcmc/'
output_dir = mcmc_dir

labels_all = { 'P(k)+T0_Boss':'eBoss',
               'P(k)+T0_Walther': 'Walther',
               'P(k)+T0_Boss_Walther':  'eBoss + Walther',
               'P(k)+T0_Walther_Boera': 'Walther + Boera',
               'P(k)+T0_Walther_Viel':  'Walther + Viel' }


# data_set_names = [ 'P(k)+T0_Boss', 'P(k)+T0_Walther', 'P(k)+T0_Boss_Walther'  ]

data_set_names = [ 'P(k)+T0_Walther_Viel', 'P(k)+T0_Walther_Boera',   ]

samples_all = {}
for index, data_set_name in enumerate(data_set_names):
  print( f'Loading Data Set: {data_set_name}' )
  mcmc_results_dir = mcmc_dir + 'fit_results_' + data_set_name + '/'
  label = labels_all[data_set_name]
  samples_ps     = Load_Pickle_Directory( mcmc_results_dir + 'samples_power_spectrum.pkl' )
  samples_param  = Load_Pickle_Directory( mcmc_results_dir + 'samples_mcmc.pkl' )
  samples_fields = Load_Pickle_Directory( mcmc_results_dir + 'samples_fields.pkl' )
  samples_all[index] = { 'parameters':samples_param, 'fields':samples_fields, 'power_spectrum':samples_ps, 'name':data_set_name, 'label':label }


labels = { 'scale_He':r'$\beta_{\mathrm{He}}$', 'scale_H':r'$\beta_{\mathrm{H}}$', 'deltaZ_He':r'$\Delta z_{\mathrm{He}}$', 'deltaZ_H':r'$\Delta z_{\mathrm{H}}$'    }

Plot_Sampling_Corner( samples_all, labels, output_dir, figure_name='corner.png', system=system, alpha=0.7, lower_mask_factor=25 )
Plot_Sampling_T0( samples_all, output_dir, figure_name='fig_T0_sampling.png', system=system )

Plot_Sampling_Power_Spectrum( samples_all, 'small',  ps_data_dir, output_dir, figure_name='flux_ps_sampling', system=system )
Plot_Sampling_Power_Spectrum( samples_all, 'large',  ps_data_dir, output_dir, figure_name='flux_ps_sampling', system=system )
Plot_Sampling_Power_Spectrum( samples_all, 'middle', ps_data_dir, output_dir, figure_name='flux_ps_sampling', system=system )



