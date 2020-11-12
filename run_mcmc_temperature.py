import os, sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from parameters_UVB_rates import param_UVB_Rates
from simulation_grid import Simulation_Grid
from simulation_parameters import *
from mcmc_functions import *
from mcmc_data_functions import *
from data_thermal_history import *



SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Load_Grid_Analysis_Data()

comparable_T0_Gaikwad = Get_Comparable_T0_Gaikwad()
comparable_T0_grid = Get_Comparable_T0_from_Grid( comparable_T0_Gaikwad, SG )

param_to_fit = 0

nIter = 100000 
nBurn = nIter / 5
nThin = 1
model = mcmc_model_1D( param_to_fit, comparable_T0_Gaikwad, comparable_T0_grid, SG )
MDL = pymc.MCMC( model )  
MDL.sample( iter=nIter, burn=nBurn, thin=nThin )
stats = MDL.stats()


scale_He_stats = stats['scale_He']
scale_He_mean = scale_He_stats['mean']
scale_He_sigma = scale_He_stats['standard deviation']
T0_stats = stats['mcmc_model_1D']
T0_mean_mcmc = T0_stats['mean']
T0_sigma_mcmc = T0_stats['standard deviation']


print( f'\n scale_He: {scale_He_mean}  +- {scale_He_sigma}' )


n_samples = 10000
T0_all = []
for i in range ( n_samples ):
  scale_He = np.random.normal( scale_He_mean, scale_He_sigma )
  T0_interp, z_interp = Interpolate_Observable_1D( param_to_fit, 'T0', scale_He,  SG )
  T0_all.append( T0_interp )

T0_all = np.array( T0_all ).T
T0_mean  = [ T0_vals.mean() for T0_vals in T0_all ]
T0_sigma = []
for i in range( len(T0_all) ):
  T0_sigma_val = np.sqrt( (( T0_all[i] - T0_mean[i] )**2 ).mean() )
  T0_sigma.append( T0_sigma_val )
T0_mean  = np.array( T0_mean )
T0_sigma = np.array( T0_sigma )


output_dir = root_dir + 'fit_parameters/'
create_directory( output_dir )


data_sets = [ data_thermal_history_Gaikwad_2020a, data_thermal_history_Gaikwad_2020b ]
error_colors = [ 'C9', 'C1']
font_size = 15

nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))

c0 = 'C0'
alpha = 0.6
label = 'scale_He= {0:.2}  +- {1:.2}'.format( scale_He_mean, scale_He_sigma )
ax.fill_between( z_interp, T0_mean + T0_sigma, T0_mean - T0_sigma, alpha=alpha, color=c0)
ax.plot( z_interp, T0_mean, c=c0, label=label )


for i, data_set in enumerate( data_sets ):
  data_name = data_set['name']
  data_x = data_set['z']
  data_mean = data_set['T0'].astype(np.float) 
  data_error_p = data_set['T0_sigma_plus']
  data_error_m = data_set['T0_sigma_minus']
  data_error = np.array([ data_error_m, data_error_p ]).astype(np.float) 
  ax.errorbar( data_x, data_mean, yerr=data_error, fmt='none',  alpha=0.8, ecolor= error_colors[i])
  ax.scatter( data_x, data_mean, label=data_name, alpha=0.8, color= error_colors[i]) 

ax.set_ylabel( r'$T_0$', fontsize=font_size  )
ax.set_xlabel( r'$z$', fontsize=font_size )
leg = ax.legend(loc=1, frameon=False, fontsize=font_size)
ax.set_xlim( 2, 12 )
ax.set_ylim( 3000, 18000)
figure_name = output_dir + 'mcmc_phase_diagram.png'
fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
print( f'Saved Figure: {figure_name}' )
# 
# 
# 
# 
# 
# cwd = os.getcwd()
# os.chdir( output_dir )
# 
# out_file_name = output_dir + 'fit_mcmc.pkl'
# f = open( out_file_name, "wb")
# pickle.dump( stats, f)
# print ( f'Saved File: {out_file_name}' )
# 
# pymc.Matplot.plot(MDL)  
# 
# os.chdir( cwd )  
# 