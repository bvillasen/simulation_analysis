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



SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Load_Grid_Analysis_Data()

comparable_tau = Get_Comparable_Tau()
comparable_tau_grid = Get_Comparable_Tau_from_Grid( comparable_tau, SG )

param_to_fit = 1

nIter = 100000 
nBurn = nIter / 5
nThin = 1
model = mcmc_model_1D( param_to_fit, comparable_tau, comparable_tau_grid, SG )
MDL = pymc.MCMC( model )  
MDL.sample( iter=nIter, burn=nBurn, thin=nThin )
stats = MDL.stats()


scale_H_stats = stats['scale_H']
scale_H_mean = scale_H_stats['mean']
scale_H_sigma = scale_H_stats['standard deviation'] 
tau_stats = stats['mcmc_model_1D']
tau_mean_mcmc = tau_stats['mean']
tau_sigma_mcmc = tau_stats['standard deviation']


print( f'\n scale_H: {scale_H_mean}  +- {scale_H_sigma}' )


n_samples = 10000
tau_all = []
for i in range ( n_samples ):
  scale_H = np.random.normal( scale_H_mean, scale_H_sigma )
  tau_interp, z_interp = Interpolate_Observable_1D( param_to_fit, 'tau', scale_H,  SG )
  tau_all.append( tau_interp )

tau_all = np.array( tau_all ).T
tau_mean  = [ tau_vals.mean() for tau_vals in tau_all ]
tau_sigma = []
for i in range( len(tau_all) ):
  tau_sigma_val = np.sqrt( (( tau_all[i] - tau_mean[i] )**2 ).mean() )
  tau_sigma.append( tau_sigma_val )
tau_mean  = np.array( tau_mean )
tau_sigma = np.array( tau_sigma )

 
output_dir = root_dir + 'fit_parameters/'
create_directory( output_dir )

error_colors = [ 'C9', 'C1']
font_size = 15

nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))

c0 = 'C0'
alpha = 0.6
label = 'scale_H= {0:.2}  +- {1:.2}'.format( scale_H_mean, scale_H_sigma )
ax.fill_between( z_interp, tau_mean + tau_sigma, tau_mean - tau_sigma, alpha=alpha, color=c0)
ax.plot( z_interp, tau_mean, c=c0, label=label )


data_name = 'tau_data'
data_x = comparable_tau['z']
data_mean = comparable_tau['mean']
data_error = comparable_tau['sigma'] 
ax.errorbar( data_x, data_mean, yerr=data_error, fmt='none',  alpha=0.8, ecolor= error_colors[0])
ax.scatter( data_x, data_mean, label=data_name, alpha=0.8, color= error_colors[0]) 

ax.set_ylabel( r'$\tau_{eff}$', fontsize=font_size  )
ax.set_xlabel( r'$z$', fontsize=font_size )
leg = ax.legend(loc=1, frameon=False, fontsize=font_size)
ax.set_xlim( 2, 6 )
ax.set_yscale('log')
ax.set_ylim( 0.1, 8 )

figure_name = output_dir + 'mcmc_tau.png'
fig.savefig( figure_name, bbox_inches='tight', dpi=300 )
print( f'Saved Figure: {figure_name}' )





cwd = os.getcwd()
os.chdir( output_dir )

out_file_name = output_dir + 'fit_mcmc.pkl'
f = open( out_file_name, "wb")
pickle.dump( stats, f)
print ( f'Saved File: {out_file_name}' )

pymc.Matplot.plot(MDL)  

os.chdir( cwd )  
