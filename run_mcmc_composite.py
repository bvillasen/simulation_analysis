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

output_dir = root_dir + 'fit_results_composite/'
create_directory( output_dir )

load_mcmc_stats = True


SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )
SG.Load_Grid_Analysis_Data()
sim_ids = SG.sim_ids

comparable_data = Get_Comparable_Composite_T0_tau()
comparable_grid = Get_Comparable_Composite_T0_tau_from_Grid( comparable_data, SG )


field = 'T0+tau'
stats_file = output_dir + 'fit_mcmc.pkl'


if load_mcmc_stats:
  
  print( f'Loading File: {stats_file}')
  stats = pickle.load( open( stats_file, 'rb' ) )


else:
  
  nIter = 200000 
  nBurn = nIter / 5
  nThin = 1
  model, params_mcmc = mcmc_model_4D( comparable_data, comparable_grid, field, 'mean', SG )
  MDL = pymc.MCMC( model )  
  MDL.sample( iter=nIter, burn=nBurn, thin=nThin )
  stats = MDL.stats()

  cwd = os.getcwd()
  os.chdir( output_dir )

  f = open( stats_file, 'wb' )
  pickle.dump( stats, f)
  print ( f'Saved File: {out_file_name}' )
  pymc.Matplot.plot(MDL)  

  labels, samples = [], []
  for p_id in params_mcmc.keys():
    param = params_mcmc[p_id]
    labels.append( param['name'] )
    samples.append( param['sampler'].trace())
  samples = np.array( samples ).T

  import corner
  corner_fig = corner.corner(samples[:,:], labels=labels )
  corner_fig.savefig( 'corner_fig.png' )  

  os.chdir( cwd )  












