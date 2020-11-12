import sys, os, time
import numpy as np
import h5py as h5
import pymc
import matplotlib.pyplot as plt
from tools import *
from mcmc_data_functions import Interpolate_Comparable_1D



def mcmc_model_1D( param_to_fit, comparable_data, comparable_grid, SG):
  parameters = SG.parameters
  param_name = parameters[param_to_fit]['name']
  param_vals = parameters[param_to_fit]['values']
  print(f' Fitting: {param_name}  {param_vals}')
  param_min = min(param_vals)
  param_max = max(param_vals)
  param_mid = ( param_max + param_min ) / 2.
  param_mcmc  = pymc.Uniform(param_name, param_min, param_max, value=param_mid )
  @pymc.deterministic( plot=False )
  def mcmc_model_1D( param_to_fit=param_to_fit, comparable_grid=comparable_grid, SG=SG, param_value=param_mcmc   ):
    mean_interp = Interpolate_Comparable_1D( param_to_fit, param_value,  comparable_grid, SG ) 
    return mean_interp
  densObsrv = pymc.Normal('T0', mu=mcmc_model_1D, tau=1./(comparable_data['sigma']**2), value=comparable_data['mean'], observed=True)
  return locals()
 
  

