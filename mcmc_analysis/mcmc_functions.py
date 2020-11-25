import sys, os, time
import numpy as np
import h5py as h5
import pymc
import pickle
import matplotlib.pyplot as plt
from tools import *
from mcmc_data_functions import Interpolate_Comparable_1D, Interpolate_MultiDim


def mcmc_model_4D( comparable_data, comparable_grid, field, sub_field, SG):
  print( '\nRunning MCMC Sampler')
  parameters = SG.parameters
  param_ids = parameters.keys()
  params_mcmc = {}
  for param_id in param_ids:
    param_name = parameters[param_id]['name']
    param_vals = parameters[param_id]['values']
    print(f' Fitting: {param_name}  {param_vals}')
    param_min = min(param_vals)
    param_max = max(param_vals)
    param_mid = ( param_max + param_min ) / 2.
    param_mcmc = pymc.Uniform(param_name, param_min, param_max, value=param_mid )
    params_mcmc[param_id] = {}
    params_mcmc[param_id]['sampler'] = param_mcmc
    params_mcmc[param_id]['name'] = param_name
  @pymc.deterministic( plot=False )
  def mcmc_model_4D( comparable_grid=comparable_grid, SG=SG, p0=params_mcmc[0]['sampler'], p1=params_mcmc[1]['sampler'], p2=params_mcmc[2]['sampler'], p3=params_mcmc[3]['sampler']   ):
    mean_interp = Interpolate_MultiDim( p0, p1, p2, p3, comparable_grid, field, sub_field, SG ) 
    return mean_interp
  densObsrv = pymc.Normal(field, mu=mcmc_model_4D, tau=1./(comparable_data[field]['sigma']**2), value=comparable_data[field]['mean'], observed=True)
  return locals(), params_mcmc
   
  


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
 
  


