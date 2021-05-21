import os, sys
import numpy as np


def get_phase_diagram_bins( density, temperature, bins_dens, bins_temp  ):
  density = density / density.mean()
  density = np.log10(density).flatten()
  temperature = np.log10(temperature).flatten()
  bins_dens = np.log10(bins_dens)
  bins_temp = np.log10(bins_temp)
  if density.min() < bins_dens[0]:  print("ERROR: Density out of range")
  if density.max() > bins_dens[-1]: print("ERROR: Density out of range")
  if temperature.min() < bins_temp[0]:  print("ERROR: Temperature out of range")
  if temperature.max() > bins_temp[-1]: print("ERROR: Temperature out of range")

  phase, yedges, xedges  = np.histogram2d( density, temperature, bins=[bins_dens, bins_temp] )
  phase = phase.T
  phase = phase.astype(np.float)
  phase /= phase.sum()
  xcenters = (xedges[:-1] + xedges[1:])/2
  ycenters = (yedges[:-1] + yedges[1:])/2
  return ycenters, xcenters, phase

def fit_thermal_parameters_mcmc( n_file, values_to_fit, fit_dir, save_file=True ):
  import pymc
  import pickle

  def linear_model( overdensity_line, temp_mean, temp_sigma ):
    T0_mc  = pymc.Uniform('T0', 0, 5, value=3 )
    gamma_mc    = pymc.Uniform('gamma', -1, 1, value=0  )
    @pymc.deterministic( plot=False )
    def linear_model( overdensity_line= overdensity_line, T0=T0_mc, gamma=gamma_mc,   ):
     temperature_line  = T0 + gamma*overdensity_line
     return temperature_line
    densObsrv = pymc.Normal('line', mu=linear_model, tau=1./(temp_sigma**2), value=temp_mean, observed=True)
    return locals()

  overdensity_values = values_to_fit['density']
  temp_mean_values = values_to_fit['temperature']
  temp_sigma_values = values_to_fit['temperature_sigma']

  nIter = 100000
  nBurn = nIter / 5
  nThin = 1
  
  print( f'Fitting MCMC n_file: {n_file}')

  model = linear_model( overdensity_values, temp_mean_values, temp_sigma_values )
  linear_MDL = pymc.MCMC( model )
  linear_MDL.sample( iter=nIter, burn=nBurn, thin=nThin )
  mean_T0 = linear_MDL.stats()['T0']['mean']
  sigma_T0 = linear_MDL.stats()['T0']['standard deviation']
  mean_gamma = linear_MDL.stats()['gamma']['mean']
  sigma_gamma = linear_MDL.stats()['gamma']['standard deviation'] 
  print( f'\nFit:   T0: {mean_T0:.3f} +- {sigma_T0:.3f}      gamma:{mean_gamma:.3f} +- {sigma_gamma:.3f}' )
  # plot(linear_MDL)
  if save_file:
    outFileName = fit_dir + f'fit_{n_file}.pkl'
    f = open( outFileName, "wb")
    pickle.dump( linear_MDL.stats(), f)
    print( f'Saved File: {outFileName}' )
  return linear_MDL.stats()



def get_density_temperature_values_to_fit( pd_data, delta_min=-1, delta_max=1, n_samples_line=50, fraction_enclosed=0.70 ):
  pd = pd_data['data']
  n_dens = pd_data['n_dens']
  n_temp = pd_data['n_temp']
  temp_max = pd_data['temp_max']
  temp_min = pd_data['temp_min']
  dens_max = pd_data['dens_max']
  dens_min = pd_data['dens_min']
  log_temp_max = np.log10(temp_max)
  log_temp_min = np.log10(temp_min) 
  log_dens_max = np.log10(dens_max)
  log_dens_min = np.log10(dens_min) 
  log_temp_vals = np.linspace( log_temp_min, log_temp_max, n_temp )
  log_dens_vals = np.linspace( log_dens_min, log_dens_max, n_dens )

  indx_l = np.searchsorted( log_dens_vals, delta_min )
  indx_r = np.searchsorted( log_dens_vals, delta_max )
  
  while pd[:,indx_l].sum() <= 1e-5:   indx_l += 1
  while pd[:,indx_r].sum() <= 1e-5:   indx_r -= 1
  # print( indx_l, indx_r )
    
  overdensity_indices_all = np.linspace( indx_l, indx_r, n_samples_line ).astype(np.int)
  # overdensity_values = log_dens_vals[ overdensity_indices]

  density_values = []
  temp_mean_values = []
  temp_sigma_values = []
  temp_sigma_l_values = []
  temp_sigma_r_values = []
  
  for index,overdensity_index in enumerate(overdensity_indices_all):
    temp_slice = pd[:, overdensity_index].copy()
    slice_sum = temp_slice.sum()
    if slice_sum == 0: continue
    # if index == 0:print( temp_slice.max() )
    max_val, sigma, sigma_l, sigma_r = get_max_and_sigma( fraction_enclosed, temp_slice, log_temp_vals,  method='asymmetric'   )
    dens_val = log_dens_vals[overdensity_index]
    density_values.append( dens_val )
    temp_mean_values.append( max_val )
    temp_sigma_values.append( sigma )
    temp_sigma_l_values.append( sigma_l )
    temp_sigma_r_values.append( sigma_r )
    
  density_values = np.array(density_values)
  temp_mean_values = np.array( temp_mean_values )
  temp_sigma_values = np.array( temp_sigma_values )
  temp_sigma_l_values = np.array( temp_sigma_l_values )
  temp_sigma_r_values = np.array( temp_sigma_r_values )
  
  values_to_fit = {}
  values_to_fit['density'] = density_values
  values_to_fit['temperature'] = temp_mean_values
  values_to_fit['temperature_sigma'] = temp_sigma_values
  values_to_fit['temperature_sigma_l'] = temp_sigma_l_values
  values_to_fit['temperature_sigma_p'] = temp_sigma_r_values
  return values_to_fit

def get_max_and_sigma( fraction_enclosed, data, centers_data,  method='asymmetric' ):

  # delta_overdensity = (centers_dens[1:] - centers_dens[:-1])[0]
  delta = (centers_data[1:] - centers_data[:-1])[0]


  if method == 'symmetric': max_index, id_l, id_r, sum_fraction = get_indices_enclosed_symmetric(data, fraction_enclosed)   
  if method == 'asymmetric': max_index, id_l, id_r, sum_fraction = get_indices_enclosed_asymmetric(data, fraction_enclosed)    

  max_val = centers_data[max_index]
  sigma_l = ( max_index - id_l ) * delta 
  sigma_r = ( id_r - max_index ) * delta
  sigma = 0.5 * ( sigma_r + sigma_l )
  return max_val, sigma, sigma_l, sigma_r

def get_indices_enclosed_symmetric( data, fraction_enclosed):
  max_index = np.where(data == data.max())[0][0]
 
  for i in range( 500 ):
    # sum_fraction += data[max_index + i ]
    # sum_fraction += data[max_index - i ]
    id_l = max_index - i
    id_r = max_index + i
    sum_fraction = data[id_l:id_r+1].sum()
    # print( "{0}: {1} ".format( i, sum_fraction) )
    if sum_fraction >= fraction_enclosed:
      interval_size = i
      break
  return max_index, id_l, id_r, sum_fraction
  
def get_indices_enclosed_asymmetric( data, fraction_enclosed):
  # print(data.sum())
  data = data / data.sum()
  max_index = np.where(data == data.max())[0][0]
  sum_fraction = 0
  # print ( data[max_index], max_index)

  id_l, id_r = max_index-1, max_index+1
  val_l, val_r = data[id_l], data[id_r]
  sum_fraction += val_l + val_r
  sum_fraction = data[id_l:id_r+1].sum()

  while sum_fraction < fraction_enclosed:
    val_l, val_r = data[id_l], data[id_r]
    moved_l, moved_r  = False, False
    if val_l < val_r: 
      id_r += 1
      moved_r = True
    elif val_r < val_r: 
      id_l -= 1
      moved_l = True
    else:
      id_l -= 1
      id_r += 1
      moved_l, moved_r = True, True
    if id_l < 0: 
      id_l = 0
      id_r += 1
    sum_fraction = data[id_l:id_r+1].sum()
    # print( id_l, id_r, sum_fraction )
    
  # print('Exit')
  return max_index, id_l, id_r, sum_fraction
