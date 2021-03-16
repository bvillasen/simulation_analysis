import numpy as np


def compute_distribution( values, n_bins, log=False ):
  if log: values = np.log10( values )
  val_min, val_max = values.min(), values.max()
  edges = np.linspace( val_min, val_max, n_bins+1 )
  hist, edges = np.histogram( values, bins=edges )
  # print( edges )
  hist = hist.astype( np.float )
  distribution = hist / hist.sum()
  if not log: centers = 0.5*( edges[:-1] + edges[1:] )
  else: 
    edges = 10**edges
    centers = np.sqrt( edges[:-1] * edges[1:] )
  if centers[0] > centers[-1]:
    centers = centers[::-1]
    distribution = distribution[::-1]
  return distribution, centers


def get_highest_probability_interval( bin_centers, distribution, fill_sum, log=False, n_interpolate=None):
  if log: bin_centers = np.log10( bin_centers )
  # print( f'Original {bin_centers}')

  min, max = bin_centers.min(), bin_centers.max()
  if n_interpolate: 
    bin_centers_interp = np.linspace( min, max, n_interpolate )
    distribution = np.interp( bin_centers_interp, bin_centers, distribution )
    bin_centers = bin_centers_interp
    # print( f'interpolated {bin_centers}')
  distribution = distribution / distribution.sum()
  n = len( distribution )
  v_max = distribution.max()
  id_max = np.where( distribution == v_max )[0]
  sum_val = distribution.sum()
  # if len( id_max ) > 1:
  #   print('ERROR: Unable to find unique maximum in distribution')
  #   exit(-1)
  id_max  = id_max[0]
  id_l, id_r = id_max - 1, id_max + 1
  # print( id_l, id_r )
  if id_l == -1: id_l = 0
  if id_r == n:  id_r = n-1
  # print( distribution.sum() )
  while distribution[id_l:id_r].sum() < fill_sum* sum_val:
    # print( id_l, id_r, distribution[id_l:id_r].sum()*sum_val)
    if  id_r == n-1: id_l -= 1
    elif  id_l == 0: id_r += 1
    elif distribution[id_l] < distribution[id_r]: id_r += 1
    elif distribution[id_l] > distribution[id_r]: id_l -= 1
    elif distribution[id_l] == distribution[id_r]: 
      id_l -= 1
      id_r += 1
    if id_l < 0 or id_r > n-1: break 
  if log: bin_centers = 10**bin_centers
  v_l, v_r, v_max = bin_centers[id_l], bin_centers[id_r], bin_centers[id_max]
  return v_l, v_r, v_max,  distribution[id_l:id_r].sum() 
