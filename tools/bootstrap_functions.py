import numpy as np
from tools import print_line_flush

def get_sample_mean( n_in_sample, data ):
  n_dim = data.ndim
  n_total = data.shape[0]
  sample_ids = np.random.randint( 0, n_total, n_in_sample  )
  sample = data[sample_ids]
  sample_mean = sample.mean(axis=0)
  return sample_mean


def bootstrap_sample_mean( n_iterations, n_in_sample, data, print_out ):
  sample_mean_all = []
  for i in range( n_iterations ):
    if i%(n_iterations//10)==0: print_line_flush( f' Sampling Iteration: {i}/{n_iterations}   {i/n_iterations*100}%' )
    sample_mean = get_sample_mean( n_in_sample, data )
    sample_mean_all.append( sample_mean )
  if (i+1)%(n_iterations//10)==0: print_line_flush( f' Sampling Iteration: {i+1}/{n_iterations}   {(i+1)/n_iterations*100}%' )  
  sample_mean_all = np.array( sample_mean_all )
  return sample_mean_all