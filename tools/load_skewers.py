import os, sys
import h5py as h5
import numpy as np


def load_skewers_single_axis(  n_skewers, skewer_axis,  nSnap, input_dir, set_random_seed=True, print_out=True ):
  inFileName = input_dir + f'skewers_{skewer_axis}_{nSnap}.h5'
  inFile = h5.File( inFileName, 'r' )
  n_total = inFile.attrs['n']
  current_z = inFile.attrs['current_z']

  if set_random_seed:   np.random.seed(12345)
  skewer_ids = np.random.randint(0, n_total, n_skewers)
  if print_out: print(f" Loading {n_skewers} skewers {skewer_axis} axis")

  skewers_dens, skewers_temp, skewers_HI, skewers_vel = [], [], [], []
  for skewer_id in skewer_ids:
    skewer_data = inFile[str(skewer_id)]
    density = skewer_data['density'][...]
    HI_density = skewer_data['HI_density'][...]
    temperature = skewer_data['temperature'][...]
    velocity = skewer_data['velocity'][...]
    skewers_dens.append( density )
    skewers_HI.append( HI_density )
    skewers_temp.append( temperature )
    skewers_vel.append(velocity)

  inFile.close() 
  data_out  = {}
  data_out['current_z'] = current_z
  data_out['n_skewers'] = n_skewers
  data_out['density']     = np.array( skewers_dens )
  data_out['HI_density']  = np.array( skewers_HI )
  data_out['temperature'] = np.array( skewers_temp )
  data_out['velocity']    = np.array( skewers_vel )
  return data_out







dataDir = '/data/groups/comp-astro/bruno/'
input_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/skewers_pchw18/'

nSnap = 169

skewer_axis = 'x'
n_skewers = 100
data_skewers = load_skewers_single_axis(n_skewers, skewer_axis,  nSnap, input_dir, set_random_seed=True, print_out=True )
current_z = data_skewers['current_z']
