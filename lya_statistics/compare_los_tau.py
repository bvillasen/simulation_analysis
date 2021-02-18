import os, sys
import numpy as np
import h5py as h5
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)
from flux_power_spectrum import get_skewer_flux_power_spectrum
from load_skewers import load_skewers_multiple_axis
from spectra_functions import compute_optical_depth
from tools import *


uvb = 'pchw18'


# dataDir = '/home/bruno/Desktop/ssd_0/data/'
dataDir = '/raid/bruno/data/'
# dataDir = '/data/groups/comp-astro/bruno/'
simulation_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/'
input_dir = simulation_dir + 'skewers_{0}_HeII/'.format(uvb)

# Box parameters
Lbox = 50000.0 #kpc/h
nPoints = 2048
nx = nPoints
ny = nPoints
nz = nPoints
ncells = nx * ny * nz
box = {'Lbox':[ Lbox, Lbox, Lbox ] }

# Cosmology parameters
cosmology = {}
cosmology['H0'] = 67.66 
cosmology['Omega_M'] = 0.3111
cosmology['Omega_L'] = 0.6889

n_snap = 169


# Load skewer data
n_skewers = [ 10, 0, 0 ]
axis_list = [ 'x' ]
skewer_dataset = load_skewers_multiple_axis( axis_list, n_skewers, n_snap, input_dir, load_HeII=True, set_random_seed=False, print_out=True )
current_z = skewer_dataset['current_z']
cosmology['current_z'] = current_z
los_density = skewer_dataset['density']
los_HI_density = skewer_dataset['HI_density']
los_HeII_density = skewer_dataset['HeII_density']
los_velocity = skewer_dataset['velocity']
los_temperature = skewer_dataset['temperature']


skewer_id = 0
skewer_data = {}  
skewer_data['HI_density']  = los_HI_density[skewer_id]
skewer_data['temperature'] = los_temperature[skewer_id]
skewer_data['velocity']    = los_velocity[skewer_id]
skewer_data['HeII_density']  = los_HeII_density[skewer_id]

tau_los_data = compute_optical_depth( cosmology, box, skewer_data, space='redshift', method='error_function', chem_type='HI' )
vel_hubble_HI = tau_los_data['vel_Hubble']
tau_HI = tau_los_data['tau']

tau_los_data = compute_optical_depth( cosmology, box, skewer_data, space='redshift', method='error_function', chem_type='HeII' )
vel_hubble_HeII = tau_los_data['vel_Hubble']
tau_HeII = tau_los_data['tau']

dens_HI = skewer_data['HI_density']
dens_HeII = skewer_data['HeII_density']

dens_fracc = dens_HI / dens_HeII 
tau_fracc = tau_HI / tau_HeII
