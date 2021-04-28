import os, sys
from os import listdir
from os.path import isfile, join
import h5py as h5
import numpy as np
import subprocess
import yt
#Extend path to inclide local modules
root_dir = os.path.dirname(os.getcwd())
sub_directories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(sub_directories)
from tools import *
from constants_cgs import Msun, kpc
from ics_particles import generate_ics_particles
from ics_grid import expand_data_grid_to_cholla

# Box Size
Lbox = 50000.0    #kpc
nPoints = 256
nBoxes  = 128

current_z = 100.
H0 = 67.66
Omega_b = 0.0497

h = H0 / 100.
current_a = 1./(current_z + 1)

data_dir = '/raid/bruno/data/'
# data_dir = '/data/groups/comp-astro/bruno/'
# data_dir = '/home/bruno/Desktop/ssd_0/data/'
input_dir = data_dir + f'cosmo_sims/ics/enzo/{nPoints}_hydro_50Mpc/raw/'
output_dir = data_dir + f'cosmo_sims/{nPoints}_hydro_50Mpc/ics_{nBoxes}/'
print(f'Input Dir: {input_dir}' )
print(f'Output Dir: {output_dir}' )
# create_directory( output_dir )

hydro = True
particles = True


field_name = 'GridDensity'

file_name = input_dir + field_name
print( f'Loading File: {file_name}' )
file =  h5.File( file_name, 'r' )
file_header = file.attrs
field_data = file[field_name]
field_header = field_data.attrs
# data = field_data[...]*current_a**3/Msun*kpc**3/h**2
data_raw = field_data[...]


input_dir = data_dir + f'cosmo_sims/ics/enzo/{nPoints}_hydro_50Mpc/ics/'

nSnap = 0
snapKey = '{0:03}'.format(nSnap)
inFileName = 'DD0{0}/data0{0}'.format( snapKey)
ds = yt.load( input_dir + inFileName )
data = ds.all_data()
h = ds.hubble_constant
current_z = np.float(ds.current_redshift)
current_a = 1./(current_z + 1)

data_enzo = { 'dm':{}, 'gas':{} }
data_enzo['current_a'] = current_a
data_enzo['current_z'] = current_z

# if hydro:
print( 'Loading Hydro Data')
data_grid = ds.covering_grid( level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions )
# gas_dens = data_grid[ ('gas', 'density')].in_units('msun/kpc**3').v*current_a**3/h**2
gas_dens = data_grid[ ('gas', 'density')].v


