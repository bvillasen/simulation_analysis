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
from ics_particles import generate_ics_particles
from ics_grid import expand_data_grid_to_cholla

# data_dir = '/raid/bruno/data/'
data_dir = '/data/groups/comp-astro/bruno/'
input_dir = data_dir + 'cosmo_sims/enzo/512_hydro_50Mpc/ics/'
output_dir = data_dir + 'cosmo_sims/512_hydro_50Mpc/ics_8/'
print(f'Input Dir: {input_dir}' )
print(f'Output Dir: {output_dir}' )
create_directory( output_dir )

nSnap = 0
snapKey = '{0:03}'.format(nSnap)
inFileName = 'DD0{0}/data0{0}'.format( snapKey)
ds = yt.load( input_dir + inFileName )
data = ds.all_data()
h = ds.hubble_constant
current_z = np.float(ds.current_redshift)
current_a = 1./(current_z + 1)


data_grid = ds.covering_grid( level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions )
gas_dens = data_grid[ ('gas', 'density')].in_units('msun/kpc**3').v*current_a**3/h**2
gas_vel_x = data_grid[('gas','velocity_x')].in_units('km/s').v
gas_vel_y = data_grid[('gas','velocity_y')].in_units('km/s').v
gas_vel_z = data_grid[('gas','velocity_z')].in_units('km/s').v
gas_u = data_grid[('gas', 'thermal_energy' )].v * 1e-10 * gas_dens #km^2/s^2
gas_E = 0.5 * gas_dens * ( gas_vel_x*gas_vel_x + gas_vel_y*gas_vel_y + gas_vel_z*gas_vel_z ) + gas_u




p_mass = data[('all', 'particle_mass')].in_units('msun')*h
p_pos_x = data[('all', 'particle_position_x')].in_units('kpc')/current_a*h
p_pos_y = data[('all', 'particle_position_y')].in_units('kpc')/current_a*h
p_pos_z = data[('all', 'particle_position_z')].in_units('kpc')/current_a*h
p_vel_x = data[('all', 'particle_velocity_x')].in_units('km/s')
p_vel_y = data[('all', 'particle_velocity_y')].in_units('km/s')
p_vel_z = data[('all', 'particle_velocity_z')].in_units('km/s')

data_enzo = { 'dm':{}, 'gas':{} }
data_enzo['current_a'] = current_a
data_enzo['current_z'] = current_z

data_enzo['dm']['mass'] = p_mass
data_enzo['dm']['pos_x'] = p_pos_x
data_enzo['dm']['pos_y'] = p_pos_y
data_enzo['dm']['pos_z'] = p_pos_z
data_enzo['dm']['vel_x'] = p_vel_x
data_enzo['dm']['vel_y'] = p_vel_y
data_enzo['dm']['vel_z'] = p_vel_z

data_enzo['gas']['density'] = gas_dens
data_enzo['gas']['momentum_x'] = gas_dens * gas_vel_x
data_enzo['gas']['momentum_y'] = gas_dens * gas_vel_y
data_enzo['gas']['momentum_z'] = gas_dens * gas_vel_z
data_enzo['gas']['GasEnergy'] = gas_u
data_enzo['gas']['Energy'] = gas_E


Lbox = 50000.0
proc_grid = [ 2, 2, 2]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 512, 512, 512 ]
output_base_name = '{0}_particles.h5'.format(nSnap)
generate_ics_particles(data_enzo, output_dir, output_base_name, proc_grid, box_size, grid_size)

output_base_name = '{0}.h5'.format(nSnap)
expand_data_grid_to_cholla( proc_grid, data_enzo['gas'], output_dir, output_base_name )