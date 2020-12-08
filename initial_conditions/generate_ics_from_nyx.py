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
from internal_energy import get_internal_energy

# Box Size
Lbox = 50000.0    #kpc
nPoints = 256
nBoxes  = 8

# data_dir = '/raid/bruno/data/'
# data_dir = '/data/groups/comp-astro/bruno/'
data_dir = '/home/bruno/Desktop/ssd_0/data/'
input_dir = data_dir + f'cosmo_sims/nyx/{nPoints}_hydro_50Mpc/'
output_dir = data_dir + f'cosmo_sims/nyx/{nPoints}_hydro_50Mpc/ics_nyx_8/'
print(f'Input Dir: {input_dir}' )
print(f'Output Dir: {output_dir}' )
create_directory( output_dir )


nSnap = 0
inFileName = 'plt{0:05}'.format( nSnap )
ds = yt.load( input_dir + inFileName )
data = ds.all_data()
h = ds.hubble_constant
current_z = np.float(ds.current_redshift)
current_a = 1./(current_z + 1)

data_grid = ds.covering_grid( level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions )
gas_dens = data_grid[('boxlib', 'density')] 
gas_temp = data_grid[('boxlib', 'Temp')]    
gas_vel_x = data_grid[('boxlib', 'xmom')] / gas_dens 
gas_vel_y = data_grid[('boxlib', 'ymom')] / gas_dens
gas_vel_z = data_grid[('boxlib', 'zmom')] / gas_dens
gas_dens = gas_dens / h / h / 1e9
gas_U = get_internal_energy( gas_temp ) * gas_dens 
gas_E = 0.5 * gas_dens * ( gas_vel_x**2 + gas_vel_y**2 +gas_vel_z**2 )  + gas_U
gas_mom_x = gas_vel_x * gas_dens
gas_mom_y = gas_vel_y * gas_dens
gas_mom_z = gas_vel_z * gas_dens



p_mass = data[('all', 'particle_mass')] * h
p_pos_x = data[('all', 'particle_position_x')].in_units('kpc')/current_a * h
p_pos_y = data[('all', 'particle_position_y')].in_units('kpc')/current_a * h
p_pos_z = data[('all', 'particle_position_z')].in_units('kpc')/current_a * h
p_vel_x = data[('all', 'particle_xvel')]
p_vel_y = data[('all', 'particle_yvel')]
p_vel_z = data[('all', 'particle_zvel')]
# 
data_nyx = { 'dm':{}, 'gas':{} }
data_nyx['current_a'] = current_a
data_nyx['current_z'] = current_z
# 
data_nyx['dm']['mass'] = p_mass
data_nyx['dm']['pos_x'] = p_pos_x
data_nyx['dm']['pos_y'] = p_pos_y
data_nyx['dm']['pos_z'] = p_pos_z
data_nyx['dm']['vel_x'] = p_vel_x
data_nyx['dm']['vel_y'] = p_vel_y
data_nyx['dm']['vel_z'] = p_vel_z

data_nyx['gas']['density'] = gas_dens
data_nyx['gas']['momentum_x'] = gas_mom_x
data_nyx['gas']['momentum_y'] = gas_mom_y
data_nyx['gas']['momentum_z'] = gas_mom_z
data_nyx['gas']['GasEnergy'] = gas_U
data_nyx['gas']['Energy'] = gas_E


if nBoxes == 8: proc_grid = [ 2, 2, 2]
if nBoxes == 16: proc_grid = [ 4, 2, 2]

box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ nPoints, nPoints, nPoints ]
output_base_name = '{0}_particles.h5'.format(nSnap)
generate_ics_particles(data_nyx, output_dir, output_base_name, proc_grid, box_size, grid_size)

output_base_name = '{0}.h5'.format(nSnap)
expand_data_grid_to_cholla( proc_grid, data_nyx['gas'], output_dir, output_base_name )