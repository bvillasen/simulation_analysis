import sys, os, time
import numpy as np
import h5py as h5
from tools import *
from load_data import load_snapshot_data_distributed



data_dir = '/gpfs/alpine/csc434/scratch/bvilasen/'
input_dir = data_dir + 'cosmo_sims/rescaled_P19/2048_50Mpc/snapshot_files/'
output_dir = data_dir + 'cosmo_sims/rescaled_P19/2048_50Mpc/reduced_snapshots/'


n_snap = 0


type = 'hydro'
# type = 'particles'

fields_hydro = [ 'density', 'temperature' ]
fields_particles = [ 'density' ]

if type == 'hydro': base_file_name = '.h5.'
if type == 'particles': base_file_name = '_particles.h5.'



files_snapshot = [f for f in listdir(input_dir) if f.find(f'{n_snap}{base_file_name}') == 0 ]



