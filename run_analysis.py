import os, sys
import h5py as h5
import numpy as np
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from load_data import load_analysis_data

#Parse Command Parameters
args = sys.argv[1:]
n_args = len(args)

data_dir = '/raid/bruno/data/'
input_dir = data_dir + 'cosmo_sims/256_hydro_50Mpc/analysis_files/'


n_file = 0
data = load_analysis_data( n_file, input_dir )
