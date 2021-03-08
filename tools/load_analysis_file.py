import numpy as np
import h5py as h5



sim_dir = '/raid/bruno/data/cosmo_sims/sim_grid/1024_P19m_np4_nsim256/'
input_dir = sim_dir + 'S000_A0_B0_C0_D0/analysis_files/'


file_name = input_dir + '0_analysis.h5'

in_file = h5.File( file_name, 'r' )

pd = in_file['phase_diagram']

lya = in_file['lya_statistics']