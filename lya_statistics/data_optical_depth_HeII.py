import numpy as np


file_nane = 'lya_statistics/data/data_optical_depth_HeII.txt'
data = np.loadtxt( file_nane )
z, tau_HeII = data.T
data_tau_HeII = { 'z':z, 'tau':tau_HeII } 