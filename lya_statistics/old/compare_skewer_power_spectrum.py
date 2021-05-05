import os, sys
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from load_data import load_analysis_data
from spectra_functions import compute_optical_depth

# Box parameters
Lbox = 50000.0 #kpc/h
nPoints = 256

delta_F = np.arange( nPoints )
n = len( delta_F )
ft = 1./n * np.fft.fft( delta_F )
ft_amp2 = ft.real * ft.real + ft.imag * ft.imag

dv = 0.5 
k_vals = np.fft.fftfreq( n, d=dv )
