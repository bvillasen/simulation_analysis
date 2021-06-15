import os, sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
import h5py as h5
root_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(subDirectories)

file_name = 'CloudyData_UVB_Puchwein2019_cloudy.h5'
file = h5.File( file_name, 'r' )
rates = file['UVBRates']

z = rates['z']
rates_heating = rates['Photoheating']
heating_HI   = rates_heating['piHI']
heating_HeI  = rates_heating['piHeI']
heating_HeII = rates_heating['piHeII']
rates_ionization = rates['Chemistry']
ionization_HI   = rates_ionization['k24']
ionization_HeI  = rates_ionization['k26']
ionization_HeII = rates_ionization['k25']








