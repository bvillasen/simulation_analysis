import sys, os, time
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
analysis_dir = os.path.dirname(os.getcwd()) + '/'
sys.path.append(analysis_dir + 'tools')
from tools import *
#Append analysis directories to path
extend_path()
from data_thermal_history import data_thermal_history_Gaikwad_2020a, data_thermal_history_Gaikwad_2020b

