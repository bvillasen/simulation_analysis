
# Parameters for changing the H and He Photoionization and Photoheating Rates

param_UVB_Rates = {}


param_UVB_Rates[0] = {}
param_UVB_Rates[0]['key'] = 'A'
param_UVB_Rates[0]['name'] = 'scale_He'
param_UVB_Rates[0]['values'] = [ 0.1, 0.3, 0.53, 0.76, 1.0 ]

param_UVB_Rates[1] = {}
param_UVB_Rates[1]['key'] = 'B'
param_UVB_Rates[1]['name'] = 'scale_H'
param_UVB_Rates[1]['values'] = [ 0.6, 0.73, 0.86, 1.0 ]

param_UVB_Rates[2] = {}
param_UVB_Rates[2]['key'] = 'C'
param_UVB_Rates[2]['name'] = 'deltaZ_He'
param_UVB_Rates[2]['values'] = [ -0.1, 0.2, 0.5, 0.8 ]

param_UVB_Rates[3] = {}
param_UVB_Rates[3]['key'] = 'D'
param_UVB_Rates[3]['name'] = 'deltaZ_H'
param_UVB_Rates[3]['values'] = [ -0.6, -0.4, -0.2, 0.0, 0.2 ]

