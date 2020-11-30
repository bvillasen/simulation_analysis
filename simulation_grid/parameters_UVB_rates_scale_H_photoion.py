
# Parameters for changing the H and He Photoionization and Photoheating Rates

param_UVB_Rates = {}

param_UVB_Rates[0] = {}
param_UVB_Rates[0]['key'] = 'A'
param_UVB_Rates[0]['name'] = 'scale_H_photoion'
param_UVB_Rates[0]['values'] = [ 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95 ]

param_UVB_Rates[1] = {}
param_UVB_Rates[1]['key'] = 'B'
param_UVB_Rates[1]['name'] = 'scale_H_photoheat'
param_UVB_Rates[1]['values'] = [ 0.87 ]

param_UVB_Rates[2] = {}
param_UVB_Rates[2]['key'] = 'C'
param_UVB_Rates[2]['name'] = 'scale_He_photoion'
param_UVB_Rates[2]['values'] = [ 0.6   ]

param_UVB_Rates[3] = {}
param_UVB_Rates[3]['key'] = 'D'
param_UVB_Rates[3]['name'] = 'scale_He_photoheat'
param_UVB_Rates[3]['values'] = [ 0.6 ]
