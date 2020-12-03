
# Parameters for changing the H and He Photoionization and Photoheating Rates

param_UVB_Rates = {}

param_UVB_Rates[0] = {}
param_UVB_Rates[0]['key'] = 'A'
param_UVB_Rates[0]['name'] = 'scale_He'
# param_UVB_Rates[0]['values'] = [  0.40, 0.70,  1.0  ]
param_UVB_Rates[0]['values'] = [  1.0  ]

param_UVB_Rates[1] = {}
param_UVB_Rates[1]['key'] = 'B'
param_UVB_Rates[1]['name'] = 'scale_H'
# param_UVB_Rates[1]['values'] = [ 0.7, 0.85, 1.0 ]
param_UVB_Rates[1]['values'] = [ 1.0 ]

param_UVB_Rates[2] = {}
param_UVB_Rates[2]['key'] = 'C'
param_UVB_Rates[2]['name'] = 'deltaZ_He'
# param_UVB_Rates[2]['values'] = [ -0.3, 0,  0.3   ]
param_UVB_Rates[2]['values'] = [  0,    ]

param_UVB_Rates[3] = {}
param_UVB_Rates[3]['key'] = 'D'
param_UVB_Rates[3]['name'] = 'deltaZ_H'
# param_UVB_Rates[3]['values'] = [ -0.2, -0.1, 0.0,  ]
param_UVB_Rates[3]['values'] = [ -0.4, -0.2,  0.0, 0.2, 0.4  ]