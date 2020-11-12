import numpy as np


data_keating = np.array( [
[ 4.2, 0.386215, 0.401744, 0.373139 ],
[ 4.4, 0.302689, 0.317515, 0.287006 ],
[ 4.6, 0.235014, 0.25179,  0.218802 ],
[ 4.8, 0.1695,   0.186972, 0.154833 ],
[ 5.0, 0.129563, 0.140772, 0.120686 ],
[ 5.2, 0.109078, 0.115458, 0.102487 ],
[ 5.4, 0.079947, 0.084354, 0.074827 ],
[ 5.6, 0.046133, 0.052229, 0.042031 ],
[ 5.8, 0.020964, 0.023283, 0.017795 ],
[ 6.0, 0.010206, 0.013067, 0.006533 ] ]).T
z = data_keating[0]
F = data_keating[1]
tau = -np.log(data_keating[1])
tau_p = -np.log(data_keating[2])
tau_m = -np.log(data_keating[3])
tau_error = ( tau_m - tau_p )/2 
data_optical_depth_Keating_2020 = {
'name':'Keating et al. (2020)',
'z': z,
'tau': tau,
'tau_sigma': tau_error,
'tau_sigma_p': tau_p - tau,
'tau_sigma_m': tau - tau_m
}




tau_sigma_p = np.array([ 0.04, 0.08, 0.1 ])
tau_sigma_m = np.array([ 0.04, 0.09, 0.11 ])
data_optical_depth_Boera_2019 = {
'name': 'Boera et al. (2019)',
'z': np.array([ 4.2, 4.6, 5.0 ]) ,
'tau': np.array([ 1.02, 1.41, 1.69 ]) ,
'tau_sigma_p': tau_sigma_p ,
'tau_sigma_m': tau_sigma_m ,
'tau_sigma': ( tau_sigma_p + tau_sigma_m ) / 2.0
}


data_becker = np.array( [[2.15,  0.8806 , 0.0103],
[2.25,  0.8590 , 0.0098],
[2.35,  0.8304 , 0.0093],
[2.45,  0.7968 , 0.0089],
[2.55,  0.7810 , 0.0090],
[2.65,  0.7545 , 0.0088],
[2.75,  0.7371 , 0.0088],
[2.85,  0.7167 , 0.0086],
[2.95,  0.6966 , 0.0084],
[3.05,  0.6670 , 0.0082],
[3.15,  0.6385 , 0.0080],
[3.25,  0.6031 , 0.0079],
[3.35,  0.5762 , 0.0074],
[3.45,  0.5548 , 0.0071],
[3.55,  0.5325 , 0.0071],
[3.65,  0.4992 , 0.0069],
[3.75,  0.4723 , 0.0068],
[3.85,  0.4470 , 0.0072],
[3.95,  0.4255 , 0.0071],
[4.05,  0.4030 , 0.0071],
[4.15,  0.3744 , 0.0074],
[4.25,  0.3593 , 0.0075],
[4.35,  0.3441 , 0.0102],
[4.45,  0.3216 , 0.0094],
[4.55,  0.3009 , 0.0104],
[4.65,  0.2881 , 0.0117],
[4.75,  0.2419 , 0.0201],
[4.85,  0.2225 , 0.0151]]).T
z = data_becker[0]
F = data_becker[1]
tau = -np.log(data_becker[1])
tau_error = 1/F * data_becker[2] 
data_optical_depth_Becker_2013 = {
'name': 'Becker et al. (2013)',
'z': z,
'tau': tau,
'tau_sigma': tau_error,
'tau_sigma_p': tau_error,
'tau_sigma_m': tau_error
}




data_bosman = np.array([[5.0, 0.135, 0.012 ],
[5.2, 0.114, 0.006 ],
[5.4, 0.084, 0.005 ],
[5.6, 0.050, 0.005 ],
[5.8, 0.023, 0.004 ],
[6.0, 0.0072, 0.0018 ]]).T 
z = data_bosman[0]
F = data_bosman[1]
tau = -np.log(data_bosman[1])
tau_error = 1/F * data_bosman[2] 
data_optical_depth_Bosman_2018 = {
'name':'Bosman et al. (2018)',
'z': z,
'tau': tau,
'tau_sigma': tau_error,
'tau_sigma_p': tau_error,
'tau_sigma_m': tau_error
}

# 
# data_optical_depth_Gaikward_2020 = np.array([
#     [ 2.0,  0.8690, 0.0214],
#     [ 2.2,  0.8261, 0.0206],
#     [ 2.4,  0.7919, 0.0210],
#     [ 2.6,  0.7665, 0.0216],
#     [ 2.8,  0.7398, 0.0212],
#     [ 3.0,  0.7105, 0.0213],
#     [ 3.2,  0.6731, 0.0223],
#     [ 3.4,  0.5927, 0.0247],
#     [ 3.6,  0.5320, 0.0280],
#     [ 3.8,  0.4695, 0.0278] ])