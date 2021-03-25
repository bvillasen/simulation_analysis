import numpy as np

# 
# file_nane = 'lya_statistics/data/data_optical_depth_HeII.txt'
# data = np.loadtxt( file_nane )
# z, tau_HeII = data.T
# data_tau_HeII = { 'z':z, 'tau':tau_HeII, 'name':'Worseck et al. 2019'  } 



data_worsec = np.array([ [ 2.30, 2.54, 1.27, 0.10, 0.06, 1.20, 1.51 ],
                         [ 2.54, 2.66, 1.43, 0.17, 0.10, 1.12, 1.78 ],
                         [ 2.66, 2.74, 1.95, 0.09, 0.06, 1.58, 2.54 ],
                         [ 2.74, 2.82, 2.21, 0.10, 0.15, 1.68, 3.78 ],
                         [ 2.82, 2.94, 2.53, 0.10, 0.20, 1.92, 3.74 ] ]).T
                        
z_0, z_1, tau_median, tau_sigma_p, tau_sigma_m, tau_percentile_16, tau_percentile_84 = data_worsec 
z_bin_center = 0.5 * ( z_0 + z_1 )
# tau_sigma = 0.5 * ( tau_sigma_p + tau_sigma_m )
tau_sigma = 0.5 * ( tau_percentile_84 - tau_percentile_16   )
data_tau_HeII_Worserc_2019 = {}
data_tau_HeII_Worserc_2019['name'] = 'Worseck et al. 2019'
data_tau_HeII_Worserc_2019['z'] = z_bin_center
data_tau_HeII_Worserc_2019['tau'] = tau_median
data_tau_HeII_Worserc_2019['tau_sigma'] = tau_sigma
data_tau_HeII_Worserc_2019['tau_sigma_p'] = tau_percentile_84
data_tau_HeII_Worserc_2019['tau_sigma_m'] = tau_percentile_16




