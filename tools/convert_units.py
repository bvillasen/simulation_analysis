from constants_cgs import kpc, Msun


h = 0.6766

rho_cosmo = 5.5853e-1 * 3
rho = rho_cosmo * h**2
rho_cgs = rho_cosmo * Msun / kpc**3

