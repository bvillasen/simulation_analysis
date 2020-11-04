

# system = 'Lux'
system = 'Shamrock'



if system == 'Lux':
  root_dir   = '/data/groups/comp-astro/bruno/cosmo_sims/sim_grid/'
  cholla_dir = '/home/brvillas/cholla/'    

if system == 'Shamrock':
  root_dir   = '/raid/bruno/data/cosmo_sims/sim_grid/'
  cholla_dir = '/home/bruno/cholla/'    


n_points = 512
Lbox = 50000.0 # kpc

sim_params = {}
sim_params['nx'] = n_points
sim_params['ny'] = n_points
sim_params['nz'] = n_points
sim_params['tout'] = 10000
sim_params['outstep'] = 1000
sim_params['gamma'] = 1.66666667
sim_params['H0'] = 67.66
sim_params['Omega_M'] = 0.3111
sim_params['Omega_b'] = 0.0497
sim_params['Omega_L'] = 0.6889
sim_params['End_redshift'] = 2.0
sim_params['xmin'] = 0.0
sim_params['ymin'] = 0.0
sim_params['zmin'] = 0.0 
sim_params['xlen'] = Lbox
sim_params['ylen'] = Lbox
sim_params['zlen'] = Lbox
sim_params['xl_bcnd'] = 1
sim_params['xu_bcnd'] = 1
sim_params['yl_bcnd'] = 1
sim_params['yu_bcnd'] = 1
sim_params['zl_bcnd'] = 1
sim_params['zu_bcnd'] = 1
sim_params['lya_skewers_stride'] = 8
sim_params['init'] = 'Read_Grid'
sim_params['nfile'] = 1
sim_params['indir'] = root_dir + 'ics/512_50Mpc/ics_8_z20/'
sim_params['scale_outputs_file'] = cholla_dir + 'scale_output_files/outputs_single_output_z2.txt'
sim_params['analysis_scale_outputs_file'] = cholla_dir + 'scale_output_files/outputs_cosmo_analysis_150.txt'


job_params = {}
job_params['partition'] = 'comp-astro'
job_params['n_mpi'] = 8
job_params['n_nodes'] = 8
job_params['n_tasks_per_node'] = 1
job_params['time'] = '2:00:00'
job_params['output'] = 'output'
job_params['command'] = cholla_dir + 'cholla'
job_params['command_params'] = 'param.txt'