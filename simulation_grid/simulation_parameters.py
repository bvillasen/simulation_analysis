

system = 'Lux'
# system = 'Shamrock'

# n_points = 256
# n_points = 512
n_points = 1024


# grid_name = f'{n_points}_P19m_np3_nsim64'
# grid_name = f'{n_points}_P19m_np2_HeII'
# grid_name = f'{n_points}_P19'
# grid_name = f'{n_points}_P19m'
# grid_name = f'{n_points}_P19m_np3_nsim8'
grid_name = f'{n_points}_P19m_np3_nsim8_v2'


if system == 'Lux':
  root_dir   = f'/data/groups/comp-astro/bruno/cosmo_sims/sim_grid/{grid_name}/'
  ics_dir    = f'/data/groups/comp-astro/bruno/cosmo_sims/sim_grid/ics/'
  cholla_dir = '/home/brvillas/cholla/'    

if system == 'Shamrock':
  root_dir   = f'/raid/bruno/data/cosmo_sims/sim_grid/{grid_name}/'
  ics_dir    = f'/raid/bruno/data/cosmo_sims/sim_grid/ics/'
  cholla_dir = '/home/bruno/cholla/'    



figures_dir = root_dir + 'figures/'


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
if n_points == 512:  sim_params['lya_skewers_stride'] = 8
if n_points == 1024: sim_params['lya_skewers_stride'] = 16
sim_params['lya_Pk_d_log_k'] = 0.1
sim_params['init'] = 'Read_Grid'
sim_params['nfile'] = 1
# sim_params['nfile'] = 1
# sim_params['indir'] = ics_dir + f'{n_points}_50Mpc/ics_8_z100/'
if n_points == 512:  sim_params['indir'] = ics_dir + f'512_50Mpc/ics_8_z20/'
if n_points == 1024: sim_params['indir'] = ics_dir + f'1024_50Mpc/ics_16_z20/'
sim_params['scale_outputs_file'] = cholla_dir + 'scale_output_files/outputs_single_output_z2.txt'
sim_params['analysis_scale_outputs_file'] = cholla_dir + 'scale_output_files/outputs_cosmo_analysis_150.txt'


job_params = {}
job_params['exclude'] = ['gpu001', 'gpu002', 'gpu003', 'gpu004' ] 
job_params['partition'] = 'gpu'
# job_params['partition'] = 'comp-astro'
if n_points == 512:
  job_params['n_mpi'] = 8
  job_params['n_nodes'] = 4
if n_points == 1024:
  job_params['n_mpi'] = 16
  job_params['n_nodes'] = 8
job_params['n_tasks_per_node'] = 2
job_params['time'] = '24:00:00'
job_params['output'] = 'output'
job_params['command'] = cholla_dir + 'cholla'
job_params['command_params'] = 'param.txt'

