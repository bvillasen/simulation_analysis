import os, sys
from scaling_parameter_file import Generate_Parameter_File
from configurations import Get_Configuration
analysis_dir = os.path.dirname(os.getcwd()) + '/'
sys.path.append(analysis_dir + 'tools')
sys.path.append(analysis_dir + 'simulation_grid')
from tools import *
from submit_job_scripts import Create_Submit_Job_Script_Summit

GPUS_PER_NODE = 6

n_per_gpu = 128
# n_per_gpu = 256


n_mpi_list = [ 8, 64, 512, 1024  ]
# n_mpi_list = [ 2048, 4096  ]
# n_mpi_list = [ 8192, 16384  ]

for n_mpi_total in n_mpi_list:
  
  print( f'\nN MPI: {n_mpi_total} ' )

  n_nodes = ( n_mpi_total - 1 ) // GPUS_PER_NODE + 1

  time = '0:05'

  n_mpi_x, n_mpi_y, n_mpi_z = Get_Configuration( n_mpi_total )
  nx, ny, nz = n_mpi_x * n_per_gpu, n_mpi_y * n_per_gpu, n_mpi_z * n_per_gpu

  # Box size
  tile_length = 25000.0   #Lenght of box for tiling 25Mpc/h
  xlen, ylen, zlen = n_mpi_x*tile_length, n_mpi_y*tile_length, n_mpi_z*tile_length  

  # Directories
  # root_dir = '/gpfs/alpine/csc434/proj-shared/cholla/scaling_2021/'
  root_dir = '/gpfs/alpine/csc434/scratch/bvilasen/scaling_2021/'
  simulation_dir = root_dir + f'scale_{n_per_gpu}/n_gpus_{n_mpi_total}'
  input_dir      = root_dir + f'scale_{n_per_gpu}/ics/'
  output_dir     = root_dir + f'scale_{n_per_gpu}/output_files/'

  print( f'MPI Domain:  total:{n_mpi_total}  nx:{n_mpi_x}  ny:{n_mpi_y}  nz:{n_mpi_z}')
  print( f'Grid Size:  nx:{nx}  ny:{ny}  nz:{nz}')
  print( f'Box Size: xlen:{xlen} ylen:{ylen} zlen:{zlen} '  )
  print( f'Simulation Dir: {simulation_dir}')

  params = { 'nx':nx, 'ny':ny, 'nz':nz, 
    'n_mpi_x':n_mpi_x, 'n_mpi_y':n_mpi_y, 'n_mpi_z':n_mpi_z, 
    'xlen':xlen, 'ylen':ylen, 'zlen':zlen,
    'tile_length': tile_length,
    'indir': input_dir, 'outdir':output_dir
  }

  job_name = f'S_{n_per_gpu}_{n_mpi_total}'

  params_job =  {
    'summit_project': 'CSC434',
    'name': job_name,
    'n_mpi': n_mpi_total,
    'n_nodes': n_nodes,
    'time': time,
    'sim_directory': simulation_dir,
    'root_dir':root_dir 
  }


  create_directory( simulation_dir )

  parameter_file_str = Generate_Parameter_File( params )
  params_file_name = f'{simulation_dir}/param.txt'
  params_file = open( params_file_name, 'w' )
  params_file.write( parameter_file_str )
  params_file.close()
  print( f'Saved File: {params_file_name}' )


  Create_Submit_Job_Script_Summit( params_job, save_file=True, file_name='submit_job_summit.lsf' )