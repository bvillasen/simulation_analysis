import os, sys
from scaling_parameter_file import Generate_Parameter_File
from configurations import Get_Configuration
analysis_dir = os.path.dirname(os.getcwd()) + '/'
sys.path.append(analysis_dir + 'tools')
from tools import *

n_per_gpu = 128

n_mpi_total = 8


n_mpi_x, n_mpi_y, n_mpi_z = Get_Configuration( n_mpi_total )
nx, ny, nz = n_mpi_x * n_per_gpu, n_mpi_y * n_per_gpu, n_mpi_z * n_per_gpu

# Box size
tile_length = 25000.0   #Lenght of box for tiling 25Mpc/h
xlen, ylen, zlen = n_mpi_x*tile_length, n_mpi_y*tile_length, n_mpi_z*tile_length  

# Directories
simulation_dir = f'/gpfs/alpine/csc434/scratch/bvilasen/scaling_2021/scale_{n_per_gpu}/n_gpus_{n_mpi_total}'
input_dir  = f'/gpfs/alpine/csc434/scratch/bvilasen/scaling_2021/scale_{n_per_gpu}/ics/'
output_dir = f'/gpfs/alpine/csc434/scratch/bvilasen/scaling_2021/scale_{n_per_gpu}/output_files/'

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


create_directory( simulation_dir )

parameter_file_str = Generate_Parameter_File( params )
params_file_name = f'{simulation_dir}/params.txt'
params_file = open( params_file_name, 'w' )
params_file.write( parameter_file_str )
params_file.close()
print( f'Saved File: {params_file_name}' )

