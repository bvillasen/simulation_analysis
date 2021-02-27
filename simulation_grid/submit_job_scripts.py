




def Create_Submit_Job_Script_Summit( job_params, save_file=True, file_name='submit_job_summit.lsf' ):
  
  job_name = job_params['name']
  partition = job_params['partition']
  n_mpi_tasks = job_params['n_mpi']
  n_nodes = job_params['n_nodes']
  n_tasks_per_node = job_params['n_tasks_per_node']
  time = job_params['time']
  output = job_params['output']
  command = job_params['command']
  command_params = job_params['command_params']
  sim_directory = job_params['sim_directory']
  
  submit_str = f"""#!/bin/bash          
#BSUB -P {summit_project}       
#BSUB -W {time}          
#BSUB -nnodes {n_nodes}               
#BSUB -J {job_name}    
#BSUB -o {job_name}.output
#BSUB -e {job_name}.error
#BSUB -alloc_flags "smt4"

module load gcc/6.4.0
module load hdf5/1.10.4
module load cuda/10.1.243

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/ccs/home/bvilasen/code/fftw/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/ccs/home/bvilasen/code/grackle/lib

export CHOLLA_HOME=/ccs/home/bvilasen/cholla
export WORK_DIR={sim_directory}
cd {sim_directory}


date
export OMP_NUM_THREADS=7
jsrun -n {n_mpi_tasks} -a 1 -c 7 -g 1 -l CPU-CPU -d packed -b packed:7 $CHOLLA_HOME/cholla $WORK_DIR/param.txt > $WORK_DIR/run_output.log |sort

"""
  
  if save_file:
    
    if sim_directory[-1] != '/': sim_directory += '/'
    file_name = sim_directory + file_name
    file = open( file_name, 'w' )
    file.write( submit_str )
    file.close()
    print(f' Saved File: {file_name}')

  return submit_str




def Create_Submit_Job_Script_Lux( job_params, save_file=True, file_name='submit_job_lux' ):
  
  job_name = job_params['name']
  partition = job_params['partition']
  n_mpi_tasks = job_params['n_mpi']
  n_nodes = job_params['n_nodes']
  n_tasks_per_node = job_params['n_tasks_per_node']
  time = job_params['time']
  output = job_params['output']
  command = job_params['command']
  command_params = job_params['command_params']
  sim_directory = job_params['sim_directory']
  
  submit_str = f"""#!/bin/bash          
#SBATCH --job-name={job_name}    
#SBATCH --partition={partition}       
#SBATCH --ntasks={n_mpi_tasks}             
#SBATCH --nodes={n_nodes}               
#SBATCH --ntasks-per-node={n_tasks_per_node}     
#SBATCH --time={time}          
#SBATCH --output={output}.log   

module load hdf5/1.10.6
module load cuda10.2/10.2

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/brvillas/code/grackle/lib
cd {sim_directory}

set OMP_NUM_THREADS=10
mpirun -N {n_mpi_tasks} --map-by ppr:{n_tasks_per_node}:node {command} {command_params} 
"""
  
  if save_file:
    
    if sim_directory[-1] != '/': sim_directory += '/'
    file_name = sim_directory + file_name
    file = open( file_name, 'w' )
    file.write( submit_str )
    file.close()
    print(f' Saved File: {file_name}')

  return submit_str



