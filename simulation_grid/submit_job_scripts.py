





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
  
  submit_str = f"""
#!/bin/bash          
#SBATCH --job-name={job_name}    
#SBATCH --partition={partition}       
#SBATCH --ntasks={n_mpi_tasks}             
#SBATCH --nodes={n_nodes}               
#SBATCH --ntasks-per-node={n_tasks_per_node}     
#SBATCH --time={time}          
#SBATCH --output={output}_%j.log   

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



