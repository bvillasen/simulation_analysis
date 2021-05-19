import os, sys, time
import subprocess
#Extend path to inclide local modules
root_dir = os.path.dirname(os.getcwd())
sub_directories = [x[0] for x in os.walk(root_dir)]
sys.path.extend(sub_directories)
from tools import split_indices

# fields_particles_list = [ 'mass', 'pos_x', 'pos_y', 'pos_z', 'vel_x', 'vel_y', 'vel_z' ]
fields_particles_list = [ 'vel_z'  ]

fields_hydro_list = [ 'density', 'velocity_x', 'velocity_y', 'velocity_z', 'thermal_energy'  ]


type = 'particles'
# type = 'hydro'
L_MPC = 50

merge = True

if type == 'particles': fields_list = fields_particles_list
if type == 'hydro': fields_list = fields_hydro_list
if merge: fields_list = ['merge']

use_mpi = True
if use_mpi:
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  n_procs = comm.Get_size()
else:
  rank = 0
  n_procs = 1

n_fields = len( fields_list )
indices_local = split_indices( range(n_fields), rank, n_procs )
fields_local = [ fields_list[i] for i in indices_local ]
print( f'proc_id: {rank}  fields:{fields_local}')   

  
  
for field in fields_local:
  command = f'python generate_ics_from_enzo_distributed.py type={type} field={field} L_MPC={L_MPC}'
  print(f'Command: {command}')
  process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
  for line in process.stdout:
    print(line)
  process.wait()
  print(process.returncode)

