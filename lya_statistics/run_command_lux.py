import subprocess


n_snap = 169
command = f'python compute_los_tau.py n_snap={n_snap}'

n_procs = 40
n_procs_per_node = 40

# lux_command = f'mpirun -n {n_procs} --map-by ppr:{n_procs_per_node}:node --oversubscribe {command}'
lux_command = f'mpirunlux {n_procs} {n_procs_per_node} {command}'
print('Command: {0}'.format( lux_command ))
process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
for line in process.stdout:
  print(line)
process.wait()
print(process.returncode)