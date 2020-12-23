import subprocess

snaps = [ 83, 90,  96, 102,  119, 124, 130, 136, 143, 151, 159, 169, ]
snaps_boss = [  96,  102, 106, 110,  114, 119, 124, 130, 136, 143, 151, 159 ]

snapshots = list( set( snaps_boss ).union(set(snaps)))
snapshots.sort()
print(snapshots)

uvb_list = ['pchw18', 'hm12']

for uvb in uvb_list:
  for n_snap in snapshots:
    command = 'python'
    params = f'compute_los_tau.py n_snap={n_snap} uvb={uvb}'

    n_procs = 160
    n_procs_per_node = 40

    lux_command = f'mpirun -n {n_procs} --map-by ppr:{n_procs_per_node}:node --oversubscribe {command} {params} '
    # lux_command = f'mpirunlux {n_procs} {n_procs_per_node} {command} {params}'
    print('Command: {0}'.format( lux_command ))
    process = subprocess.Popen(lux_command, shell=True, stdout=subprocess.PIPE)
    for line in process.stdout:
      print(line)
    process.wait()
    print(process.returncode)
    
for uvb in uvb_list:
  command = 'python'
  params = f'compute_flux_power_spectrum.py uvb={uvb}'
  lux_command = f'{command} {params} '
  print('Command: {0}'.format( lux_command ))
  process = subprocess.Popen(lux_command, shell=True, stdout=subprocess.PIPE)
  for line in process.stdout:
    print(line)
  process.wait()
  print(process.returncode)


    
for uvb in uvb_list:
  command = 'python'
  params = f'bootstrap_flux_power_spectrum.py n_snap={n_snap} uvb={uvb}'
  lux_command = f'{command} {params} '
  print('Command: {0}'.format( lux_command ))
  process = subprocess.Popen(lux_command, shell=True, stdout=subprocess.PIPE)
  for line in process.stdout:
    print(line)
  process.wait()
  print(process.returncode)
  
    
    

    
for uvb in uvb_list:
  command = 'python'
  params = f'compute_power_spectrum_statistics.py n_snap={n_snap} uvb={uvb}'
  lux_command = f'{command} {params} '
  print('Command: {0}'.format( lux_command ))
  process = subprocess.Popen(lux_command, shell=True, stdout=subprocess.PIPE)
  for line in process.stdout:
    print(line)
  process.wait()
  print(process.returncode)
    
    
    
    