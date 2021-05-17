import subprocess


fields_particles_list = [ 'pos_y', 'pos_z', 'vel_x', 'vel_y', 'vel_z' ]

fields_hydro_list = [ 'density', 'velocity_x', 'velocity_y', 'velocity_z', 'thermal_energy'  ]


type = 'particles'
# type = 'hydro'
L_MPC = 50

merge = False

if type == 'particles': fields_list = fields_particles_list
if type == 'hydro': fields_list = fields_hydro_list
if merge: fields_list = ['merge']

for field in fields_list:
  command = f'python generate_ics_from_enzo_distributed.py type={type} field={field} L_MPC={L_MPC}'
  print(f'Command: {command}')
  process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
  for line in process.stdout:
    print(line)
  process.wait()
  print(process.returncode)
 
