import os, sys
import numpy as np
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from parameters_UVB_rates import param_UVB_Rates
from simulation_grid import Simulation_Grid
from simulation_parameters import *
from plot_UVB_Rates import Plot_Grid_UVB_Rates

SG = Simulation_Grid( parameters=param_UVB_Rates, sim_params=sim_params, job_params=job_params, dir=root_dir )

reduced_snaps_dir = SG.root_dir + 'reduced_snapshot_files/'
output_root_dir = SG.root_dir + 'selected_snapshot_files/'
create_directory( output_root_dir )

fields = [ 'temperature' ]

params = { 'scale_He':None, 'deltaZ_He':None, 'scale_H':0.86, 'deltaZ_H':0.0 }
print( f'Selecting: {params} ' )
selected_sims = SG.Select_Simulations( params, tolerance=5e-3 )
n_sims = len( selected_sims )
print( f'N Selected Sims: {n_sims} ' )



sim_id = selected_sims[0]
data_sim = SG.Grid[sim_id]
sim_key = data_sim['key']
input_dir = reduced_snaps_dir + f'{sim_key}/'

in_dir_files = os.listdir( input_dir)
snaps = [ file.split('.')[0] for file in in_dir_files  ]
boxes = [ file.split('.')[-1] for file in in_dir_files  ]
snaps = list( set( snaps ) )
boxes = list( set( boxes ) )
snaps = [ int(snap) for snap in snaps ]
boxes = [ int(box) for box in boxes ]
snaps.sort()
boxes.sort()

print( f' N Snapshots: {len(snaps)} ')
print( f' N Boxes: {len(boxes)} ' )

box = boxes[0]
z_vals = []
for snap in snaps: 
  file_name = input_dir + f'{snap}.h5.{box}'
  file = h5.File( file_name, 'r' )
  z = file.attrs['Current_z'][0]
  z_vals.append( z )
  file.close()
z_vals = np.array( z_vals )


z_val = 2.8 
z_diff = np.abs( z_vals - z_val )
z_indx = np.where( z_diff == z_diff.min() )[0][0]

snap = snaps[z_indx]


sim_indx = 0
sim_id = selected_sims[sim_indx]
data_sim = SG.Grid[sim_id]
sim_key = data_sim['key']
input_dir = reduced_snaps_dir + f'{sim_key}/'

output_dir = output_root_dir + f'sim_{sim_indx}/'
create_directory( output_dir )
Write_Pickle_Directory(  data_sim, output_dir + 'data_sim.pkl' )

for box in boxes:
  in_file_name = input_dir + f'{snap}.h5.{box}'
  out_file_name = output_dir + f'{snap}.h5.{box}'
  
  infile  = h5.File( in_file_name,  'r' )
  outfile = h5.File( out_file_name, 'w' )
  
  for key in infile.attrs:
    outfile.attrs[key] = infile[key]
    
  for field in fields:
    data = infile[field][...]
    outfile.create_dataset( field, data=data )
  
  infile.close()
  outfile.close()
  print( f'Saved File: {out_file_name}' )
  
  
  

     
  

  









