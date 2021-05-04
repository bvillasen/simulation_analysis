import os, sys, time
import h5py as h5
import numpy as np


def get_yt_field_hydro( field, data_grid, current_a, h ):
  if field == 'density': 
    data_field = data_grid[ ('gas', 'density')].in_units('msun/kpc**3').v*current_a**3/h**2
  if field == 'velocity_x':
    data_field = data_grid[('gas','velocity_x')].in_units('km/s').v
  if field == 'velocity_y':
    data_field = data_grid[('gas','velocity_y')].in_units('km/s').v
  if field == 'velocity_z':
    data_field = data_grid[('gas','velocity_z')].in_units('km/s').v
  if field == 'thermal_energy':
    data_field = data_grid[('gas', 'thermal_energy' )].v * 1e-10  #km^2/s^2
  return data_field  
  

def generate_ics_hydro_distributed_single_field( field, proc_grid, output_dir, ds, data_grid ):
  
  h = ds.hubble_constant
  current_z = np.float(ds.current_redshift)
  current_a = 1./(current_z + 1)

  n_proc_z, n_proc_y, n_proc_x = proc_grid

  start = time.time()
  print(f'Loading field: {field} ' )
  data = get_yt_field_hydro( field, data_grid, current_a, h )
  print(f' Writing field: {field}  {data.shape}')
  nz_total, ny_total, nx_total = data.shape
  nz, ny, nx = nz_total//n_proc_z, ny_total//n_proc_y, nx_total//n_proc_x

  for pz in range( n_proc_z ):
    zStr, zEnd = pz*nz, (pz+1)*nz
    for py in range( n_proc_y ):
      yStr, yEnd = py*ny, (py+1)*ny
      for px in range( n_proc_x ):
        xStr, xEnd = px*nx, (px+1)*nx
        proc_id = pz + py*n_proc_z + px*n_proc_z*n_proc_y
        data_local = data[zStr:zEnd, yStr:yEnd, xStr:xEnd ]
        out_file_name = output_dir + f'temp_gas_{field}_{proc_id}.h5'
        out_file = h5.File( out_file_name, 'w' )
        print(f' proc_id: {proc_id}   shape: {data_local.shape}   File: {out_file_name}' )
        out_file.create_dataset( field , data=data_local.astype(np.float64) )
        out_file.close()

  end = time.time()
  print(( ' Elapsed Time: {0:.2f} min'.format((end - start)/60.) ))


def Merge_Hydro_Fileds( field_list, proc_grid, output_dir, output_base_name='h5', n_snapshot=0 ):
  n_proc_z, n_proc_y, n_proc_x = proc_grid
  n_procs = n_proc_x * n_proc_y * n_proc_z
  
  for proc_id in range( n_procs ):
    print( f'\nproc_id: {proc_id}')
    data_local = {}
    for field in field_list:
      in_file_name = output_dir + f'temp_gas_{field}_{proc_id}.h5'
      print( f' Loading File: {in_file_name}' )
      in_file = h5.File( in_file_name, 'r' )
      data = in_file[field][...]
      in_file.close()
      data_local[field] = data
    
    density = data_local['density']
    momentum_x = data_local['velocity_x'] * density
    momentum_y = data_local['velocity_y'] * density
    momentum_z = data_local['velocity_z'] * density
    GasEnergy  = data_local['thermal_energy'] * density
    Energy = 0.5 * density * ( data_local['velocity_x']**2 + data_local['velocity_y']**2 + data_local['velocity_z']**2 ) + GasEnergy
    
    out_file_name = f'{output_dir}{n_snapshot}.{output_base_name}.{proc_id}'
    out_file = h5.File( out_file_name, 'w' )
    out_file.attrs['gamma'] = 5./3
    out_file.attrs['t'] = 0
    out_file.attrs['dt'] = 1e-10
    out_file.attrs['n_step'] = 0
    out_file.create_dataset( 'density', data=density )
    out_file.create_dataset( 'momentum_x', data=momentum_x )
    out_file.create_dataset( 'momentum_y', data=momentum_y )
    out_file.create_dataset( 'momentum_z', data=momentum_z )
    out_file.create_dataset( 'Energy', data=Energy )
    out_file.create_dataset( 'GasEnergy', data=GasEnergy )
    out_file.close()
    print( f'Saved File: {out_file_name} ' )
    





def expand_data_grid_to_cholla( proc_grid, inputData, outputDir, outputBaseName ):
  nProc_z, nProc_y, nProc_x = proc_grid
  nProc = nProc_x * nProc_y * nProc_z


  gamma = 5./3
  t = 0
  dt = 1e-10
  n_step = 0

  print( '\nGenerating ICs: Grid' )
  outFiles = {}
  for pId in range( nProc ):
    outFileName = '{0}.{1}'.format(outputBaseName, pId)
    if nProc == 1: outFileName = outputBaseName
    outFiles[pId] = h5.File( outputDir + outFileName, 'w' )
    outFiles[pId].attrs['gamma'] = gamma
    outFiles[pId].attrs['t'] = t
    outFiles[pId].attrs['dt'] = dt
    outFiles[pId].attrs['n_step'] = n_step

  fields = list(inputData.keys())
  for field in fields:
    print('Writing field: {0}'.format(field))
    data = inputData[field]
    print(data.shape)
    nz_total, ny_total, nx_total = data.shape
    nz, ny, nx = nz_total//nProc_z, ny_total//nProc_y, nx_total//nProc_x

    for pz in range( nProc_z ):
      zStr, zEnd = pz*nz, (pz+1)*nz
      for py in range( nProc_y ):
        yStr, yEnd = py*ny, (py+1)*ny
        for px in range( nProc_x ):
          xStr, xEnd = px*nx, (px+1)*nx
          pId = pz + py*nProc_z + px*nProc_z*nProc_y
          print(' File: {0}'.format(pId))
          data_local = data[zStr:zEnd, yStr:yEnd, xStr:xEnd ]
          print(data_local.shape)
          outFiles[pId].create_dataset( field , data=data_local.astype(np.float64) )

  for pId in range( nProc ):
    outFiles[pId].close()

  print('Files Saved: {0}'.format(outputDir))
  return 0
