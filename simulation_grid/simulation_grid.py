import os, sys
from pathlib import Path
import numpy as np
sys.path.append('tools')
from tools import *
import subprocess
#Append analysis directories to path
extend_path()
from parameters_UVB_rates import param_UVB_Rates
from submit_job_scripts import Create_Submit_Job_Script_Lux, Create_Submit_Job_Script_Summit
from generate_grackle_uvb_file import Generate_Modified_Rates_File
from load_data import load_analysis_data
from phase_diagram_functions import fit_thermal_parameters_mcmc, get_density_temperature_values_to_fit
from simulation_parameters import system, load_reduced_files


def Combine_List_Pair( a, b ):
  output = []
  for a_i in a:
    for b_i in b:
      if type(b_i) == list:
        add_in = [a_i] + b_i
      else:
        add_in = [ a_i, b_i]
      output.append( add_in )
  return output
  

 
class Simulation_Grid:
  n_paramters   = 0
  n_simulations = 0
  parameter_names = []
  root_dir = ''
  sim_param_indx_grid = None
  parameters = None
  Grid = None
  job_parameters = None
  simulation_parameters = None
  sim_ids = None
  
  def __init__( self, job_params=None, sim_params=None, parameters=None, dir=None ):
    print("Initializing Simulation Grid:")
    root_dir = dir
    n_param = len(parameters.keys())
    param_names = [ parameters[i]['name'] for i in range(n_param) ]
    n_sims = np.prod( [ len(parameters[i]['values']) for i in range(n_param) ] )  
    self.parameter_names = param_names
    self.n_paramters = n_param
    self.n_simulations = n_sims
    self.parameters = parameters
    self.root_dir = dir
    self.snapshots_dir = dir + 'snapshot_files/'
    self.job_parameters = job_params
    self.simulation_parameters = sim_params
    print( f" n_paramters: {self.n_paramters}")
    print( f" Paramters: {self.parameter_names}")
    print( f" n_simulations: {self.n_simulations}")
    print( f" Root Dir: {self.root_dir}")
    
    create_directory( self.root_dir )
    create_directory( self.snapshots_dir )
    
    param_keys = []
    indices_list = []
    for i in range(n_param):
      param_id = n_param - 1 - i
      param_keys.append( parameters[param_id]['key'] )
      n_vals = len( parameters[param_id]['values'] )
      indices_list.append( [ x for x in range(n_vals)] )
    
    sim_param_indx_grid = indices_list[0]
    for i in range( n_param-1 ):
      sim_param_indx_grid = Combine_List_Pair( indices_list[i+1], sim_param_indx_grid )
    assert( len(sim_param_indx_grid) == n_sims ), "N_simulations doesn't match Simulation Grid"
        
    sim_grid  = {}
    for sim_id in range( n_sims ):
      sim_grid[sim_id] = {}
      sim_grid[sim_id]['param_indices'] = sim_param_indx_grid[sim_id]
      name = 'S{0:03}'.format( sim_id )
      coords = ''
      for param_id in range(n_param):
        param = parameters[param_id]
        param_key = param['key']
        name += f'_{param_key}{sim_param_indx_grid[sim_id][param_id]}'
        coords += f'_{param_key}{sim_param_indx_grid[sim_id][param_id]}'
      sim_grid[sim_id]['key'] = name
      sim_grid[sim_id]['coords'] = coords[1:] 

    
    for sim_id in range( n_sims ):
      sim_parameters = {}
      param_indices = sim_grid[sim_id]['param_indices']
      for param_id in range(n_param):
        param = parameters[param_id]
        param_name = param['name']
        param_indx = param_indices[param_id]
        param_val = parameters[param_id]['values'][param_indx]
        sim_parameters[param_name] = param_val
      sim_grid[sim_id]['parameters'] = sim_parameters
    
    coords = {}  
    for sim_id in range( n_sims ):
      coords[sim_grid[sim_id]['coords']] = sim_id
      
    
    for sim_id in range( n_sims ):
      parameter_values = []
      for p_id in range(n_param):
        p_name = parameters[p_id]['name']
        p_val = sim_grid[sim_id]['parameters'][p_name]
        parameter_values.append( p_val )
      sim_grid[sim_id]['parameter_values'] = np.array( parameter_values )
      
        
    self.Grid = sim_grid
    self.sim_ids = self.Grid.keys()
    self.coords = coords
    
  def Find_Closest_Simulation( self, params_search ):
    diff_min = np.inf
    id_max = None 
    for sim_id in self.sim_ids:
      params_sim = self.Grid[sim_id]['parameter_values']
      diff = np.abs( params_search - params_sim ).sum()
      if diff < diff_min:
        diff_min = diff
        id_max = sim_id
    return id_max
        
  def Create_Grid_Directory_Structure( self ):
    n_sims = self.n_simulations
    for sim_id in range( n_sims ):
      sim_name = self.Grid[sim_id]['key']
      dir_name = self.root_dir + sim_name
      create_directory( dir_name )
  
  def Write_UVB_Parameters( self, sim_id ):
      sim_dir = self.Get_Simulation_Directory( sim_id )
      sim_params = self.Grid[sim_id]['parameters']
      file_name = sim_dir + 'uvb_params.txt'
      file = open( file_name, 'w' )
      for key in sim_params.keys():
        string = f'{key}={sim_params[key]} \n'
        file.write( string )
      file.close()
      print( f' Saved File: {file_name}' )
      
  def Create_Submit_Job_Script( self, sim_id, save_file=True, partition='gpuq' ):
    
    root_dir = self.root_dir
    if root_dir[-1] != '/': root_dir += '/'
    simulation = self.Grid[sim_id]
    name = simulation['key']
    job_params = self.job_parameters.copy()
    job_params['name'] = name
    job_params['sim_directory'] = root_dir + name
    job_params['partition'] = partition
    if system == 'Lux': Create_Submit_Job_Script_Lux( job_params, save_file=save_file )
    if system == 'Summit': Create_Submit_Job_Script_Summit( job_params, save_file=save_file )

  def Create_All_Submit_Job_Scripts( self, save_file=True ):
    
    print("Creating Submit Job Scripts:")
    for sim_id in self.Grid.keys():
      self.Create_Submit_Job_Script( sim_id, save_file=save_file )
  
  def Get_Simulation_Directory( self, sim_id ):
    root_dir = self.root_dir
    if root_dir[-1] != '/': root_dir += '/'
    simulation = self.Grid[sim_id]
    name = simulation['key']
    sim_dir = root_dir + name + '/'
    reduced_dir = root_dir + 'reduced_files/'
    if load_reduced_files: sim_dir = reduced_dir + name + '/'
    return sim_dir

  def Create_Directories_for_Simulation( self, sim_id ):
    
    sim_dir = self.Get_Simulation_Directory( sim_id )
    sim_name = self.Grid[sim_id]['key']
    analysis_dir = sim_dir + 'analysis_files'
    snapshots_dir = self.snapshots_dir + sim_name  
        
    # directories = ['analysis_files', 'snapshot_files' ]
    directories = [ analysis_dir, snapshots_dir ]
    for dir in directories:
      create_directory( dir )
      
  def Create_Directories_for_Simulations( self ):
    for sim_id in self.Grid.keys():
      self.Create_Directories_for_Simulation( sim_id )    
  
  def Create_Simulation_Parameter_File( self, sim_id, save_file=True ):
    
    sim_dir = self.Get_Simulation_Directory( sim_id )
    sim_params = self.simulation_parameters.copy()
    sim_params['UVB_rates_file'] = sim_dir + 'UVB_rates.h5'
    sim_params['analysisdir'] = sim_dir + 'analysis_files/'
    # sim_params['outdir'] = sim_dir + 'snapshot_files/'
    root_dir = self.root_dir
    if root_dir[-1] != '/': root_dir += '/'
    simulation = self.Grid[sim_id]
    name = self.Grid[sim_id]['key']
    outdir = root_dir + 'snapshot_files/' + name + '/'
    sim_params['outdir'] = outdir
    
    if save_file:
      file_name = sim_dir + 'param.txt'
      file = open( file_name, 'w' )
      for key in sim_params.keys():
        string = f'{key}={sim_params[key]} \n'
        file.write( string )
      file.close()
      print( f' Saved File: {file_name}' )
  
  def Create_All_Parameter_Files( self, save_file=True ):
    
    print("Creating Parameter Files:")
    for sim_id in self.Grid.keys():
      self.Create_Simulation_Parameter_File( sim_id, save_file=save_file )
      self.Write_UVB_Parameters( sim_id )
  
  def Get_Simulation_Parameter_Values( self, sim_id ):
    param = self.parameters
    n_param = self.n_paramters
    simulation = self.Grid[sim_id]
    parameter_indices = simulation['param_indices']
    param_values = {}
    for param_id in range( n_param ):
      param_name = param[param_id]['name']
      param_indx = parameter_indices[param_id]
      param_values[param_name] = param[param_id]['values'][param_indx]
    return param_values
    
  def Split_Hydro_Particles_Sanposhots( self, sim_id ):
    simulation = self.Grid[sim_id]
    sim_key = simulation['key']
    snapshots_dir = self.root_dir + 'snapshot_files/'
    if not os.path.isdir( snapshots_dir ):
      print( f'Directory Doesnt Exists: {snapshots_dir}')
      return None
    snapshots_parts_dir = self.root_dir + 'snapshot_files_particles/' 
    if not os.path.isdir( snapshots_parts_dir ):
      os.mkdir( snapshots_parts_dir )
      print( f'Created Directory: {snapshots_parts_dir}')
    
    snaps_dir = snapshots_dir + f'{sim_key}/'
    parts_dir = snapshots_parts_dir + f'{sim_key}/'
    os.mkdir( snapshots_parts_dir )
    parts_files = [ f for f in listdir(inDir) if ( isfile(join(inDir, f)) and f.find('particles')>0 )  ]
    print( parts_files )
    
    
    
    
  def Create_UVB_Rates_File( self, sim_id, max_delta_z=0.1 ):
    simulation = self.Grid[sim_id]
  
    grackle_in_file_name =  'rates_uvb/CloudyData_UVB_Puchwein2019_cloudy.h5' 
    param_values =  self.Get_Simulation_Parameter_Values( sim_id )
    sim_dir = self.Get_Simulation_Directory( sim_id )
     
    out_file_name = sim_dir + 'UVB_rates.h5'
    Generate_Modified_Rates_File( grackle_in_file_name, out_file_name, param_values, max_delta_z=max_delta_z  )
          
  def Create_UVB_Rates_Files( self, max_delta_z=0.1 ):
    print("Creating UVB Rates Files:")
    for sim_id in self.Grid.keys():
      self.Create_UVB_Rates_File( sim_id, max_delta_z=max_delta_z )
      
  def Submit_Simulation_Job( self, sim_id, partition=None ):
    sim_dir = self.Get_Simulation_Directory( sim_id )
    job = self.job_parameters
    if system == 'Summit':
      self.Create_Submit_Job_Script( sim_id, save_file=True )
      cwd = os.getcwd()
      os.chdir( sim_dir )
      command = f'bsub submit_job_summit.lsf'
    
    if system == 'Lux':
      if partition == None: partition = job['partition']
      partition_key = partition
      self.Create_Submit_Job_Script( sim_id, save_file=True, partition=partition )
      print( f'Submiting job to queue: {partition}')
      cwd = os.getcwd()
      os.chdir( sim_dir )
      if partition == 'comp-astro': partition_key = 'comp'
      if partition == 'gpuq':       partition_key = 'gpu'
      exclude_comand = '' 
      for node in job['exclude']:
        exclude_comand += node + ','
      if exclude_comand != '': exclude_comand = exclude_comand[:-1]
      command = f'submit_script {partition_key} submit_job_lux {exclude_comand}'
    print( f'Changed Directory to: {sim_dir}')
    print( f' Submitting: {command}' )
    os.system( command )
    f = open("run_output.log", "a")
    f.write('Job Submitted.\n')
    f.close()
    os.chdir( cwd ) 
  
  def Fit_Simulation_Phase_Diagram( self, sim_id ):
    print( f' Fitting Simulation: {sim_id}')
    sim_dir = self.Get_Simulation_Directory( sim_id )
    input_dir = sim_dir + 'analysis_files/'
    fit_dir = input_dir + 'fit_mcmc/'
    create_directory( fit_dir )
    
    files = [f for f in listdir(input_dir) if (isfile(join(input_dir, f)) and ( f.find('_analysis') > 0) ) ]
    indices = [ '{0:03}'.format( int(file.split('_')[0]) ) for file in files ]
    indices.sort()
    n_files = len( files )
    print( f' N_Analysis_Files: {n_files}' )
    
    for n_file in indices:
      n_file = int(n_file)
      fit_file = fit_dir + f'fit_{n_file}.pkl'
      file_path = Path(fit_file)
      if file_path.is_file():
        print( f' Skiping File: {n_file} ') 
        continue
      data = load_analysis_data( n_file, input_dir )
      values_to_fit = get_density_temperature_values_to_fit( data['phase_diagram'], delta_min=-1, delta_max=1, n_samples_line=50, fraction_enclosed=0.70 )
      fit_values = fit_thermal_parameters_mcmc( n_file, values_to_fit, fit_dir )
      
  def Fit_Grid_Phase_Diagram( self ):
    print("Fitting Phase Diagram:")
    for sim_id in self.Grid.keys():
      self.Fit_Simulation_Phase_Diagram( sim_id )


  def Fit_Simulation_Phase_Diagram_MPI( self, sim_id, n_mpi=30,  n_nodes=1  ):
    print( f' Fitting Simulation: {sim_id}')
    sim_dir = self.Get_Simulation_Directory( sim_id )
    input_dir = sim_dir + 'analysis_files/'
    # fit_dir = input_dir + 'fit_mcmc/'
    # create_directory( fit_dir )
    cwd = os.getcwd()
    run_file = cwd + '/phase_diagram/fit_phase_diagram_mpi.py'
    parameters = sim_dir + 'analysis_files/'
    n_per_node = n_mpi // n_nodes + 1
    command = f'mpirun -n {n_mpi} --map-by ppr:{n_per_node}:node --oversubscribe python {run_file} {parameters}'
    print( f' Submitting: {command}' )
    os.system( command )
    
  def Fit_Grid_Phase_Diagram_MPI( self, n_mpi=30, n_nodes=1 ):
    print("Fitting Phase Diagram:")
    for sim_id in self.Grid.keys():
      self.Fit_Simulation_Phase_Diagram_MPI( sim_id, n_mpi=n_mpi, n_nodes=n_nodes )
      
  def Load_Simulation_Analysis_Data( self, sim_id, load_fit=True  ):
    str = f' Loading Simulation Analysis: {sim_id}' 
    print_line_flush( str )
    
    sim_dir = self.Get_Simulation_Directory( sim_id )
    input_dir = sim_dir + 'analysis_files/'
    files = [f for f in listdir(input_dir) if (isfile(join(input_dir, f)) and ( f.find('_analysis') > 0) ) ]
    indices = [ '{0:03}'.format( int(file.split('_')[0]) ) for file in files ]
    indices.sort()
    n_files = len( files )
    sim_data = {}
    
    sim_data['z']      = []
    sim_data['T0']     = []
    sim_data['gamma']  = []
    sim_data['F_mean'] = []
    sim_data['tau'] = []
    sim_data['tau_HeII'] = []
    sim_data['ps_mean']  = []
    sim_data['ps_kvals'] = []
    z_power_spectrum = []
    data_ps_mean = []
    data_kvals = []
    data_kmin, data_kmax = [], []
    ps_available_indices = []
    
    
    for n_file in indices:
      n_file = int(n_file)
      data = load_analysis_data( n_file, input_dir, phase_diagram=False, lya_statistics=True, load_skewer=False, load_fit=load_fit )
      z = data['cosmology']['current_z']
      if load_fit:
        T0 =    data['phase_diagram']['fit']['T0']
        gamma = data['phase_diagram']['fit']['gamma']
      F_mean = data['lya_statistics']['Flux_mean']
      tau = data['lya_statistics']['tau']
      tau_HeII = data['lya_statistics']['tau_HeII']
      k_vals  = data['lya_statistics']['power_spectrum']['k_vals']
      ps_mean = data['lya_statistics']['power_spectrum']['ps_mean']
      if ps_mean is not None and z < 5.5:
        ps_available_indices.append(n_file)
        # z_power_spectrum.append( z )
        # data_ps_mean.append( ps_mean )
        # data_kvals.append( k_vals )
        # data_kmin.append( k_vals.min() )
        # data_kmax.append( k_vals.max() )
      sim_data['ps_kvals'].append( k_vals )
      sim_data['ps_mean'].append( ps_mean )
      sim_data['z'].append(z)
      if load_fit:
        sim_data['T0'].append(T0)
        sim_data['gamma'].append(gamma)
      sim_data['F_mean'].append(F_mean)
      sim_data['tau'].append(tau)
      sim_data['tau_HeII'].append(tau_HeII)
    sim_data['z'] = np.array( sim_data['z'] )
    if load_fit:
      sim_data['T0'] = np.array( sim_data['T0'] )
      sim_data['gamma'] = np.array( sim_data['gamma'] )
    sim_data['F_mean'] = np.array( sim_data['F_mean'] )
    sim_data['tau'] = np.array( sim_data['tau'] )
    sim_data['tau_HeII'] = np.array( sim_data['tau_HeII'] )  
    sim_data['ps_available_indices'] = ps_available_indices
    # z_power_spectrum = np.array( z_power_spectrum )
    # data_kmin, data_kmax = np.array( data_kmin ), np.array( data_kmax )
    # data_ps = { 'z':z_power_spectrum, 'k_min':data_kmin, 'k_max':data_kmax, 'k_vals':data_kvals, 'ps_mean':data_ps_mean, 'available_indices':ps_available_indices }
    # sim_data['power_spectrum'] = data_ps
    self.Grid[sim_id]['analysis'] = sim_data
    
  def Load_Simulation_Power_Spectum_Data( self, sim_id, indices, load_normalized_ps=False, ps_norm=None  ):
    
    sim_dir = self.Get_Simulation_Directory( sim_id )
    input_dir = sim_dir + 'analysis_files/'
    indices.sort()
    sim_data = {}
    
    z_vals = []
    data_ps_mean = []
    data_kvals = []
    data_kmin, data_kmax = [], [] 
    
    
    for n_file in indices:
      n_file = int(n_file)
      data = load_analysis_data( n_file, input_dir, phase_diagram=False, lya_statistics=True, load_skewer=False, load_fit=False )
      z = data['cosmology']['current_z']
      k_vals  = data['lya_statistics']['power_spectrum']['k_vals']
      ps_mean = data['lya_statistics']['power_spectrum']['ps_mean']
      z_vals.append(z)
      data_kvals.append( k_vals )
      data_ps_mean.append( ps_mean )
      data_kmin.append( k_vals.min() )
      data_kmax.append( k_vals.max() )
    z_vals = np.array( z_vals )
    data_kmin, data_kmax = np.array( data_kmin ), np.array( data_kmax )
    data_ps = { 'z':z_vals, 'k_min':data_kmin, 'k_max':data_kmax, 'k_vals':data_kvals, 'ps_mean':data_ps_mean }
    self.Grid[sim_id]['analysis']['power_spectrum'] = data_ps
    
    if load_normalized_ps:
      root_dir = self.root_dir
      sim_name = self.Grid[sim_id]['key']
      ps_dir = root_dir + f'flux_power_spectrum_files/{sim_name}/'
      normalizations = [ 'Simulation', 'Becker' ]
      types = [ 'F_mean', 'tau_eff' ] 
      normalization = ps_norm['normalization']  
      type = ps_norm['type']  
      # type = types[0]
      # normalization = 'Simulation'
      # for normalization in normalizations:
      ps_key = f'power_spectrum_norm_{normalization}_{type}'  

      z_vals, data_ps_mean, data_kvals, data_kmax, data_kmin = [], [], [], [], []
      
      for n_file in indices:
        file_name = ps_dir + f'flux_ps_{n_file}.h5'
        file = h5.File( file_name, 'r' )
        current_z = file.attrs['current_z']
        # print( file.keys() )
        k_vals = file['k_vals'][...]
        ps_data = file[normalization]
        tau_eff = ps_data.attrs['tau_eff']
        ps_mean = ps_data[type][...]
        z_vals.append( current_z )
        data_kvals.append( k_vals )
        data_ps_mean.append( ps_mean )
        data_kmax.append( k_vals.max() )
        data_kmin.append( k_vals.min() )
      z_vals = np.array( z_vals )
      data_kmin, data_kmax = np.array( data_kmin ), np.array( data_kmax )
      data_ps = { 'z':z_vals, 'k_min':data_kmin, 'k_max':data_kmax, 'k_vals':data_kvals, 'ps_mean':data_ps_mean, 'normalization_key':ps_key }
      self.Grid[sim_id]['analysis']['power_spectrum_normalized'] = data_ps
      
      
      
  

  def Load_Grid_Analysis_Data( self, sim_ids=None, load_fit=True, load_normalized_ps=False, ps_norm=None  ):
    if sim_ids == None:  
      sim_ids = self.Grid.keys()
      indx_0 = list( sim_ids )[0]
    else: indx_0 = sim_ids[0]
    
    for sim_id in sim_ids:
      self.Load_Simulation_Analysis_Data( sim_id, load_fit=load_fit  )
       
    indices = self.Grid[indx_0]['analysis']['ps_available_indices']
    available_indices = []
    for n in indices:
      available = True
      for sim_id in sim_ids:
        if n not in self.Grid[sim_id]['analysis']['ps_available_indices']: available = False
      if available: available_indices.append( n )
    
    for sim_id in sim_ids:
      self.Load_Simulation_Power_Spectum_Data( sim_id, available_indices, load_normalized_ps=load_normalized_ps, ps_norm=ps_norm )
    
    print('\n')
  
  def Get_Power_Spectrum_Range( self, sim_id=0, kmin=None, kmax=None ):
    
    data_ps = self.Grid[sim_id]['analysis']['power_spectrum']
    z = data_ps['z']
    k_min = data_ps['k_min']
    k_max = data_ps['k_max']
    # print( k_max )
    
    if kmin != None: k_min = [ max( k_min_i, kmin ) for k_min_i in k_min ] 
    if kmax != None: k_max = [ min( k_max_i, kmax ) for k_max_i in k_max ]
    # print( k_max )
    
    ps_range = { 'z':z, 'k_min':k_min, 'k_max':k_max }
    return ps_range
    
  
  def Load_Simulation_UVB_Rates( self, sim_id ):
    sim_dir = self.Get_Simulation_Directory( sim_id )
    file_name = sim_dir + 'UVB_rates.h5'
    print( f' Loading File: {file_name}')
    file = h5.File( file_name, 'r' )
    rates = file['UVBRates']
    rates_out = {}
    for root_key in rates.keys():
      rates_out[root_key] = {}
      data_group = rates[root_key]
      if root_key in [ 'Chemistry', 'Photoheating']:
        for key in data_group.keys():
          rates_out[root_key][key] = data_group[key][...]
      else:
        rates_out[root_key] = data_group[...]
    self.Grid[sim_id]['UVB_rates'] = rates_out
  
  def Load_Grid_UVB_Rates( self ):
    print( 'Loading UVB Rates Files')
    sim_ids = self.Grid.keys()
    for sim_id in sim_ids:
      self.Load_Simulation_UVB_Rates( sim_id )
    
  def Get_Simulation_Status( self, sim_id ):
    sim_dir = self.Get_Simulation_Directory( sim_id )
    file_name = sim_dir + 'run_output.log'
    file_path = Path( file_name)
    if file_path.is_file():
      file = open( file_name, 'r' )
      lines = file.readlines()
      if len(lines) == 0:
        status = 'failed'
      else:
        last_line = lines[-1]
        if last_line.find('Job Submitted') >= 0:
          status = 'submitted'
        elif last_line.find('Starting calculations') >= 0: 
          status = 'running'
        elif last_line.find('Simulation completed successfully.') >= 0: 
          status = 'finished'
        else: status = 'error'
    else:
      status = 'not submitted'
    return status
    
  def Get_Grid_Status( self, check_queue=True ):
    print( '\nGrid Status: ')
    if check_queue:  queue = self.Get_Queue_Staus()
    sim_ids = self.Grid.keys()
    n = len(sim_ids)
    submitted = 0
    running = 0
    error = 0
    finished = 0
    failed = 0
    for sim_id in sim_ids:
      queue_line = ''
      status = self.Get_Simulation_Status( sim_id )
      self.Grid[sim_id]['status'] = status
      if status == 'submitted': 
        submitted += 1
        sim_in_queue, q_line = self.Find_Simulation_In_Queue( sim_id, queue )
        if sim_in_queue:
          elem = q_line[0]
          while elem == " ":
            q_line = q_line[1:]
            elem = q_line[0]
          queue_line = q_line
        else:
          status = 'failed'
          failed += 1
      if status == 'running': 
        submitted += 1
        sim_in_queue, q_line = self.Find_Simulation_In_Queue( sim_id, queue )
        if sim_in_queue:
          running += 1
          elem = q_line[0]
          while elem == " ":
            q_line = q_line[1:]
            elem = q_line[0]
          queue_line = q_line
        else:
          status = 'failed'
          failed += 1
      if status == 'finished': 
        submitted += 1
        finished += 1
      if status == 'error':
        error += 1
      print( f' id: {sim_id}    status: {status}   {queue_line}')
      self.Grid[sim_id]['status'] = status
    print( f'Submitted: {submitted} / {n}' )
    print( f'Running:   {running} / {n}' )
    print( f'Finished:  {finished} / {n}' )
    print( f'Failed:    {failed} / {n}' )
    print( f'Error:     {error} / {n}' )
    
  def Cancel_Simulation_Job( self, sim_id ):
    status = self.Get_Simulation_Status( sim_id )
    queue = self.Get_Queue_Staus()
    sim_in_queue, q_line = self.Find_Simulation_In_Queue( sim_id, queue )
    if sim_in_queue:
      elem = q_line[0]
      while elem == " ":
        q_line = q_line[1:]
        elem = q_line[0]
      job_id = q_line.split( ' ')[0]
      if system == 'Lux': command = f'scancel {job_id}'
      if system == 'Summit': command = f'bkill {job_id}'
      print( command )
      os.system( command )
      
  def Cancel_Grid_Jobs( self, avoid=[] ):
    sim_ids = self.sim_ids
    for sim_id in sim_ids:
      if sim_id in avoid: continue
      self.Cancel_Simulation_Job( sim_id )
    
  def Submit_Grid_Jobs( self, n_submit=None, partition=None ):
    sim_ids = self.sim_ids
    self.Get_Grid_Status()
    n_submitted = 0
    for sim_id in sim_ids:
      status = self.Grid[sim_id]['status']
      # print( f' {sim_id}: {status}')
      if n_submit != None:
        if n_submitted >= n_submit: continue
      if status in ['failed', 'error', 'not submitted']:
        print( f'Submiting: {sim_id}')
        self.Submit_Simulation_Job( sim_id, partition=partition )
        n_submitted += 1
    print( f'Jobs Submitted: {n_submitted} ')
  
    
  def Delete_Simulation_Snapshots( self, sim_id ):
    sim_dir = self.Get_Simulation_Directory( sim_id )
    command = f'rm {sim_dir}snapshot_files/*'
    os.system( command )
    print ( f'Deleted {sim_dir}snapshot_files')
    
  
  def Delete_Simulation_Fit_Files( self, sim_id ):
    sim_dir = self.Get_Simulation_Directory( sim_id )
    command = f'rm {sim_dir}analysis_files/fit_mcmc/*'
    os.system( command )
    print ( f'Deleted {sim_dir}analysis_files/fit_mcmc/')
    
  
  def Delete_Grid_Snapshots( self ):
    sim_ids = self.Grid.keys()
    for sim_id in sim_ids:
      self.Delete_Simulation_Snapshots( sim_id )
  
  def Delete_Simulation_Output_Files( self, sim_id ):
    sim_dir = self.Get_Simulation_Directory( sim_id )
    command = f'rm {sim_dir}output.log'
    os.system( command )
    print ( f'Deleted {sim_dir}output.log')   
    command = f'rm {sim_dir}run_output.log'
    os.system( command )
    print ( f'Deleted {sim_dir}run_output.log')
   
  def Delete_Grid_Output_files( self ):
    sim_ids = self.Grid.keys()
    for sim_id in sim_ids:
      self.Delete_Simulation_Output_Files( sim_id )
      
  def Get_Queue_Staus( self ):
    if system == 'Lux': command = [ 'squeue', '--user=brvillas' ]
    if system == 'Summit': command = [ 'bjobs' ]
    if system == 'Shamrock': return ''
    queue = str( subprocess.check_output( command ) )
    queue = queue.split('\\n')
    return queue
    
  def Find_Simulation_In_Queue( self, sim_id, queue ):
    sim_key = 'S{0:03}'.format(sim_id)
    sim_in_queue = False
    queue_line = None
    for line in queue:
      indx = line.find( sim_key )
      if indx >= 0:
        sim_in_queue =  True
        queue_line = line
        break
    return sim_in_queue, queue_line
    
  def Reduce_Simulation_Analysis_Files( self, sim_id, reduced_dir ):
    root_dir = self.root_dir
    sim_name = self.Grid[sim_id]['key']
    # sim_dir = self.Get_Simulation_Directory( sim_id )
    sim_dir = root_dir + sim_name + '/'
    analysis_dir = sim_dir + 'analysis_files/'
    reduced_dir = reduced_dir + f'/{sim_name}/'
    print( f'Origin  Dir: {analysis_dir}' )
    create_directory( reduced_dir )
    reduced_dir = reduced_dir + 'analysis_files/'
    print( f'Destiny Dir: {reduced_dir}' )
    create_directory( reduced_dir )
    
    input_files = [f for f in listdir(analysis_dir) if (isfile(join(analysis_dir, f)) and ( f.find('_analysis') > 0) ) ]
    n_files = len( input_files )
    
    for i in range(n_files):
      file_name = f'{i}_analysis.h5'
      # Continue if reduced file already exists:
      file_path = Path(reduced_dir + file_name)
      if file_path.is_file():
        print( f'Skipping File: {reduced_dir  + file_name}' )
        continue
      in_file  = h5.File( analysis_dir + file_name, 'r' )
      out_file = h5.File( reduced_dir  + file_name, 'w' )
      
      for key in in_file.attrs:
        out_file.attrs[key] = in_file.attrs[key]
        
      file_keys = in_file.keys()
      for key in file_keys:
        if key == 'lya_statistics': continue
        in_group  = in_file[key]
        out_group = out_file.create_group( key )
        header = in_group.attrs
        for h_key in header.keys():
          out_group.attrs[h_key] = header[h_key]
        for d_key in in_group.keys():
          out_group.create_dataset( d_key, data=in_group[d_key][...] )
      
      key = 'lya_statistics'  
      lya_data = in_file[key]
      out_group = out_file.create_group( key )
      for h_key in lya_data.attrs.keys():
        out_group.attrs[h_key] = lya_data.attrs[h_key]
      group_keys = lya_data.keys()
      for group_key in group_keys:
        if 'skewers' in group_key:  continue
        group_data = lya_data[group_key]
        group = out_group.create_group( group_key )
        for d_key in group_data.keys():
          group.create_dataset( d_key, data=group_data[d_key][...] )
          
      
      in_file.close()
      out_file.close()
      print( f'Saved File: {reduced_dir + file_name}' )
      
        
  def Reduce_Grid_Analysis_Files( self, reduced_dir='reduced_files' ):
    root_dir = self.root_dir
    reduced_dir = root_dir + reduced_dir
    print( f'Copying into: {reduced_dir}' )
    create_directory( reduced_dir )
    
    sim_ids = self.sim_ids
    for sim_id in sim_ids:
      self.Reduce_Simulation_Analysis_Files( sim_id, reduced_dir )
    


