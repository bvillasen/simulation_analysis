import os, sys
from pathlib import Path
import numpy as np
sys.path.append('tools')
from tools import *
import subprocess
#Append analysis directories to path
extend_path()
from parameters_UVB_rates import param_UVB_Rates
from submit_job_scripts import Create_Submit_Job_Script_Lux
from generate_grackle_uvb_file import Generate_Modified_Rates_File
from load_data import load_analysis_data
from phase_diagram_functions import fit_thermal_parameters_mcmc, get_density_tyemperature_values_to_fit


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
    self.job_parameters = job_params
    self.simulation_parameters = sim_params
    print( f" n_paramters: {self.n_paramters}")
    print( f" Paramters: {self.parameter_names}")
    print( f" n_simulations: {self.n_simulations}")
    print( f" Root Dir: {self.root_dir}")
    
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
      for param_id in range(n_param):
        param = parameters[param_id]
        param_key = param['key']
        name += f'_{param_key}{sim_param_indx_grid[sim_id][param_id]}'
      sim_grid[sim_id]['key'] = name 

    
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
      
    self.Grid = sim_grid
    
  def Create_Grid_Directory_Structure( self ):
    n_sims = self.n_simulations
    for sim_id in range( n_sims ):
      sim_name = self.Grid[sim_id]['key']
      dir_name = self.root_dir + sim_name
      create_directory( dir_name )
      
  def Create_Submit_Job_Script( self, sim_id, save_file=True ):
    
    root_dir = self.root_dir
    if root_dir[-1] != '/': root_dir += '/'
    simulation = self.Grid[sim_id]
    name = simulation['key']
    job_params = self.job_parameters.copy()
    job_params['name'] = name
    job_params['sim_directory'] = root_dir + name
    Create_Submit_Job_Script_Lux( job_params, save_file=save_file )
  
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
    return sim_dir

  def Create_Directories_for_Simulation( self, sim_id ):
    sim_dir = self.Get_Simulation_Directory( sim_id )
    
    directories = ['analysis_files', 'snapshot_files' ]
    for dir in directories:
      create_directory( sim_dir + dir )
      
  def Create_Directories_for_Simulations( self ):
    for sim_id in self.Grid.keys():
      self.Create_Directories_for_Simulation( sim_id )    
  
  def Create_Simulation_Parameter_File( self, sim_id, save_file=True ):
    
    sim_dir = self.Get_Simulation_Directory( sim_id )
    sim_params = self.simulation_parameters.copy()
    sim_params['UVB_rates_file'] = sim_dir + 'UVB_rates.h5'
    sim_params['outdir'] = sim_dir + 'snapshot_files/'
    sim_params['analysisdir'] = sim_dir + 'analysis_files/'
    
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
    
    
  def Create_UVB_Rates_File( self, sim_id ):
    simulation = self.Grid[sim_id]
  
    grackle_in_file_name =  'rates_uvb/CloudyData_UVB_Puchwein2019_cloudy.h5' 
    param_values =  self.Get_Simulation_Parameter_Values( sim_id )
    sim_dir = self.Get_Simulation_Directory( sim_id )
    
    scale_HI    = param_values['scale_H']
    scale_HeII  = param_values['scale_He']
    deltaZ_HI   = param_values['deltaZ_H']
    deltaZ_HeII = param_values['deltaZ_He']
    out_file_name = sim_dir + 'UVB_rates.h5'
    Generate_Modified_Rates_File( grackle_in_file_name, out_file_name, scale_HI, scale_HeII, deltaZ_HI, deltaZ_HeII  )
          
  def Create_UVB_Rates_Files( self ):
    print("Creating UVB Rates Files:")
    for sim_id in self.Grid.keys():
      self.Create_UVB_Rates_File( sim_id )
      
  def Submit_Simulation_Job( self, sim_id ):
    sim_dir = self.Get_Simulation_Directory( sim_id )
    job = self.job_parameters
    partition = job['partition']
    partition_key = partition
    cwd = os.getcwd()
    os.chdir( sim_dir )
    if partition == 'comp-astro': partition_key = 'comp'
    command = f'submit_script {partition_key} submit_job_lux'
    print( f'Changed Directory to: {sim_dir}')
    print( f' Submitting: {command}' )
    # subprocess.call( command.split() )
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
      values_to_fit = get_density_tyemperature_values_to_fit( data['phase_diagram'], delta_min=-1, delta_max=1, n_samples_line=50, fraction_enclosed=0.70 )
      fit_values = fit_thermal_parameters_mcmc( n_file, values_to_fit, fit_dir )
      
  def Fit_Grid_Phase_Diagram( self ):
    print("Fitting Phase Diagram:")
    for sim_id in self.Grid.keys():
      self.Fit_Simulation_Phase_Diagram( sim_id )
      
  def Load_Simulation_Analysis_Data( self, sim_id, load_fit=False  ):
    print( f' Loading Simulation Analysis: {sim_id}' )
    
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
    
    for n_file in indices:
      n_file = int(n_file)
      data = load_analysis_data( n_file, input_dir, phase_diagram=False, lya_statistics=True, load_skewer=False, load_fit=True )
      z = data['cosmology']['current_z']
      T0 = data['phase_diagram']['fit']['T0']
      gamma = data['phase_diagram']['fit']['gamma']
      F_mean = data['lya_statistics']['Flux_mean']
      sim_data['z'].append(z)
      sim_data['T0'].append(T0)
      sim_data['gamma'].append(gamma)
      sim_data['F_mean'].append(F_mean)
    sim_data['z'] = np.array( sim_data['z'] )
    sim_data['T0'] = np.array( sim_data['T0'] )
    sim_data['gamma'] = np.array( sim_data['gamma'] )
    sim_data['F_mean'] = np.array( sim_data['F_mean'] )
    self.Grid[sim_id]['analysis'] = sim_data

  def Load_Grid_Analysis_Data( self, sim_ids=None, load_fit=False  ):
    if sim_ids == None:  sim_ids = self.Grid.keys()
    
    for sim_id in sim_ids:
      self.Load_Simulation_Analysis_Data( sim_id, load_fit=load_fit  )
  
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
    
  def Get_Simulation_Status( self, sim_id ):
    sim_dir = self.Get_Simulation_Directory( sim_id )
    file_name = sim_dir + 'run_output.log'
    file_path = Path( file_name)
    if file_path.is_file():
      file = open( file_name, 'r' )
      lines = file.readlines()
      last_line = lines[-1]
      if last_line.find('Job Submitted') >= 0:
        status = 'submitted'
      elif last_line.find('Starting calculations') >= 0: 
        status = 'running'
      else: status = 'error'
    else:
      status = 'not submitted'
    return status
    
  def Get_Grid_Status( self ):
    print( 'Grid Status: ')
    sim_ids = self.Grid.keys()
    n = len(sim_ids)
    submitted = 0
    running = 0
    error = 0
    finisehd = 0
    for sim_id in sim_ids:
      status = self.Get_Simulation_Status( sim_id )
      self.Grid[sim_id]['status'] = status
      if status == 'submitted': submitted += 1
      if status == 'running': 
        submitted += 1
        running += 1
      if status == 'error':
        error += 1
      print( f' id: {sim_id}    status: {status}')
    print( f'Submitted: {submitted} / {n}' )
    print( f'Running:   {running} / {n}' )
    print( f'Finished:  {finisehd} / {n}' )
    print( f'Error:     {error} / {n}' )
  
    
    
    


