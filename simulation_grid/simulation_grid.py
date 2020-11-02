import os, sys
import numpy as np
sys.path.append('tools')
from tools import *
#Append analysis directories to path
extend_path()
from parameters_UVB_rates import param_UVB_Rates
from submit_job_scripts import Create_Submit_Job_Script_Lux
from generate_grackle_uvb_file import Generate_Modified_Rates_File

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
          
