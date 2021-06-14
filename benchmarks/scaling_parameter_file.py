


def Generate_Parameter_File( params ):
  params_str = f"""nx={params['nx']}
ny={params['ny']}
nz={params['nz']}
n_proc_x={params['n_mpi_x']}
n_proc_y={params['n_mpi_y']}
n_proc_z={params['n_mpi_z']}
tout=1000
outstep=1000
gamma=1.66666667
init=Read_Grid
nfile=0
H0=67.74
Omega_M=0.3089
Omega_L=0.6911
tile_length={params['tile_length']} 
scale_outputs_file=/ccs/home/bvilasen/cholla/scale_output_files/outputs_single_output_z0.txt 
xmin=0.0
ymin=0.0
zmin=0.0
xlen={params['xlen']}
ylen={params['ylen']}
zlen={params['zlen']}
xl_bcnd=1
xu_bcnd=1
yl_bcnd=1
yu_bcnd=1
zl_bcnd=1
zu_bcnd=1
indir={params['indir']}
outdir={params['outdir']}
"""
  return params_str
