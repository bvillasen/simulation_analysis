from PIL import Image
import PIL.ImageOps    
from tools import *

data_dir = '/home/bruno/Desktop/ssd_0/data/'
# inDir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/figures/phase_diagram_black/'
# outDir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/figures/phase_diagram_black/anim/'

in_dir = data_dir + '/cosmo_sims/256_hydro_50Mpc/phase_diagram_nyx_dpi200/'
out_dir = data_dir + '/cosmo_sims/256_hydro_50Mpc/phase_diagram_nyx_dpi200_new/'
create_directory( out_dir )

image_name = 'phase_diagram_nyx'

n_image = 0 
for n_image in range( 500 ):
  in_image_name = in_dir + "{1}_{0}.png".format(n_image, image_name)
  out_image_name = out_dir + "{1}_{0}.png".format(n_image, image_name)

  in_image = Image.open(in_image_name)
  nx, ny = in_image.size
  nx_out, ny_out = (nx//2)*2, (ny//2)*2

  out_image = in_image.resize(( nx_out, ny_out ))
  out_image.save(out_image_name)