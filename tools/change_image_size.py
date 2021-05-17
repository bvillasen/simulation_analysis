from PIL import Image
import PIL.ImageOps    
from tools import *


in_dir  = '/home/bruno/Desktop/data/cosmo_sims/rescaled_P19/figures/phase_diagram_2048_100Mpc/'
out_dir = '/home/bruno/Desktop/data/cosmo_sims/rescaled_P19/figures/phase_diagram_2048_100Mpc_new/'
create_directory( out_dir )

image_name = 'phase_diagram'

for n_image in range( 56 ):
  in_image_name = in_dir + "{1}_{0}.png".format(n_image, image_name)
  out_image_name = out_dir + "{1}_{0}.png".format(n_image, image_name)

  in_image = Image.open(in_image_name)
  nx, ny = in_image.size
  nx_out, ny_out = (nx//2)*2, (ny//2)*2

  out_image = in_image.resize(( nx_out, ny_out ))
  out_image.save(out_image_name)