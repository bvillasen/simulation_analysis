import sys, time, os
import numpy as np
import matplotlib.pyplot as plt
from subprocess import call
from shutil import copyfile


input_dir = '/home/bruno/Desktop/ssd_0/data/cosmo_sims/rescaled_P19/figures/phase_diagram_2048_100Mpc/fit_new/'

output_dir = '/home/bruno/Desktop/'

image_name = 'pd_fit'

out_anim_name = 'phase_diagram_rescaled_P19_fit'

start_frame = 3
frame_rate = 1


cmd = f'ffmpeg -framerate {frame_rate} -start_number {start_frame}  '
# cmd += ' -start_number 20'
cmd += ' -i {0}{1}_%d.png '.format( input_dir, image_name )
cmd += ' -pix_fmt yuv420p '
cmd += ' -vcodec libx264 '
# cmd += '-b 9100k '
cmd += '{0}{1}.mp4'.format( output_dir, out_anim_name )
# cmd += ' -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2"'
# cmd += ' -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2"'
cmd += ' -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"'
os.system( cmd )


