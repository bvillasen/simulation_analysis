import sys, time, os
import numpy as np
import matplotlib.pyplot as plt
from subprocess import call
from shutil import copyfile


input_dir = '/home/bruno/Desktop/ssd_0/data/cosmo_sims/256_hydro_50Mpc/phase_diagram_nyx/'

output_dir = '/home/bruno/Desktop/'

image_name = 'phase_diagram_nyx'

out_anim_name = 'phase_diagram_nyx_cholla'


cmd = 'ffmpeg -framerate 5  '
# cmd += ' -start_number 20'
cmd += ' -i {0}{1}_%d.png '.format( input_dir, image_name )
# cmd += ' -pix_fmt yuv420p '
cmd += ' -vcodec libx264 '
# cmd += '-b 9100k '
cmd += '{0}{1}.mp4'.format( output_dir, out_anim_name )
# cmd += ' -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2"'
# cmd += ' -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2"'
# cmd += ' -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"'
os.system( cmd )


# n_files = 4
#
# n_repeat = 15
# for counter in range(1,n_repeat+1):
#   index_add = counter * n_files
#   for i in range( n_files ):
#     index_new = i + index_add
#     # print i, index_new
#     src = image_name + '_{0}.png'.format( i )
#     dst = image_name + '_{0}.png'.format( index_new )
#     print src, dst
#     copyfile(inDir + src, inDir + dst)


# inDir = '/home/bruno/Desktop/anim/'
# outDir = '/home/bruno/Desktop/anim/'
# image_name = 'gpu_model'
#
#
# out_anim_name = 'gpu_model'
