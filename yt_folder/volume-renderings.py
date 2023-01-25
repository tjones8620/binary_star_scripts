# from plotting_binaries_yt import *
import sys
import pathlib
home = str(pathlib.Path.home())
import os
sys.path.insert(0, os.path.join(home, "code/project/scripts/"))

from misc.make_movies import make_movies


data_path = "/mnt/massive-stars/data/thomas_simulations/wr140-sims/compton_tests/compton_mhd_n256"
img_dir = "/home/visitor_ap4/code/project/scripts/images/volume-renderings/compton_mhd_n256"
start=-1.24e7
# trajectory_file = "/home/visitor_ap4/code/project/scripts/binary_system_calculations/trajectory_phasept95.txt"
# yt_plot = YTPlotFunction(data_path, img_dir, start_time=start)
# yt_plot.volume_rendering(trajectory_file=trajectory_file)
make_movies(img_dir, img_dir, "wr140_mhd_compton_n256_density.mp4")