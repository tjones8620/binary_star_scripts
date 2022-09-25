# Author: Sam Green, Created: 11-10-20

# Script to call all PyPion classes to plot 2d + 3d Silo data. Works for uniform and nested-grid data.

import pypion
from pypion import Plotting_Classes
from pypion import argparse_command
#-------------------------------
import matplotlib
# Using this to stop matplotlib from using a $DISPLAY environment variable.
# i.e. This now works over ssh without the need for a x-server.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#-------------------------------
# The way I plot 3 arrays onto the same grid throws up a warning.
# This just surpresses the warning.
import warnings
warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)
#-------------------------------
line = argparse_command.InputValues()
time_dicts = line.time_dicts
dimen = line.dimen
#-------------------------------

for files in time_dicts:

	arr = time_dicts[files]

	# Change var1 to modify specific settings for your plot.
	# var1 = ["parameter", max, min, "colour scale", "axis(y or z)", ngrid/2 - 1 (so for 128^3 this is 63, 256^3 it's 127, etc)]
	var1 = ["Density", -13, -18, "viridis", 'y', 63]
	# var1 = ["Temperature", 8, 3, "inferno", 'log', 'y', 127]

	fig = plt.figure()

	# To plot 2D data use the following line to plot the data:
	#a = Plotting2d(arr).plotsilo_2d(var1[0], fig, var1)
	
	# To plot 3D data use the following line to plot the data:
	a = Plotting_Classes.Plotting3d(arr).XZXYslice(var1[0], fig, var1)

	imagefile = "%s%s_%s.png" % (line.img_path, line.img_file, time_dicts[files][0][len(time_dicts[files][0]) - 13:len(time_dicts[files][0]) - 5])
	plt.savefig(imagefile, bbox_inches='tight', dpi=300)

