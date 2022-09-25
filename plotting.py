import os
import re
import numpy as np
from astropy import units
import matplotlib.pyplot as plt
from matplotlib import ticker
import glob
from pypion.ReadData import ReadData


# EDIT THIS BLOCK #######################################################
# EDIT HERE: Insert the location of simulation data (silo files)
data_path = '/mnt/local/thomas/WR140_test2'
#data_path = '/home/tony/Desktop/Simulations/Wind2DTest/SimulationData/'

# EDIT HERE: Name of your generated silo file
#filename = 'ResStudy_HD_l3n0128'
filename = 'WR140_d3l6n080'

# EDIT HERE:
# Fluid Quantity/Parameters and the corresponding tolerance
# maximum and minimum
Quantity = 'Density'
min_tolerance = np.log10(1.0e-18)
max_tolerance = np.log10(2.0e-13)
Tolerance = [min_tolerance, max_tolerance]

# If 3D, specify the coordinate plane
# Options: XY, XZ, ZY
Surface = 'XY'
# END OF THE BLOCK ********************************************************



############################################################################
# Cataloging silo files
os.chdir(data_path)
file_list = glob.glob('*.silo', recursive=True)
# file_list = sorted(glob.glob('*.silo', recursive=True))
level_list = []
files = []

for file in file_list:
    level = re.search('_level(.*)_', file)
    if level == None:
        pass
    else:
        level = level.group(1)
        if not level in level_list:
            level_list.append(level)

level_list.sort()
#print('level_list =', level_list)

# Display number of levels
if len(level_list) == 1:
    print('Simulation Info: Single level')
    catalog = []
    files = sorted(glob.glob(filename + '_0000.*.silo'))
    catalog.append(files)
else:
    print(f'Simulation Info: {len(level_list)} levels')
    catalog = []
    for i in range(len(level_list)):
        files = sorted(glob.glob(filename + '_level' + level_list[i] + '_0000.*.silo'))
        catalog.append(files)


evolution = np.array(catalog).T.tolist()

info = ReadData(evolution[0])
N_grids = info.ngrid()
N_dims = 0
for j in range(len(N_grids)):
    if N_grids[j] > 1: N_dims += 1
print(f'Simulation Info: {N_dims}D System')


os.chdir('../')
ImageDir = f'SimulationPlots/{filename}/XY-3'
if os.path.exists(ImageDir):
    print("True")