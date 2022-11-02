from pypion.ReadData import ReadData
from pypion.argparse_command import InputValues
import yt

import numpy as np
import math
from astropy import units as u
import os
import glob
import re



path = "/mnt/massive-stars/data/thomas_simulations/wr140-sims/covertex_start/orig_res/wr140-hydro-cool-n064/"
fluid_quantity = "Density"
tolerance = [-22, -27]
plane = "XY"


data_path = path
if os.path.exists(path):
    print(f"Chosen datapath: {path}")
else:
    print(f"Directory does not exist: {path}")
    exit()
    
########### Base name of SILO files ############################
filename=sorted(os.listdir(data_path))[2].replace('_level00_0000.00000000.silo','')
print(filename)

########### Define desired fluid quantity ######################
Quantity = fluid_quantity
print(f"Fluid quantity: {Quantity}")
Tolerance = tolerance
print(f"vmin: {Tolerance[0]} \nvmax: {Tolerance[1]}")

########### If 3D, specify the coordinate plane ################
surface_options = ['XY', 'XZ', 'YZ']
Surface = plane
if surface_options.count(Surface) > 0:
    print(f"Chosen surface: {Surface}")
else:
    print(f"Invalid Surface Choice: {Surface}")
    exit()


def make_snapshots(data_path, filename):
    ########## Cataloging silo files ###############################
    os.chdir(data_path)
    file_list = glob.glob('*.silo', recursive=True)
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

    ########## Categorizing data files into levels #################
    if len(level_list) == 1: 
        print('Simulation Info: Single level')
        catalog = []
        files = sorted(glob.glob(filename + '_0000.*.silo'))
        catalog.append(files)
    else:
        print(f'Simulation Info: {len(level_list)} levels')
        catalog = []
        for i in range(len(level_list)):
            files = sorted(glob.glob(f"{filename}_level{level_list[i]}_0000.*.silo"))
            catalog.append(files)
            
    # Bundle silo files of different levels of same time instant into a snapshot.
    evolution = np.array(catalog).T
    print(f"Shape of evolution array: {evolution.shape}")
    return evolution

########### Make snapshots of the simulation ###################
evolution = make_snapshots(data_path, filename)

data = ReadData(evolution[100])
N_level = data.nlevels()
N_grids = data.ngrid()
sim_time = data.sim_time()
domain_size = data.dom_size()

baseline_data = data.get_3Darray(Quantity)
fluid_parameter = baseline_data['data']

print(N_grids, N_level, domain_size, sim_time)

dims = N_grids 

arr = np.array(evolution[100])
data_den = ReadData(arr).get_3Darray("Density")['data']
data_temp = ReadData(arr).get_3Darray("Temperature")['data']
data_velx = ReadData(arr).get_3Darray("VelocityX")['data']
data_vely = ReadData(arr).get_3Darray("VelocityY")['data']
data_velz = ReadData(arr).get_3Darray("VelocityZ")['data']
data_ngmask = ReadData(arr).get_3Darray("NG_Mask")['data']

grid_data = [dict(left_edge=[0.5-0.5**(n+1)]*len(dims), right_edge=[0.5+0.5**(n+1)]*len(dims), level=n, dimensions=dims) for n in range(N_level)]


i = 0
for g in grid_data:
    g[("gas", "density")] = (data_den[i], "g/cm**3")
    g[("gas", "temperature")] = (data_temp[i], "K")
    g[("gas", "velocity_x")] = (data_velx[i], "cm/s")
    g[("gas", "velocity_y")] = (data_vely[i], "cm/s")
    g[("gas", "velocity_z")] = (data_velz[i], "cm/s")
    # g[("gas", "pressure")] = (data_pressure[i], "g/cm**3")
    i += 1

type(grid_data)

# ds = yt.load_amr_grids(grid_data, [128, 128, 128], length_unit="7e14 * cm", refine_by=1024)

# ds.index.num_grids, ds.index.grid_dimensions

yt.enable_parallelism()

ds = yt.load_amr_grids(grid_data, [128, 128, 128], length_unit="1e15 * cm", geometry=("cartesian", ("z","y","x")),sim_time=sim_time)


sc = yt.create_scene(ds, ('gas', 'density'))
cam = sc.camera
cam.resolution = (1024, 1024) # resolution of each frame
cam.switch_orientation(normal_vector=[1, 0, 0], north_vector=[0, 0, 1])

sc.annotate_domain(ds, color=[1, 1, 1, 0.005]) # draw the domain boundary [r,g,b,alpha] sc.annotate_grids(ds, alpha=0.005) # draw the grid boundaries
sc.save('/mnt/local/thomas/frame0000.png', sigma_clip=4)
nspin = 5

for i in cam.iter_rotate(np.pi, nspin): # rotate by 180 degrees over nspin frames
    sc.save('/mnt/local/thomas/frame%04d.png' % (i+1), sigma_clip=4)